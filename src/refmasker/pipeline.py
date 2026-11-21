"""The `refmasker run` workflow: detect, combine, summarise, report."""

from __future__ import annotations

import json
import logging
import os
import platform
import shlex
import shutil
import stat
import sys
import time
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from refmasker import __version__, external
from refmasker.annotation import FeatureIndex, detect_format, load_annotation
from refmasker.detectors import Context, Result, categories
from refmasker.detectors.plasmids import plasmid_ids
from refmasker.detectors.prophages import genomad_db, genomad_db_version
from refmasker.detectors.repeats import SENSITIVITY_IDENTITIES, SENSITIVITY_LENGTHS
from refmasker.regions import (
    MergedRegion,
    Region,
    merge,
    pad_and_clip,
    sort_by_record,
    write_annotated_bed,
    write_mask_bed,
)
from refmasker.resources import limit_library_threads, resolve_threads
from refmasker.safeio import open_replacing
from refmasker.seqio import Record, checksums, read_fasta, write_fasta
from refmasker.summary import (
    SensitivityRow,
    SummaryRow,
    repeat_sensitivity,
    summarise,
    write_summary_tsv,
)

log = logging.getLogger(__name__)

# A large prophage can overlap hundreds of genes; list at most this many per region.
MAX_GENES_LISTED = 50


@dataclass
class RunData:
    """Everything the reports and run_info.json need, collected by `_run`."""

    reference: Path
    reference_md5: str
    reference_sha256: str
    records: list[Record]
    categories: list[str]
    regions: list[Region]
    merged: list[MergedRegion]
    summary: list[SummaryRow]
    sensitivity: list[SensitivityRow]
    params: dict[str, Any]
    versions: dict[str, str]
    tool_paths: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    notes: dict[str, list[str]] = field(default_factory=dict)
    plasmids: list[str] = field(default_factory=list)
    annotation: Path | None = None
    redacted: bool = False
    command: str = ""
    started: str = ""
    runtime_s: float = 0.0

    @property
    def genome_length(self) -> int:
        """Total length of all records in bp."""
        return sum(r.length for r in self.records)


def resolve_categories(
    spec: str, extra_bed: bool, mask_plasmids: bool, gene_regex: bool = False
) -> list[str]:
    """Turn '--categories default,crispr' into an ordered list of category names.

    'default' expands to the default categories and 'all' to every category except 'user'
    and 'genes', which need extra inputs. Options that imply a category (--extra-bed,
    --mask-plasmids, --mask-genes-matching) add it. The result follows the registry order,
    whatever order the user typed.
    """
    known = categories()
    chosen: list[str] = []
    for token in (t.strip().lower() for t in spec.split(",") if t.strip()):
        if token == "default":
            names = [c.name for c in known.values() if c.default]
        elif token == "all":
            names = [c.name for c in known.values() if c.name not in {"user", "genes"}]
        elif token in known:
            names = [token]
        else:
            raise ValueError(
                f"Unknown category '{token}'. Choose from: {', '.join(known)}, default, all"
            )
        chosen.extend(n for n in names if n not in chosen)  # each name once
    # Options that imply a category switch it on...
    if extra_bed and "user" not in chosen:
        chosen.append("user")
    if mask_plasmids and "plasmid" not in chosen:
        chosen.append("plasmid")
    if gene_regex and "genes" not in chosen:
        chosen.append("genes")
    # ...and a category chosen without the option it needs fails now, before any work.
    if "user" in chosen and not extra_bed:
        raise ValueError("Category 'user' needs --extra-bed")
    if "genes" in chosen and not gene_regex:
        raise ValueError("Category 'genes' needs --mask-genes-matching")
    order = list(known)
    return sorted(chosen, key=order.index)


# Files a previous `refmasker run` leaves behind. --force only clears folders that have one,
# so pointing -o at an unrelated folder (e.g. a Nextflow project, whose work/ we would delete)
# fails instead of destroying data.
RUN_MARKERS = ("run_info.json", "refmasker.log")

# Folders a run creates inside OUTDIR; --force deletes them before running again.
RUN_FOLDERS = ("work", "figures", "logs")


def check_outdir_owner(outdir: Path) -> None:
    """Refuse an existing output folder that belongs to someone else or any user can write to.

    Another user who can write into the folder could plant a symbolic link there (e.g.
    mask.bed -> one of your files, or work -> a folder they can read) and make this run
    overwrite your file or copy the reference to them.
    """
    st = outdir.stat()  # follows a symlinked OUTDIR, which is fine if the target is yours
    if not stat.S_ISDIR(st.st_mode):
        raise FileExistsError(f"Output path {outdir} exists and is not a folder")
    if st.st_uid != os.getuid():
        raise PermissionError(
            f"Output folder {outdir} belongs to another user. Use a folder you created yourself."
        )
    if st.st_mode & stat.S_IWOTH:
        raise PermissionError(
            f"Output folder {outdir} can be written to by every user, so others could tamper "
            "with the results. Use a private folder, or remove that permission (chmod o-w)."
        )


def prepare_outdir(outdir: Path, force: bool) -> None:
    """Create the output folder, refusing anything unsafe or that is not a previous refmasker run.

    A new folder gets no access for other users; group access follows the umask. An existing
    folder must pass check_outdir_owner and, when it is not empty, needs --force, must look like
    a refmasker output folder and must not contain symbolic links. With --force, the old work/,
    figures/ and logs/ folders are deleted; the other output files are replaced as the new run
    writes them.
    """
    if not outdir.exists():
        # The umask is applied to 0o770: "others" never get access, the group only if the
        # umask allows it (e.g. umask 027 gives rwxr-x---, umask 077 gives rwx------).
        outdir.mkdir(mode=0o770, parents=True)
        return
    check_outdir_owner(outdir)
    if any(outdir.iterdir()):
        if not force:
            raise FileExistsError(f"Output folder {outdir} exists and is not empty (use --force)")
        if not any((outdir / m).is_file() for m in RUN_MARKERS):
            raise FileExistsError(
                f"Output folder {outdir} is not empty and does not look like a refmasker output "
                f"folder (no {' or '.join(RUN_MARKERS)}), so --force will not overwrite it. "
                "Choose a new or empty folder."
            )
        links = sorted(p.name for p in outdir.iterdir() if p.is_symlink())
        if links:
            raise FileExistsError(
                f"Output folder {outdir} contains symbolic links ({', '.join(links[:5])}). "
                "refmasker does not write through links in its output folder: remove them or "
                "choose a new folder."
            )
        for name in RUN_FOLDERS:
            shutil.rmtree(outdir / name, ignore_errors=True)


def add_file_logging(outdir: Path) -> logging.Handler:
    """Copy every log message (INFO and above) to OUTDIR/refmasker.log, with timestamps."""
    log_path = outdir / "refmasker.log"
    # Delete the previous run's log, then create the new one exclusively (mode "x"), so the log
    # is never written through a file or link that appeared in its place.
    log_path.unlink(missing_ok=True)
    handler = logging.FileHandler(log_path, mode="x")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    logging.getLogger().addHandler(handler)
    return handler


def input_paths(opts) -> list[str]:
    """Every path the user gave on the command line, as typed."""
    paths = [opts.reference, opts.outdir, opts.annotation, opts.genomad_db]
    paths += opts.extra_bed or []
    return [str(p) for p in paths if p]


def redactor(paths: list[str]):
    """Return a function that replaces each given path (as typed or absolute) by its file name.

    Used by --redact so shared reports do not reveal folder layouts, usernames or project
    names. Longest paths are replaced first so a folder never clobbers a file inside it.
    """
    mapping: dict[str, str] = {}
    for p in paths:
        name = Path(p).name or p
        mapping[p] = name
        mapping[str(Path(p).resolve())] = name
    keys = sorted((k for k, v in mapping.items() if k != v), key=len, reverse=True)

    def redact(text: str) -> str:
        for k in keys:
            text = text.replace(k, mapping[k])
        return text

    return redact


def annotate_genes(regions: list[Region], index: FeatureIndex | None) -> list[Region]:
    """Attach the labels of the annotated genes each region overlaps (no-op without annotation)."""
    if index is None:
        return regions
    out = []
    for r in regions:
        feats = index.overlapping(r.chrom, r.start, r.end)
        labels = [f.label for f in feats]
        if len(labels) > MAX_GENES_LISTED:
            labels = labels[:MAX_GENES_LISTED] + [f"+{len(labels) - MAX_GENES_LISTED} more"]
        out.append(replace(r, genes=tuple(labels)))
    return out


def run(opts) -> RunData:
    """Run the whole workflow for one reference (the `refmasker run` command).

    Sets up the output folder and its log file, then calls `_run`. Whether the run succeeds or
    fails, work/ is deleted afterwards (unless --keep-work) and the log file handler removed.
    """
    t0 = time.time()
    started = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    outdir = Path(opts.outdir)
    prepare_outdir(outdir, opts.force)
    handler = add_file_logging(outdir)
    succeeded = False
    try:
        data = _run(opts, outdir, t0, started)
        succeeded = True
        return data
    finally:
        if not opts.keep_work:
            clean_workdir(outdir / "work", outdir / "logs", keep_logs=not succeeded)
        logging.getLogger().removeHandler(handler)
        handler.close()


def clean_workdir(workdir: Path, logs_dir: Path, keep_logs: bool) -> None:
    """Delete work/, which holds a copy of the reference and possibly GBs of tool output.

    With `keep_logs` (after a failed run), first move each tool's log (work/<tool>/<name>.log)
    to `logs_dir`, so the failure can still be investigated.
    """
    if keep_logs and workdir.is_dir() and not workdir.is_symlink():
        tool_logs = sorted(workdir.glob("*/*.log"))
        if tool_logs:
            logs_dir.mkdir(exist_ok=True)
            for f in tool_logs:
                shutil.move(f, logs_dir / f.name)
            log.info("The run failed; tool logs are kept in %s", logs_dir)
    shutil.rmtree(workdir, ignore_errors=True)


def _run(opts, outdir: Path, t0: float, started: str) -> RunData:
    """The workflow, in order:

    1. resolve threads and categories; check the tools they need are installed
    2. read the reference; warn about draft assemblies and N bases
    3. load the annotation (optional); check its sequence names match the FASTA
    4. run each category's detector on a normalised copy of the FASTA
    5. pad, label and sort the regions; write mask.bed, mask.annotated.bed and summary.tsv
    6. collect parameters and tool versions; redact paths if asked
    7. write the reports and run_info.json (run() then deletes work/ unless --keep-work)
    """
    # ---- 1. Threads, categories and tools ----
    log.info("refmasker %s", __version__)
    # The command as typed, quoted so it can be copied and re-run.
    command = "refmasker " + " ".join(shlex.quote(a) for a in sys.argv[1:])
    log.info("Command: %s", command)

    threads = resolve_threads(opts.threads)
    # Also cap numerical libraries' thread pools, so total CPU use stays within the limit.
    limit_library_threads(threads)
    log.info("Threads: %d", threads)

    cats = resolve_categories(
        opts.categories, bool(opts.extra_bed), opts.mask_plasmids, bool(opts.mask_genes_matching)
    )
    log.info("Categories: %s", ", ".join(cats))
    registry = categories()

    # Check every tool up front rather than failing halfway through a long run.
    needed = sorted({t for c in cats for t in registry[c].tools(opts)})
    missing = external.missing_tools(needed)
    if missing:
        raise external.ToolError(
            "Missing tools on PATH: "
            + ", ".join(missing)
            + ". Install them with `conda env create -f environment.yml`, "
            "or drop the categories that need them with --categories."
        )
    # Record where each tool was found: for provenance, and to warn below about any program
    # that another user could have replaced.
    tool_paths = {t: shutil.which(t) or "" for t in needed}
    for tool, path in tool_paths.items():
        log.info("Using %s: %s", tool, path)

    # ---- 2. Reference genome ----
    reference = Path(opts.reference)
    max_bp = opts.max_genome_mb * 1_000_000 if opts.max_genome_mb else None  # 0 = no limit
    records = read_fasta(reference, max_bp)
    lengths = {r.id: r.length for r in records}
    log.info(
        "Reference: %s, %d record(s), %s bp", reference, len(records), f"{sum(lengths.values()):,}"
    )

    warnings: list[str] = []
    swappable = [t for t, p in tool_paths.items() if external.writable_by_others(p)]
    if swappable:
        warnings.append(
            "These tools are in a location any user can modify, so they could have been "
            f"replaced: {', '.join(swappable)}. Use a conda environment that only you or your "
            "group can write to."
        )
    plasmids = plasmid_ids(records, opts.plasmid_ids)
    chromosomal = [r for r in records if r.id not in plasmids]
    # A complete genome has one chromosome (plus plasmids); more records suggests a draft.
    if len(chromosomal) > 1:
        warnings.append(
            f"The reference has {len(chromosomal)} non-plasmid records. If this is a draft "
            "assembly, repeats and IS elements are often collapsed or split at contig ends, so "
            "they will be under-detected. A complete genome is recommended."
        )
    n_count = sum(r.seq.upper().count("N") for r in records)
    if n_count:
        warnings.append(f"The reference contains {n_count:,} N bases.")

    # ---- 3. Annotation (optional) ----
    features = None
    index = None
    annotation = Path(opts.annotation) if opts.annotation else None
    opts.annotation_format = None
    if annotation is not None:
        opts.annotation_format = detect_format(annotation)
        features = load_annotation(annotation)
        index = FeatureIndex(features)
        log.info(
            "Annotation: %s (%s), %d features", annotation, opts.annotation_format, len(features)
        )
        # Sequence names that differ between annotation and FASTA would silently give empty
        # gene lists, so warn about them.
        unmatched = index.seqids - set(lengths)
        if not index.seqids & set(lengths):
            warnings.append(
                "None of the sequence names in the annotation match the FASTA, so no genes are "
                f"listed. Annotation has: {', '.join(sorted(index.seqids)[:5])}."
            )
        elif unmatched:
            warnings.append(
                f"Annotation sequences not in the FASTA were ignored: {', '.join(sorted(unmatched))}."
            )

    # ---- 4. Detectors ----
    # Tools run on a copy of the reference with ID-only headers, so they all report the same
    # sequence names as the FASTA IDs used in mask.bed.
    # work/ is private to the user (mode 0o700) because it holds a copy of the reference.
    workdir = outdir / "work"
    workdir.mkdir(mode=0o700)
    fasta = workdir / "reference.fna"
    write_fasta(fasta, records)
    ctx = Context(fasta, records, workdir, threads, opts, features, annotation)

    raw: list[Region] = []
    notes: dict[str, list[str]] = {}
    results: dict[str, Result] = {}
    for name in cats:
        t = time.time()
        log.info("Detecting %s ...", name)
        res = registry[name].run(ctx)
        results[name] = res
        # Regions on sequences not in the FASTA cannot go into mask.bed; drop them with a warning.
        bad = [r for r in res.regions if r.chrom not in lengths]
        if bad:
            names = sorted({r.chrom for r in bad})
            warnings.append(
                f"{name}: {len(bad)} region(s) on sequences not in the FASTA were dropped "
                f"({', '.join(names[:5])})."
            )
        raw.extend(r for r in res.regions if r.chrom in lengths)
        if res.notes:
            notes[name] = res.notes
        log.info("%s: %d regions in %.0f s", name, len(res.regions), time.time() - t)

    # ---- 5. Combine regions and write the mask files ----
    # Pad each region (--pad, 0 by default), list the genes it overlaps, sort in FASTA order.
    regions = pad_and_clip(raw, opts.pad, lengths)
    regions = annotate_genes(regions, index)
    order = [r.id for r in records]
    regions = sort_by_record(regions, order)
    merged = sort_by_record(merge(regions), order)

    # mask.bed holds the merged intervals for snippy; mask.annotated.bed every region and why.
    write_mask_bed(outdir / "mask.bed", merged)
    write_annotated_bed(outdir / "mask.annotated.bed", regions)
    summary = summarise(regions, records, cats)
    write_summary_tsv(outdir / "summary.tsv", summary)

    # How the repeat mask would change at other thresholds, from the same nucmer alignments.
    sensitivity: list[SensitivityRow] = []
    if "repeats" in results:
        sensitivity = repeat_sensitivity(
            results["repeats"].extra["hits"],
            sum(lengths.values()),
            opts.repeat_min_len,
            opts.repeat_min_id,
            SENSITIVITY_LENGTHS,
            SENSITIVITY_IDENTITIES,
        )

    # ---- 6. Provenance: tool versions and parameters ----
    versions = {"refmasker": __version__, "python": platform.python_version()}
    for tool in needed:
        versions[tool] = external.tool_version(tool)
    if "prophage" in cats and opts.phage_tool == "genomad":
        db = genomad_db(opts.genomad_db)
        versions["genomad_db"] = genomad_db_version(db)
        # No path in the message: it may come from an environment variable --redact cannot hide.
        if db is not None and external.writable_by_others(db):
            warnings.append(
                "The geNomad database folder can be modified by any user, so its contents may "
                "not be what was downloaded. Keep it in a read-only shared location."
            )

    # Options that belong to a category are recorded as None when that category did not run.
    params = {
        "categories": cats,
        "read_length": opts.read_length,
        "repeat_min_len": opts.repeat_min_len,
        "repeat_min_id": opts.repeat_min_id,
        "pad": opts.pad,
        "rrna_operon_gap": opts.rrna_operon_gap if "rrna" in cats else None,
        "isescan_partial": opts.isescan_partial if "is" in cats else None,
        "mask_genes_matching": opts.mask_genes_matching,
        "phage_tool": opts.phage_tool if "prophage" in cats else None,
        "homopolymer_min_len": opts.homopolymer_min_len if "lowcomplexity" in cats else None,
        "threads": threads,
    }

    # Checksums identify the exact reference file, even when --redact hides its path.
    sums = checksums(reference)
    if opts.redact:
        # Reports only: refmasker.log keeps full paths for debugging.
        redact = redactor(input_paths(opts))
        command = redact(command)
        warnings = [redact(w) for w in warnings]
        notes = {k: [redact(n) for n in v] for k, v in notes.items()}
        log.info(
            "--redact: the reports and run_info.json show file names only and omit FASTA "
            "header descriptions. refmasker.log keeps full paths, so share the report files, "
            "not the whole output folder"
        )

    data = RunData(
        reference=Path(reference.name) if opts.redact else reference,
        reference_md5=sums["md5"],
        reference_sha256=sums["sha256"],
        records=records,
        categories=cats,
        regions=regions,
        merged=merged,
        summary=summary,
        sensitivity=sensitivity,
        params=params,
        versions=versions,
        tool_paths=tool_paths,
        warnings=warnings,
        notes=notes,
        plasmids=plasmids,
        annotation=Path(annotation.name) if (opts.redact and annotation) else annotation,
        redacted=opts.redact,
        command=command,
        started=started,
    )

    for w in warnings:
        log.warning(w)

    # The genome-wide total row of the summary.
    total = next(r for r in summary if r.record == "ALL" and r.category == "total")
    log.info(
        "Masked %s of %s bp (%.2f%%) in %d regions",
        f"{total.bp:,}",
        f"{total.length:,}",
        total.pct,
        len(merged),
    )

    # ---- 7. Reports and clean-up ----
    data.runtime_s = time.time() - t0
    from refmasker.report import write_reports

    write_reports(data, outdir, opts.formats)
    write_run_info(data, outdir / "run_info.json")
    log.info("Done in %.0f s. Mask: %s", time.time() - t0, outdir / "mask.bed")
    return data


def write_run_info(data: RunData, path: Path) -> None:
    """Write run_info.json: a machine-readable record of inputs, parameters, versions and result.

    FASTA header descriptions are left out with --redact.
    """
    total = next(r for r in data.summary if r.record == "ALL" and r.category == "total")
    info = {
        "refmasker_version": __version__,
        "command": data.command,
        "started": data.started,
        "runtime_s": round(data.runtime_s, 1),
        "reference": str(data.reference),
        "reference_md5": data.reference_md5,
        "reference_sha256": data.reference_sha256,
        "annotation": str(data.annotation) if data.annotation else None,
        "redacted": data.redacted,
        "records": [
            {"id": r.id, "length": r.length}
            | ({} if data.redacted else {"description": r.description})
            for r in data.records
        ],
        "plasmids": data.plasmids,
        "params": data.params,
        "versions": data.versions,
        # Full paths of the tools used; left out with --redact (they reveal folder layouts).
        "tool_paths": None if data.redacted else data.tool_paths,
        "masked_bp": total.bp,
        "masked_pct": round(total.pct, 4),
        "warnings": data.warnings,
    }
    with open_replacing(path) as fh:
        fh.write(json.dumps(info, indent=2) + "\n")
