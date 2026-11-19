"""The `refmasker run` workflow: detect, combine, summarise, report."""

from __future__ import annotations

import json
import logging
import platform
import shlex
import shutil
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
from refmasker.seqio import Record, checksums, read_fasta, write_fasta
from refmasker.summary import (
    SensitivityRow,
    SummaryRow,
    repeat_sensitivity,
    summarise,
    write_summary_tsv,
)

log = logging.getLogger(__name__)

MAX_GENES_LISTED = 50


@dataclass
class RunData:
    """Everything the report needs."""

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
        return sum(r.length for r in self.records)


def resolve_categories(
    spec: str, extra_bed: bool, mask_plasmids: bool, gene_regex: bool = False
) -> list[str]:
    """Turn '--categories default,crispr' into an ordered list of category names."""
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
        chosen.extend(n for n in names if n not in chosen)
    if extra_bed and "user" not in chosen:
        chosen.append("user")
    if mask_plasmids and "plasmid" not in chosen:
        chosen.append("plasmid")
    if gene_regex and "genes" not in chosen:
        chosen.append("genes")
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


def prepare_outdir(outdir: Path, force: bool) -> None:
    if outdir.exists() and any(outdir.iterdir()):
        if not force:
            raise FileExistsError(f"Output folder {outdir} exists and is not empty (use --force)")
        if not any((outdir / m).is_file() for m in RUN_MARKERS):
            raise FileExistsError(
                f"Output folder {outdir} is not empty and does not look like a refmasker output "
                f"folder (no {' or '.join(RUN_MARKERS)}), so --force will not overwrite it. "
                "Choose a new or empty folder."
            )
        for name in ("work", "figures"):
            shutil.rmtree(outdir / name, ignore_errors=True)
    outdir.mkdir(parents=True, exist_ok=True)


def add_file_logging(outdir: Path) -> logging.Handler:
    handler = logging.FileHandler(outdir / "refmasker.log", mode="w")
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
    t0 = time.time()
    started = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    outdir = Path(opts.outdir)
    prepare_outdir(outdir, opts.force)
    handler = add_file_logging(outdir)
    try:
        return _run(opts, outdir, t0, started)
    finally:
        logging.getLogger().removeHandler(handler)
        handler.close()


def _run(opts, outdir: Path, t0: float, started: str) -> RunData:
    log.info("refmasker %s", __version__)
    command = "refmasker " + " ".join(shlex.quote(a) for a in sys.argv[1:])
    log.info("Command: %s", command)

    threads = resolve_threads(opts.threads)
    limit_library_threads(threads)
    log.info("Threads: %d", threads)

    cats = resolve_categories(
        opts.categories, bool(opts.extra_bed), opts.mask_plasmids, bool(opts.mask_genes_matching)
    )
    log.info("Categories: %s", ", ".join(cats))
    registry = categories()

    needed = sorted({t for c in cats for t in registry[c].tools(opts)})
    missing = external.missing_tools(needed)
    if missing:
        raise external.ToolError(
            "Missing tools on PATH: "
            + ", ".join(missing)
            + ". Install them with `conda env create -f environment.yml`, "
            "or drop the categories that need them with --categories."
        )

    reference = Path(opts.reference)
    records = read_fasta(reference)
    lengths = {r.id: r.length for r in records}
    log.info(
        "Reference: %s, %d record(s), %s bp", reference, len(records), f"{sum(lengths.values()):,}"
    )

    warnings: list[str] = []
    plasmids = plasmid_ids(records, opts.plasmid_ids)
    chromosomal = [r for r in records if r.id not in plasmids]
    if len(chromosomal) > 1:
        warnings.append(
            f"The reference has {len(chromosomal)} non-plasmid records. If this is a draft "
            "assembly, repeats and IS elements are often collapsed or split at contig ends, so "
            "they will be under-detected. A complete genome is recommended."
        )
    n_count = sum(r.seq.upper().count("N") for r in records)
    if n_count:
        warnings.append(f"The reference contains {n_count:,} N bases.")

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

    workdir = outdir / "work"
    workdir.mkdir(exist_ok=True)
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

    regions = pad_and_clip(raw, opts.pad, lengths)
    regions = annotate_genes(regions, index)
    order = [r.id for r in records]
    regions = sort_by_record(regions, order)
    merged = sort_by_record(merge(regions), order)

    write_mask_bed(outdir / "mask.bed", merged)
    write_annotated_bed(outdir / "mask.annotated.bed", regions)
    summary = summarise(regions, records, cats)
    write_summary_tsv(outdir / "summary.tsv", summary)

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

    versions = {"refmasker": __version__, "python": platform.python_version()}
    for tool in needed:
        versions[tool] = external.tool_version(tool)
    if "prophage" in cats and opts.phage_tool == "genomad":
        versions["genomad_db"] = genomad_db_version(genomad_db(opts.genomad_db))

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

    sums = checksums(reference)
    if opts.redact:
        # Reports only: refmasker.log keeps full paths for debugging.
        redact = redactor(input_paths(opts))
        command = redact(command)
        warnings = [redact(w) for w in warnings]
        notes = {k: [redact(n) for n in v] for k, v in notes.items()}
        log.info("--redact: reports show file names only and omit FASTA header descriptions")

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

    total = next(r for r in summary if r.record == "ALL" and r.category == "total")
    log.info(
        "Masked %s of %s bp (%.2f%%) in %d regions",
        f"{total.bp:,}",
        f"{total.length:,}",
        total.pct,
        len(merged),
    )

    data.runtime_s = time.time() - t0
    from refmasker.report import write_reports

    write_reports(data, outdir, opts.formats)
    write_run_info(data, outdir / "run_info.json")

    if not opts.keep_work:
        shutil.rmtree(workdir, ignore_errors=True)
    log.info("Done in %.0f s. Mask: %s", time.time() - t0, outdir / "mask.bed")
    return data


def write_run_info(data: RunData, path: Path) -> None:
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
        "masked_bp": total.bp,
        "masked_pct": round(total.pct, 4),
        "warnings": data.warnings,
    }
    path.write_text(json.dumps(info, indent=2) + "\n")
