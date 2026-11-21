"""Prophages: geNomad (default) or PhiSpy."""

from __future__ import annotations

import csv
import logging
import os
from pathlib import Path

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)

# Environment variables checked, in this order, when --genomad-db is not given.
GENOMAD_DB_ENV = ("REFMASKER_GENOMAD_DB", "GENOMAD_DB")


def genomad_db(cli_value: str | None) -> Path | None:
    """The geNomad database folder: --genomad-db, else the first variable set, else None."""
    if cli_value:
        return Path(cli_value)
    for var in GENOMAD_DB_ENV:
        if os.environ.get(var):
            return Path(os.environ[var])
    return None


def genomad_db_version(db: Path | None) -> str:
    """Version recorded in the geNomad database folder's version.txt, for provenance.

    refmasker cannot checksum the download itself (geNomad fetches it), so recording the
    version is how a run says which database it used.
    """
    if db is None:
        return "not set"
    vfile = db / "version.txt"
    if not vfile.is_file():
        return "unknown (no version.txt)"
    lines = vfile.read_text().split()
    return lines[0] if lines else "unknown (empty version.txt)"


def parse_genomad_virus_summary(path: Path, lengths: dict[str, int]) -> list[Region]:
    """Parse geNomad's *_virus_summary.tsv.

    Proviruses are named `contig|provirus_START_END`; a sequence without that suffix was
    classified as viral end to end and is masked whole.
    """
    regions = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq_name = row["seq_name"]
            score = row.get("virus_score", "")
            taxonomy = row.get("taxonomy", "")
            tax_short = [t for t in taxonomy.split(";") if t]
            # Region name: the virus score plus the most specific taxon geNomad assigned.
            name = f"score={score}" + (f" {tax_short[-1]}" if tax_short else "")
            if "|provirus_" in seq_name:
                # A provirus inside a longer sequence. Use the 'coordinates' column (START-END,
                # 1-based) if present, otherwise the coordinates in the name.
                contig = seq_name.split("|provirus_")[0]
                coords = row.get("coordinates", "")
                if "-" in coords:
                    s, e = (int(x) for x in coords.split("-"))
                else:
                    s, e = (int(x) for x in seq_name.split("|provirus_")[1].split("_")[:2])
                regions.append(Region(contig, s - 1, e, "prophage", "geNomad", name))
            elif seq_name in lengths:
                # The whole sequence was classified as viral: mask all of it.
                regions.append(
                    Region(
                        seq_name,
                        0,
                        lengths[seq_name],
                        "prophage",
                        "geNomad",
                        name + " (whole sequence)",
                    )
                )
    return regions


def parse_phispy_coordinates(path: Path) -> list[Region]:
    """Parse PhiSpy's prophage_coordinates.tsv (id, contig, start, stop, ...)."""
    regions = []
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            # Skip the header line and any malformed line (start must be a number).
            if len(cols) < 4 or not cols[2].isdigit():
                continue
            s, e = int(cols[2]), int(cols[3])
            regions.append(Region(cols[1], min(s, e) - 1, max(s, e), "prophage", "PhiSpy", cols[0]))
    return regions


def detect(ctx: Context) -> Result:
    """Run the prophage caller chosen with --phage-tool."""
    if ctx.options.phage_tool == "phispy":
        return _phispy(ctx)
    return _genomad(ctx)


def _genomad(ctx: Context) -> Result:
    """Run geNomad's end-to-end pipeline and parse its virus summary table.

    --cleanup makes geNomad delete its large intermediate files.
    """
    # Check the database before starting, so a missing one gives a clear message.
    db = genomad_db(ctx.options.genomad_db)
    if db is None or not db.is_dir():
        raise external.ToolError(
            "geNomad needs its database. Pass --genomad-db DIR or set REFMASKER_GENOMAD_DB. "
            "Download it once with `refmasker download-db DIR`, or use --phage-tool phispy."
        )
    wd = ctx.workdir / "genomad"
    wd.mkdir(parents=True, exist_ok=True)
    external.run(
        [
            "genomad",
            "end-to-end",
            "--cleanup",
            "--threads",
            str(ctx.threads),
            str(ctx.fasta),
            str(wd / "out"),
            str(db),
        ],
        wd,
        "genomad",
    )
    # The summary is in a subfolder named after the input file, so search for it.
    summaries = sorted((wd / "out").rglob("*_virus_summary.tsv"))
    if not summaries:
        raise external.ToolError(f"geNomad finished but no *_virus_summary.tsv in {wd / 'out'}")
    regions = parse_genomad_virus_summary(summaries[0], ctx.lengths)
    log.info("prophage: geNomad found %d regions", len(regions))
    return Result(regions=regions)


def _phispy(ctx: Context) -> Result:
    """Run PhiSpy, which needs an annotated GenBank file (it scores genes, not raw sequence).

    --output_choice 1 asks PhiSpy for prophage_coordinates.tsv only.
    """
    ann = ctx.annotation_path
    if ann is None or ctx.options.annotation_format != "genbank":
        raise external.ToolError(
            "PhiSpy needs an annotated GenBank file: pass one with --annotation, "
            "or use the default --phage-tool genomad."
        )
    wd = ctx.workdir / "phispy"
    wd.mkdir(parents=True, exist_ok=True)
    external.run(
        [
            "PhiSpy.py",
            str(ann),
            "-o",
            str(wd / "out"),
            "--threads",
            str(ctx.threads),
            "--output_choice",
            "1",
        ],
        wd,
        "phispy",
    )
    coords = wd / "out" / "prophage_coordinates.tsv"
    regions = parse_phispy_coordinates(coords) if coords.exists() else []
    log.info("prophage: PhiSpy found %d regions", len(regions))
    return Result(regions=regions)
