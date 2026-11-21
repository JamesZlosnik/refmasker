"""Insertion sequences and transposons: ISEScan, plus transposase genes from the annotation."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from refmasker import external
from refmasker.annotation import Feature
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)


def parse_isescan_tsv(path: Path, include_partial: bool = False) -> list[Region]:
    """Parse ISEScan's per-genome .tsv (columns seqID, family, cluster, isBegin, isEnd, ...).

    ISEScan marks each call complete ('c') or partial ('p'). Partial calls include genuine
    IS remnants but also more false positives, so they are skipped unless asked for.
    """
    regions = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Strip stray whitespace from column names and values.
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            # Skip rows without usable coordinates.
            try:
                s, e = int(row["isBegin"]), int(row["isEnd"])
            except (KeyError, ValueError):
                continue
            name = "/".join(x for x in (row.get("family", ""), row.get("cluster", "")) if x)
            if row.get("type") == "p" and not include_partial:
                continue
            kind = {"c": "complete", "p": "partial"}.get(row.get("type", ""), "")
            if kind:
                name = f"{name} ({kind})"
            # Order the coordinates in case begin > end, and convert to BED (0-based start).
            regions.append(Region(row["seqID"], min(s, e) - 1, max(s, e), "is", "ISEScan", name))
    return regions


def find_isescan_table(outdir: Path) -> Path | None:
    """Find ISEScan's per-genome results table under `outdir`, or None if there is none.

    ISEScan writes into subfolders that depend on the input file's path, so look for the .tsv
    whose header has the expected columns rather than a fixed file name.
    """
    for p in sorted(outdir.rglob("*.tsv")):
        with open(p) as fh:
            header = fh.readline()
        if "isBegin" in header and "seqID" in header:
            return p
    return None


def annotation_regions(features: list[Feature]) -> list[Region]:
    """Regions for annotated mobile elements and transposase/IS genes (see Feature.is_mobile)."""
    regions = []
    for f in features:
        if f.is_mobile:
            name = f.product or f.ftype
            if f.locus_tag:
                name = f"{f.locus_tag} {name}"
            regions.append(Region(f.seqid, f.start, f.end, "is", "annotation", name))
    return regions


def detect(ctx: Context) -> Result:
    """Run ISEScan and, with an annotation, add its transposase/IS features.

    The two sources often find the same elements; overlaps are merged when mask.bed is built.
    """
    wd = ctx.workdir / "isescan"
    wd.mkdir(parents=True, exist_ok=True)
    out = wd / "out"
    # Absolute paths, because ISEScan runs with work/isescan as its working directory.
    external.run(
        [
            "isescan.py",
            "--seqfile",
            str(ctx.fasta.resolve()),
            "--output",
            str(out.resolve()),
            "--nthread",
            str(ctx.threads),
        ],
        wd,
        "isescan",
        cwd=wd,
    )
    result = Result()
    table = find_isescan_table(out)
    if table is None:
        log.info("is: ISEScan found no insertion sequences")
        result.notes.append("ISEScan found no insertion sequences.")
    else:
        result.regions.extend(parse_isescan_tsv(table, ctx.options.isescan_partial))
    if ctx.features is not None:
        ann = annotation_regions(ctx.features)
        result.regions.extend(ann)
        log.info("is: %d transposase/IS features from the annotation", len(ann))
    log.info("is: %d regions in total", len(result.regions))
    return result
