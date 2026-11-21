"""CRISPR arrays (MinCED)."""

from __future__ import annotations

import logging

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)


def parse_minced_gff(text: str) -> list[Region]:
    """Parse MinCED's GFF: one 'repeat_region' line per CRISPR array.

    The score column (6th) is reported as the number of repeats in the array.
    """
    regions = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) != 9 or cols[2] != "repeat_region":
            continue
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        name = f"{attrs.get('ID', 'CRISPR')} ({cols[5]} repeats)"
        # GFF is 1-based inclusive; BED start is 0-based, so subtract 1 from the start only.
        regions.append(Region(cols[0], int(cols[3]) - 1, int(cols[4]), "crispr", "MinCED", name))
    return regions


def detect(ctx: Context) -> Result:
    """Run MinCED on the reference; one region per CRISPR array."""
    wd = ctx.workdir / "minced"
    wd.mkdir(parents=True, exist_ok=True)
    gff = wd / "minced.gff"
    # MinCED writes a text report and a GFF; only the GFF is used, and only if it was written.
    external.run(["minced", str(ctx.fasta), str(wd / "minced.txt"), str(gff)], wd, "minced")
    regions = parse_minced_gff(gff.read_text()) if gff.exists() else []
    log.info("crispr: %d arrays", len(regions))
    return Result(regions=regions)
