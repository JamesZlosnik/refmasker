"""CRISPR arrays (MinCED)."""

from __future__ import annotations

import logging

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)


def parse_minced_gff(text: str) -> list[Region]:
    regions = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) != 9 or cols[2] != "repeat_region":
            continue
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        name = f"{attrs.get('ID', 'CRISPR')} ({cols[5]} repeats)"
        regions.append(Region(cols[0], int(cols[3]) - 1, int(cols[4]), "crispr", "MinCED", name))
    return regions


def detect(ctx: Context) -> Result:
    wd = ctx.workdir / "minced"
    wd.mkdir(parents=True, exist_ok=True)
    gff = wd / "minced.gff"
    external.run(["minced", str(ctx.fasta), str(wd / "minced.txt"), str(gff)], wd, "minced")
    regions = parse_minced_gff(gff.read_text()) if gff.exists() else []
    log.info("crispr: %d arrays", len(regions))
    return Result(regions=regions)
