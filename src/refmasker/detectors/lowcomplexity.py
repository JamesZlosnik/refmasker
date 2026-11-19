"""Low-complexity sequence (dustmasker) and homopolymer runs."""

from __future__ import annotations

import logging
import re

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region
from refmasker.seqio import Record

log = logging.getLogger(__name__)


def parse_dust_intervals(text: str) -> list[Region]:
    """Parse `dustmasker -outfmt interval` (0-based inclusive `start - end` lines)."""
    regions = []
    chrom = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            chrom = line[1:].split()[0]
        elif " - " in line and chrom:
            s, e = (int(x) for x in line.split(" - "))
            regions.append(Region(chrom, s, e + 1, "lowcomplexity", "dustmasker", "dust"))
    return regions


def find_homopolymers(records: list[Record], min_len: int) -> list[Region]:
    if min_len <= 1:
        raise ValueError("--homopolymer-min-len must be at least 2")
    pattern = re.compile(r"([ACGTacgt])\1{" + str(min_len - 1) + ",}")
    regions = []
    for rec in records:
        for m in pattern.finditer(rec.seq):
            base = m.group(1).upper()
            regions.append(
                Region(
                    rec.id,
                    m.start(),
                    m.end(),
                    "lowcomplexity",
                    "homopolymer",
                    f"poly{base} x{m.end() - m.start()}",
                )
            )
    return regions


def detect(ctx: Context) -> Result:
    wd = ctx.workdir / "dustmasker"
    wd.mkdir(parents=True, exist_ok=True)
    out = wd / "dust.txt"
    external.run(
        ["dustmasker", "-in", str(ctx.fasta), "-outfmt", "interval"],
        wd,
        "dustmasker",
        stdout_path=out,
    )
    regions = parse_dust_intervals(out.read_text())
    homo = find_homopolymers(ctx.records, ctx.options.homopolymer_min_len)
    log.info("lowcomplexity: %d dust intervals, %d homopolymers", len(regions), len(homo))
    return Result(regions=regions + homo)
