"""Repeats: regions present in more than one copy, found by aligning the genome to itself."""

from __future__ import annotations

import logging
from dataclasses import dataclass

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)

# Sensitivity grid reported alongside the chosen thresholds.
SENSITIVITY_LENGTHS = (100, 200, 300, 500)
SENSITIVITY_IDENTITIES = (80.0, 90.0, 95.0)


@dataclass(frozen=True)
class RepeatHit:
    chrom1: str
    start1: int  # 0-based, exclusive end
    end1: int
    chrom2: str
    start2: int
    end2: int
    length: int  # shorter of the two aligned lengths
    identity: float


def parse_show_coords(text: str) -> list[RepeatHit]:
    """Parse `show-coords -T -H -r -c -l` output, dropping each sequence's self-hit."""
    hits = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 13:
            continue
        s1, e1, s2, e2 = (int(x) for x in cols[0:4])
        len1, len2 = int(cols[4]), int(cols[5])
        idy = float(cols[6])
        tag1, tag2 = cols[11], cols[12]
        if tag1 == tag2 and s1 == s2 and e1 == e2:
            continue  # the trivial alignment of a sequence to itself
        a1, b1 = min(s1, e1), max(s1, e1)
        a2, b2 = min(s2, e2), max(s2, e2)
        hits.append(RepeatHit(tag1, a1 - 1, b1, tag2, a2 - 1, b2, min(len1, len2), idy))
    return hits


def hits_to_regions(hits: list[RepeatHit], min_len: int, min_id: float) -> list[Region]:
    """Both copies of every hit passing the thresholds become regions.

    A multi-copy element produces one alignment per pair of copies, so overlapping
    intervals are merged into one region per locus, named by how many other copies
    it matched and the highest identity.
    """
    pieces: list[tuple[str, int, int, str, float]] = []
    for h in hits:
        if h.length < min_len or h.identity < min_id:
            continue
        pieces.append((h.chrom1, h.start1, h.end1, f"{h.chrom2}:{h.start2}", h.identity))
        pieces.append((h.chrom2, h.start2, h.end2, f"{h.chrom1}:{h.start1}", h.identity))
    pieces.sort()
    regions: list[Region] = []
    group: list[tuple[str, int, int, str, float]] = []
    for piece in pieces:
        if group and (piece[0] != group[0][0] or piece[1] >= max(g[2] for g in group)):
            regions.append(_locus(group))
            group = []
        group.append(piece)
    if group:
        regions.append(_locus(group))
    return regions


def _locus(group) -> Region:
    chrom = group[0][0]
    start = min(g[1] for g in group)
    end = max(g[2] for g in group)
    partners = len({g[3] for g in group})
    best = max(g[4] for g in group)
    word = "copy" if partners == 1 else "copies"
    name = f"matches {partners} other {word}, up to {best:.1f}% identity"
    return Region(chrom, start, end, "repeats", "nucmer", name)


def detect(ctx: Context) -> Result:
    opts = ctx.options
    wd = ctx.workdir / "repeats"
    wd.mkdir(parents=True, exist_ok=True)
    prefix = wd / "self"
    external.run(
        [
            "nucmer",
            "--maxmatch",
            "--nosimplify",
            "--threads",
            str(ctx.threads),
            "-p",
            str(prefix),
            str(ctx.fasta),
            str(ctx.fasta),
        ],
        wd,
        "nucmer",
    )
    coords = wd / "self.coords.tsv"
    external.run(
        ["show-coords", "-T", "-H", "-r", "-c", "-l", f"{prefix}.delta"],
        wd,
        "show-coords",
        stdout_path=coords,
    )
    hits = parse_show_coords(coords.read_text())
    regions = hits_to_regions(hits, opts.repeat_min_len, opts.repeat_min_id)
    log.info(
        "repeats: %d alignments, %d regions at >=%d bp and >=%.0f%% identity",
        len(hits),
        len(regions),
        opts.repeat_min_len,
        opts.repeat_min_id,
    )
    return Result(regions=regions, extra={"hits": hits})
