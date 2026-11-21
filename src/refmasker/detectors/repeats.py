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
    """One alignment between two copies of a repeat, in BED coordinates."""

    chrom1: str
    start1: int  # 0-based, exclusive end
    end1: int
    chrom2: str
    start2: int
    end2: int
    length: int  # shorter of the two aligned lengths
    identity: float


def parse_show_coords(text: str) -> list[RepeatHit]:
    """Parse `show-coords -T -H -r -c -l` output, dropping each sequence's self-hit.

    -T tab-separated, -H no header, -r sorted by reference, -c coverage columns, -l length
    columns. The 13 columns are S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ TAG1 TAG2,
    with 1-based inclusive coordinates.
    """
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
        # A copy aligned on the reverse strand has start > end: order the coordinates, then
        # convert to BED (0-based start).
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
    # Each passing hit gives two pieces, one per copy:
    # (chrom, start, end, location of the other copy, identity).
    pieces: list[tuple[str, int, int, str, float]] = []
    for h in hits:
        if h.length < min_len or h.identity < min_id:
            continue
        pieces.append((h.chrom1, h.start1, h.end1, f"{h.chrom2}:{h.start2}", h.identity))
        pieces.append((h.chrom2, h.start2, h.end2, f"{h.chrom1}:{h.start1}", h.identity))
    # Sort by chromosome and start, then group overlapping pieces into loci. A new locus
    # starts when the chromosome changes or a piece starts at or after the end of the group.
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
    """Collapse overlapping pieces into one region.

    The name says how many distinct other copies the locus matched, and the best identity.
    """
    chrom = group[0][0]
    start = min(g[1] for g in group)
    end = max(g[2] for g in group)
    partners = len({g[3] for g in group})
    best = max(g[4] for g in group)
    word = "copy" if partners == 1 else "copies"
    name = f"matches {partners} other {word}, up to {best:.1f}% identity"
    return Region(chrom, start, end, "repeats", "nucmer", name)


def detect(ctx: Context) -> Result:
    """Align the genome to itself with nucmer and turn the alignments into repeat regions.

    --maxmatch uses all matches, not just unique ones (repeats are by definition not unique),
    and --nosimplify keeps alignments nucmer would otherwise drop as shadows of the trivial
    self-alignment, which is what MUMmer recommends for finding repeats.
    """
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
    # Convert nucmer's .delta alignment file into a table of alignments.
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
    # The raw hits go back to the pipeline for the threshold sensitivity table.
    return Result(regions=regions, extra={"hits": hits})
