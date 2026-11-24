"""Genomic intervals: merging, padding, clipping and BED I/O.

All coordinates are BED-style: 0-based start, exclusive end.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path

from refmasker.safeio import open_replacing


@dataclass(frozen=True)
class Region:
    """One masked interval and why it was masked.

    `category` is the mask category (e.g. "repeats"), `source` the tool or file that found it
    (e.g. "nucmer"), `name` a short human-readable label, and `genes` the annotated genes it
    overlaps (filled in later by the pipeline, and ignored when comparing two Regions).
    """

    chrom: str
    start: int
    end: int
    category: str = ""
    source: str = ""
    name: str = ""
    genes: tuple[str, ...] = field(default=(), compare=False)

    @property
    def length(self) -> int:
        """Length in bp (the end is exclusive, so no +1 is needed)."""
        return self.end - self.start


@dataclass(frozen=True)
class MergedRegion:
    """An interval made by merging overlapping Regions; remembers which categories it contains."""

    chrom: str
    start: int
    end: int
    categories: tuple[str, ...]

    @property
    def length(self) -> int:
        return self.end - self.start


def merge(regions: Iterable[Region], gap: int = 0) -> list[MergedRegion]:
    """Merge overlapping (or book-ended, or within `gap` bp) regions per chromosome."""
    # Group by chromosome, skipping empty (zero-length) intervals.
    by_chrom: dict[str, list[Region]] = defaultdict(list)
    for r in regions:
        if r.end > r.start:
            by_chrom[r.chrom].append(r)
    out: list[MergedRegion] = []
    for chrom, items in by_chrom.items():
        # Sweep left to right: extend the current interval while the next one starts inside it
        # (or within `gap` bp of its end); otherwise close it and start a new one. Empty
        # category names (e.g. from a plain BED file) are left out of the category list.
        items.sort(key=lambda r: (r.start, r.end))
        cur_s, cur_e = items[0].start, items[0].end
        cats = {items[0].category}
        for r in items[1:]:
            if r.start <= cur_e + gap:
                cur_e = max(cur_e, r.end)
                cats.add(r.category)
            else:
                out.append(MergedRegion(chrom, cur_s, cur_e, tuple(sorted(c for c in cats if c))))
                cur_s, cur_e, cats = r.start, r.end, {r.category}
        out.append(MergedRegion(chrom, cur_s, cur_e, tuple(sorted(c for c in cats if c))))
    return out


def sort_by_record(items, record_order: list[str]):
    """Sort by record in FASTA order, then by position. Unknown records go last, by name."""
    rank = {name: i for i, name in enumerate(record_order)}
    return sorted(items, key=lambda r: (rank.get(r.chrom, len(rank)), r.chrom, r.start, r.end))


def union_bp(regions: Iterable[Region]) -> int:
    """Number of bases covered by at least one region (overlapping bases counted once)."""
    return sum(m.length for m in merge(regions))


def pad_and_clip(regions: Iterable[Region], pad: int, lengths: Mapping[str, int]) -> list[Region]:
    """Extend every region by `pad` bp on each side and clip to record bounds."""
    out = []
    for r in regions:
        length = lengths[r.chrom]
        s = max(0, r.start - pad)
        e = min(length, r.end + pad)
        if e > s:  # drop regions that are empty after clipping
            out.append(replace(r, start=s, end=e))
    return out


def _by_chrom(regions: Iterable[Region]) -> dict[str, list[MergedRegion]]:
    """Merge `regions` and group them by chromosome, each group sorted by start."""
    # merge() sorts each chromosome's intervals by start and emits one chromosome at a time,
    # so appending in that order keeps every group sorted.
    out: dict[str, list[MergedRegion]] = defaultdict(list)
    for m in merge(regions):
        out[m.chrom].append(m)
    return out


def overlaps(a: Iterable[Region], b: Iterable[Region]) -> Iterator[tuple[str, int, int]]:
    """Yield (chrom, start, end) for every stretch covered by both region sets.

    Each set is merged first, so overlaps within one set are not yielded twice. Lazy, so
    counting bases does not build a list.
    """
    by_a = _by_chrom(a)
    by_b = _by_chrom(b)
    for chrom, xs in by_a.items():
        ys = by_b.get(chrom, [])
        # Two-pointer sweep along the chromosome: yield the overlap of the current pair, then
        # move past whichever interval ends first (it cannot overlap anything further along).
        i = j = 0
        while i < len(xs) and j < len(ys):
            s = max(xs[i].start, ys[j].start)
            e = min(xs[i].end, ys[j].end)
            if e > s:
                yield chrom, s, e
            if xs[i].end < ys[j].end:
                i += 1
            else:
                j += 1


def intersect_bp(a: Iterable[Region], b: Iterable[Region]) -> int:
    """Number of bases covered by both region sets."""
    return sum(e - s for _chrom, s, e in overlaps(a, b))


def intersect(
    a: Iterable[Region], b: Iterable[Region], category: str = "", source: str = ""
) -> list[Region]:
    """The intervals covered by both region sets."""
    return [Region(c, s, e, category, source) for c, s, e in overlaps(a, b)]


def subtract(
    a: Iterable[Region], b: Iterable[Region], category: str = "", source: str = ""
) -> list[Region]:
    """The parts of `a` that `b` does not cover.

    This is a per-base difference, so a region of `a` that `b` covers halfway comes back
    halved. That is different from "regions of `a` that `b` misses entirely", which is what
    `compare` lists as unique regions.
    """
    by_b = _by_chrom(b)
    # One pointer per chromosome, carried between regions of `a`: both sides are sorted, so
    # it only ever moves forward and the whole sweep stays linear.
    start_at: dict[str, int] = {}
    out: list[Region] = []
    for m in merge(a):
        ys = by_b.get(m.chrom, [])
        j = start_at.get(m.chrom, 0)
        while j < len(ys) and ys[j].end <= m.start:
            j += 1
        start_at[m.chrom] = j
        pos = m.start
        k = j
        while k < len(ys) and ys[k].start < m.end:
            if ys[k].start > pos:
                out.append(Region(m.chrom, pos, ys[k].start, category, source))
            pos = max(pos, ys[k].end)
            k += 1
        if pos < m.end:
            out.append(Region(m.chrom, pos, m.end, category, source))
    return out


def overlaps_any(region: Region, others: Iterable[Region]) -> bool:
    """True if `region` overlaps any of `others`. Intervals that only touch do not overlap."""
    return any(
        o.chrom == region.chrom and o.start < region.end and o.end > region.start for o in others
    )


def read_bed(path: str | Path, category: str = "user", source: str = "bed") -> list[Region]:
    """Read a BED file (3+ columns). Column 4, if present, becomes the region name."""
    regions = []
    with open(path) as fh:
        # Count lines from 1 so error messages match the line numbers in a text editor.
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n\r")
            # Skip blank lines, comments and UCSC track/browser header lines.
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise ValueError(
                    f"{path}:{lineno}: BED lines need at least 3 tab-separated columns"
                )
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: start/end are not integers") from exc
            if start < 0 or end < start:
                raise ValueError(f"{path}:{lineno}: invalid interval {start}-{end}")
            name = cols[3] if len(cols) > 3 else ""
            regions.append(Region(cols[0], start, end, category, source, name))
    return regions


def write_mask_bed(path: str | Path, merged: Iterable[MergedRegion]) -> None:
    """Write the snippy-ready BED: 3 columns, no header (snippy-core rejects header lines)."""
    with open_replacing(path) as fh:
        for m in merged:
            fh.write(f"{m.chrom}\t{m.start}\t{m.end}\n")


# Columns of mask.annotated.bed: one row per region before merging, with why it was masked.
ANNOTATED_HEADER = ["chrom", "start", "end", "category", "source", "length", "name", "genes"]


def write_annotated_bed(path: str | Path, regions: Iterable[Region]) -> None:
    """Write every unmerged region with its category, source, name and overlapped genes.

    Empty fields are written as '.', the usual BED placeholder.
    """
    with open_replacing(path) as fh:
        fh.write("#" + "\t".join(ANNOTATED_HEADER) + "\n")
        for r in regions:
            fh.write(
                "\t".join(
                    [
                        r.chrom,
                        str(r.start),
                        str(r.end),
                        r.category,
                        r.source,
                        str(r.length),
                        r.name or ".",
                        ",".join(r.genes) or ".",
                    ]
                )
                + "\n"
            )


def _unplaceholder(value: str) -> str:
    """Turn the '.' that write_annotated_bed writes for an empty field back into ''."""
    return "" if value == "." else value


def is_annotated_bed(path: str | Path) -> bool:
    """True if `path` looks like a refmasker mask.annotated.bed (by its header line)."""
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                # The header is the first non-blank line that write_annotated_bed writes.
                cols = line.lstrip("#").rstrip("\n").split("\t")
                return cols[:4] == ANNOTATED_HEADER[:4]
    except OSError:
        return False
    return False


def read_annotated_bed(path: str | Path) -> list[Region]:
    """Read a mask.annotated.bed back, keeping each region's category, source and name.

    The length column is recomputed from the coordinates rather than trusted, and '.' means
    an empty field (as written by write_annotated_bed).
    """
    regions = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                raise ValueError(
                    f"{path}:{lineno}: an annotated BED needs at least the columns "
                    f"{', '.join(ANNOTATED_HEADER[:5])}"
                )
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: start/end are not integers") from exc
            if start < 0 or end < start:
                raise ValueError(f"{path}:{lineno}: invalid interval {start}-{end}")
            name = _unplaceholder(cols[6]) if len(cols) > 6 else ""
            genes = _unplaceholder(cols[7]) if len(cols) > 7 else ""
            regions.append(
                Region(
                    cols[0],
                    start,
                    end,
                    _unplaceholder(cols[3]),
                    _unplaceholder(cols[4]),
                    name,
                    tuple(g for g in genes.split(",") if g),
                )
            )
    return regions
