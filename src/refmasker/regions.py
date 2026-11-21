"""Genomic intervals: merging, padding, clipping and BED I/O.

All coordinates are BED-style: 0-based start, exclusive end.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
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


def intersect_bp(a: Iterable[Region], b: Iterable[Region]) -> int:
    """Number of bases covered by both region sets."""
    # Merge each set first so overlaps within a set are not counted twice, then group by
    # chromosome. merge() returns each chromosome's intervals sorted by start.
    ma = merge(a)
    mb = merge(b)
    by_a: dict[str, list[MergedRegion]] = defaultdict(list)
    by_b: dict[str, list[MergedRegion]] = defaultdict(list)
    for m in ma:
        by_a[m.chrom].append(m)
    for m in mb:
        by_b[m.chrom].append(m)
    total = 0
    for chrom, xs in by_a.items():
        ys = by_b.get(chrom, [])
        # Two-pointer sweep along the chromosome: add the overlap of the current pair, then
        # move past whichever interval ends first (it cannot overlap anything further along).
        i = j = 0
        while i < len(xs) and j < len(ys):
            s = max(xs[i].start, ys[j].start)
            e = min(xs[i].end, ys[j].end)
            if e > s:
                total += e - s
            if xs[i].end < ys[j].end:
                i += 1
            else:
                j += 1
    return total


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
