"""Compare two BED masks (e.g. refmasker vs a published mask)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from refmasker.regions import Region, intersect, intersect_bp, merge, read_bed, subtract
from refmasker.seqio import read_fasta


@dataclass(frozen=True)
class RecordComparison:
    """How the two masks compare on one sequence: a row of the per-record table."""

    record: str
    bp_a: int
    bp_b: int
    bp_both: int
    regions_a: int
    regions_b: int
    # Record length, when a reference FASTA was given. Without one the BEDs do not say how
    # long the sequence is, so percentages are not available and `extent` stands in for the
    # axis of the figures.
    length: int | None
    extent: int  # furthest masked base on this record, in either mask

    @property
    def bp_only_a(self) -> int:
        """Bases masked by A and not by B."""
        return self.bp_a - self.bp_both

    @property
    def bp_only_b(self) -> int:
        return self.bp_b - self.bp_both

    @property
    def span(self) -> int:
        """Length to draw this record over: the real length if known, else the masked extent."""
        return self.length or self.extent

    @property
    def jaccard(self) -> float:
        union = self.bp_a + self.bp_b - self.bp_both
        return self.bp_both / union if union else 1.0

    def pct(self, bp: int) -> float | None:
        """`bp` as a percentage of the record, or None when the length is unknown."""
        return 100.0 * bp / self.length if self.length else None


@dataclass
class Comparison:
    """Overlap between mask A and mask B, in bp and in merged regions.

    `only_a` / `only_b` are merged regions of one mask that share no base with the other:
    regions the other mask misses completely. `diff_a` / `diff_b` are the per-base
    differences, so a region the other mask half covers appears there halved.
    """

    bp_a: int
    bp_b: int
    bp_both: int
    regions_a: int
    regions_b: int
    only_a: list[Region]
    only_b: list[Region]
    merged_a: list[Region] = field(default_factory=list)
    merged_b: list[Region] = field(default_factory=list)
    shared: list[Region] = field(default_factory=list)
    diff_a: list[Region] = field(default_factory=list)
    diff_b: list[Region] = field(default_factory=list)
    per_record: list[RecordComparison] = field(default_factory=list)
    lengths: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    @property
    def bp_only_a(self) -> int:
        """Bases masked by A and not by B."""
        return self.bp_a - self.bp_both

    @property
    def bp_only_b(self) -> int:
        return self.bp_b - self.bp_both

    @property
    def genome_length(self) -> int | None:
        """Total length of the records in the reference, or None without one."""
        return sum(self.lengths.values()) or None

    def pct(self, bp: int) -> float | None:
        """`bp` as a percentage of the whole reference, or None without a reference."""
        total = self.genome_length
        return 100.0 * bp / total if total else None

    @property
    def jaccard(self) -> float:
        """bp in both / bp in either (1.0 = identical; two empty masks count as identical)."""
        union = self.bp_a + self.bp_b - self.bp_both
        return self.bp_both / union if union else 1.0

    @property
    def frac_b_covered(self) -> float:
        """Fraction of B's bases that A also masks (how much of B that A recovers)."""
        return self.bp_both / self.bp_b if self.bp_b else 0.0

    @property
    def frac_a_in_b(self) -> float:
        """Fraction of A's bases that B also masks (the rest is masked by A alone)."""
        return self.bp_both / self.bp_a if self.bp_a else 0.0


def _as_regions(merged) -> list[Region]:
    """Turn merged intervals back into plain Regions (without categories)."""
    return [Region(m.chrom, m.start, m.end) for m in merged]


def compare(
    a: list[Region], b: list[Region], lengths: dict[str, int] | None = None
) -> Comparison:
    """Compare two region sets. Each is merged first, so overlaps within one mask count once.

    `lengths` (record id -> bp, from a reference FASTA) is optional and only adds the
    percentages and the true figure axis; everything else is computed from the BEDs alone.
    """
    ma, mb = _as_regions(merge(a)), _as_regions(merge(b))
    shared = intersect(ma, mb, "both", "compare")
    # Regions of one mask that the other does not touch at all. Kept alongside the per-base
    # difference below because "this whole region is missing from the other mask" and "these
    # bases differ" answer different questions.
    by_b = {r.chrom for r in mb}
    only_a = [r for r in ma if r.chrom not in by_b or intersect_bp([r], mb) == 0]
    by_a = {r.chrom for r in ma}
    only_b = [r for r in mb if r.chrom not in by_a or intersect_bp([r], ma) == 0]
    lengths = dict(lengths or {})
    return Comparison(
        bp_a=sum(r.length for r in ma),
        bp_b=sum(r.length for r in mb),
        bp_both=sum(r.length for r in shared),
        regions_a=len(ma),
        regions_b=len(mb),
        only_a=only_a,
        only_b=only_b,
        merged_a=ma,
        merged_b=mb,
        shared=shared,
        diff_a=subtract(ma, mb, "only_a", "compare"),
        diff_b=subtract(mb, ma, "only_b", "compare"),
        per_record=_per_record(ma, mb, shared, lengths),
        lengths=lengths,
        warnings=_warnings(ma, mb, lengths),
    )


def _per_record(
    ma: list[Region], mb: list[Region], shared: list[Region], lengths: dict[str, int]
) -> list[RecordComparison]:
    """One row per sequence named by either mask, longest first."""
    rows = []
    names = {r.chrom for r in ma} | {r.chrom for r in mb}
    for name in names:
        ra = [r for r in ma if r.chrom == name]
        rb = [r for r in mb if r.chrom == name]
        rows.append(
            RecordComparison(
                record=name,
                bp_a=sum(r.length for r in ra),
                bp_b=sum(r.length for r in rb),
                bp_both=sum(r.length for r in shared if r.chrom == name),
                regions_a=len(ra),
                regions_b=len(rb),
                length=lengths.get(name),
                extent=max((r.end for r in ra + rb), default=0),
            )
        )
    return sorted(rows, key=lambda r: (-r.span, r.record))


def _warnings(ma: list[Region], mb: list[Region], lengths: dict[str, int]) -> list[str]:
    """Problems worth saying out loud: mismatched names, or regions past the end of a record."""
    out = []
    names_a = {r.chrom for r in ma}
    names_b = {r.chrom for r in mb}
    if names_a and names_b and not (names_a & names_b):
        out.append(
            "The two masks share no sequence names, so nothing can overlap and every number "
            f"below is a comparison of unrelated records. A names {_show(names_a)}; "
            f"B names {_show(names_b)}. Use --rename OLD=NEW if they are the same sequence "
            "under different names."
        )
    elif names_a ^ names_b:
        out.append(
            f"Sequences named by only one mask: {_show(names_a ^ names_b)}. Their regions "
            "count as unique to that mask."
        )
    if lengths:
        unknown = (names_a | names_b) - set(lengths)
        if unknown:
            out.append(
                f"Not in the reference FASTA: {_show(unknown)}. These records have no length, "
                "so they are drawn and reported without percentages."
            )
        past = sorted(
            {r.chrom for r in ma + mb if r.chrom in lengths and r.end > lengths[r.chrom]}
        )
        if past:
            out.append(
                f"Regions run past the end of {_show(set(past))} in the reference FASTA, so "
                "this is probably not the reference these masks were built for."
            )
    return out


def _show(names: set[str], limit: int = 5) -> str:
    """Format a set of sequence names for a message, listing at most `limit` of them."""
    shown = sorted(names)[:limit]
    more = f" (+{len(names) - limit} more)" if len(names) > limit else ""
    return ", ".join(shown) + more


def compare_files(
    a: Path, b: Path, rename: dict[str, str] | None = None, ref: Path | None = None
) -> Comparison:
    """Read two BED files and compare them.

    `rename` maps old to new sequence names in both files, for masks that call the same
    chromosome by different names. `ref` is the reference FASTA the masks were built for,
    which supplies record lengths: without it the comparison is still complete, but there is
    no genome size to express the masks as a percentage of.
    """
    ra, rb = read_bed(a), read_bed(b)
    if rename:
        ra = [Region(rename.get(r.chrom, r.chrom), r.start, r.end) for r in ra]
        rb = [Region(rename.get(r.chrom, r.chrom), r.start, r.end) for r in rb]
    lengths = {r.id: r.length for r in read_fasta(ref)} if ref else None
    return compare(ra, rb, lengths)


def format_report(c: Comparison, name_a: str, name_b: str) -> str:
    """Format the comparison as a text table for the terminal (file names cut to 14 chars)."""
    lines = [
        f"{'':28}{'A':>14}{'B':>14}",
        f"{'file':28}{name_a[-14:]:>14}{name_b[-14:]:>14}",
        f"{'merged regions':28}{c.regions_a:>14,}{c.regions_b:>14,}",
        f"{'bp masked':28}{c.bp_a:>14,}{c.bp_b:>14,}",
        f"{'regions not in the other':28}{len(c.only_a):>14,}{len(c.only_b):>14,}",
        "",
        f"bp masked by both:          {c.bp_both:,}",
        f"fraction of B covered by A: {c.frac_b_covered:.3f}",
        f"fraction of A inside B:     {c.frac_a_in_b:.3f}",
        f"Jaccard (bp):               {c.jaccard:.3f}",
    ]
    return "\n".join(lines)
