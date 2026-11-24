"""Compare two BED masks (e.g. refmasker vs a published mask)."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path

from refmasker.regions import (
    Region,
    intersect,
    intersect_bp,
    is_annotated_bed,
    merge,
    read_annotated_bed,
    read_bed,
    subtract,
    union_bp,
)
from refmasker.seqio import read_fasta


@dataclass(frozen=True)
class CategoryComparison:
    """How much of one mask category the other mask also covers.

    Categories overlap each other (a region can be both a repeat and an IS element), so these
    rows do not add up to the mask's total.
    """

    category: str
    regions: int
    bp: int
    bp_shared: int

    @property
    def bp_only(self) -> int:
        """Bases of this category the other mask does not cover at all."""
        return self.bp - self.bp_shared

    @property
    def pct_missed(self) -> float:
        """Percentage of this category the other mask misses."""
        return 100.0 * self.bp_only / self.bp if self.bp else 0.0


@dataclass(frozen=True)
class UniqueRegion:
    """A merged region one mask has that the other does not touch, and why it was masked.

    `categories` and `names` come from the annotated BED when one was available, and are
    empty otherwise. A merged region can span more than one category.
    """

    chrom: str
    start: int
    end: int
    categories: tuple[str, ...] = ()
    names: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return self.end - self.start


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
    # Only filled when an annotated BED was available for that mask (see annotated_source).
    categories_a: list[CategoryComparison] = field(default_factory=list)
    categories_b: list[CategoryComparison] = field(default_factory=list)
    unique_a: list[UniqueRegion] = field(default_factory=list)
    unique_b: list[UniqueRegion] = field(default_factory=list)

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
    a: list[Region],
    b: list[Region],
    lengths: dict[str, int] | None = None,
    annotated_a: list[Region] | None = None,
    annotated_b: list[Region] | None = None,
) -> Comparison:
    """Compare two region sets. Each is merged first, so overlaps within one mask count once.

    `lengths` (record id -> bp, from a reference FASTA) is optional and only adds the
    percentages and the true figure axis; everything else is computed from the BEDs alone.

    `annotated_a` / `annotated_b` are the same masks with their categories: the unmerged
    regions from a mask.annotated.bed. A plain BED carries no categories, so without them
    the comparison still works but cannot say what kind of masking each difference is.
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
        categories_a=_by_category(annotated_a, mb),
        categories_b=_by_category(annotated_b, ma),
        unique_a=_unique_regions(only_a, annotated_a),
        unique_b=_unique_regions(only_b, annotated_b),
    )


def _by_category(annotated: list[Region] | None, other: list[Region]) -> list[CategoryComparison]:
    """One row per mask category: how much of it the other mask also covers.

    Each category is merged on its own first, so bases masked twice within a category count
    once, exactly as summary.tsv counts them.
    """
    if not annotated:
        return []
    rows = []
    for name in sorted({r.category for r in annotated if r.category}):
        regions = [r for r in annotated if r.category == name]
        rows.append(
            CategoryComparison(
                category=name,
                regions=len(merge(regions)),
                bp=union_bp(regions),
                bp_shared=intersect_bp(regions, other),
            )
        )
    # Biggest difference first: that is what the reader is looking for.
    return sorted(rows, key=lambda r: -r.bp_only)


def _unique_regions(
    only: list[Region], annotated: list[Region] | None
) -> list[UniqueRegion]:
    """Label each region the other mask misses with the categories and names behind it.

    Without an annotated BED the labels are empty and the region is still listed.
    """
    if not annotated:
        return [UniqueRegion(r.chrom, r.start, r.end) for r in only]
    by_chrom: dict[str, list[Region]] = defaultdict(list)
    for r in annotated:
        by_chrom[r.chrom].append(r)
    out = []
    for region in only:
        covering = [
            r
            for r in by_chrom.get(region.chrom, [])
            if r.start < region.end and r.end > region.start
        ]
        out.append(
            UniqueRegion(
                region.chrom,
                region.start,
                region.end,
                tuple(sorted({r.category for r in covering if r.category})),
                tuple(dict.fromkeys(r.name for r in covering if r.name)),
            )
        )
    return out


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
        past = sorted({r.chrom for r in ma + mb if r.chrom in lengths and r.end > lengths[r.chrom]})
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


# A plain mask.bed has an annotated twin written beside it by the same run.
NONE = "none"


def annotated_source(mask: Path, override: str | None) -> Path | None:
    """Which file to take mask categories from, or None to go without them.

    In order: `--annotated-a/-b` if given (the literal "none" switches it off), the mask
    itself when it is already an annotated BED, then a mask.annotated.bed sitting next to it,
    which is where `refmasker run` puts one.
    """
    if override == NONE:
        return None
    if override:
        path = Path(override)
        if not path.is_file():
            raise FileNotFoundError(f"No such annotated BED: {path}")
        if not is_annotated_bed(path):
            raise ValueError(
                f"{path} is not a refmasker mask.annotated.bed (its header does not start "
                "with #chrom, start, end, category)"
            )
        return path
    if is_annotated_bed(mask):
        return mask
    sibling = mask.with_suffix(".annotated.bed")
    return sibling if sibling.is_file() and is_annotated_bed(sibling) else None


def compare_files(
    a: Path,
    b: Path,
    rename: dict[str, str] | None = None,
    ref: Path | None = None,
    annotated_a: Path | None = None,
    annotated_b: Path | None = None,
) -> Comparison:
    """Read two BED files and compare them.

    `rename` maps old to new sequence names in both files, for masks that call the same
    chromosome by different names. `ref` is the reference FASTA the masks were built for,
    which supplies record lengths: without it the comparison is still complete, but there is
    no genome size to express the masks as a percentage of. `annotated_a` / `annotated_b`
    are mask.annotated.bed files (see annotated_source), which add the categories behind
    each difference.
    """

    def load(path: Path, annotated: Path | None) -> tuple[list[Region], list[Region] | None]:
        """The mask, and the same mask with categories when one is available."""
        labelled = read_annotated_bed(annotated) if annotated else None
        # Reading the annotated twin of a mask is enough on its own: merged, it is the mask.
        plain = labelled if (annotated and annotated == path) else read_bed(path)
        return plain, labelled

    ra, labelled_a = load(a, annotated_a)
    rb, labelled_b = load(b, annotated_b)
    if rename:
        ra, rb = _rename(ra, rename), _rename(rb, rename)
        labelled_a = _rename(labelled_a, rename) if labelled_a else None
        labelled_b = _rename(labelled_b, rename) if labelled_b else None
    lengths = {r.id: r.length for r in read_fasta(ref)} if ref else None
    return compare(ra, rb, lengths, labelled_a, labelled_b)


def _rename(regions: list[Region], rename: dict[str, str]) -> list[Region]:
    """Apply --rename to a region list, keeping every other field."""
    return [replace(r, chrom=rename.get(r.chrom, r.chrom)) for r in regions]


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
