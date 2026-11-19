"""Compare two BED masks (e.g. refmasker vs a published mask)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from refmasker.regions import Region, intersect_bp, merge, read_bed


@dataclass
class Comparison:
    bp_a: int
    bp_b: int
    bp_both: int
    regions_a: int
    regions_b: int
    only_a: list[Region]
    only_b: list[Region]

    @property
    def jaccard(self) -> float:
        union = self.bp_a + self.bp_b - self.bp_both
        return self.bp_both / union if union else 1.0

    @property
    def frac_b_covered(self) -> float:
        return self.bp_both / self.bp_b if self.bp_b else 0.0

    @property
    def frac_a_in_b(self) -> float:
        return self.bp_both / self.bp_a if self.bp_a else 0.0


def _as_regions(merged) -> list[Region]:
    return [Region(m.chrom, m.start, m.end) for m in merged]


def compare(a: list[Region], b: list[Region]) -> Comparison:
    ma, mb = _as_regions(merge(a)), _as_regions(merge(b))
    both = intersect_bp(ma, mb)
    only_a = [r for r in ma if intersect_bp([r], mb) == 0]
    only_b = [r for r in mb if intersect_bp([r], ma) == 0]
    return Comparison(
        sum(r.length for r in ma),
        sum(r.length for r in mb),
        both,
        len(ma),
        len(mb),
        only_a,
        only_b,
    )


def compare_files(a: Path, b: Path, rename: dict[str, str] | None = None) -> Comparison:
    ra, rb = read_bed(a), read_bed(b)
    if rename:
        ra = [Region(rename.get(r.chrom, r.chrom), r.start, r.end) for r in ra]
        rb = [Region(rename.get(r.chrom, r.chrom), r.start, r.end) for r in rb]
    return compare(ra, rb)


def format_report(c: Comparison, name_a: str, name_b: str) -> str:
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
