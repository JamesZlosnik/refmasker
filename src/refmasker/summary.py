"""How much of the genome each category masks."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from refmasker.regions import Region, merge, union_bp
from refmasker.seqio import Record

ALL = "ALL"
TOTAL = "total"


@dataclass(frozen=True)
class SummaryRow:
    record: str
    category: str
    regions: int
    bp: int
    length: int

    @property
    def pct(self) -> float:
        return 100.0 * self.bp / self.length if self.length else 0.0


def summarise(
    regions: list[Region], records: list[Record], categories: list[str]
) -> list[SummaryRow]:
    """Rows per record and genome-wide (record 'ALL'), per category plus an overlap-aware total."""
    rows = []
    genome_len = sum(r.length for r in records)
    scopes = [(ALL, genome_len, None)] + [(r.id, r.length, r.id) for r in records]
    for name, length, chrom in scopes:
        scoped = [r for r in regions if chrom is None or r.chrom == chrom]
        for cat in categories:
            rs = [r for r in scoped if r.category == cat]
            rows.append(SummaryRow(name, cat, len(merge(rs)), union_bp(rs), length))
        rows.append(SummaryRow(name, TOTAL, len(merge(scoped)), union_bp(scoped), length))
    return rows


def write_summary_tsv(path: str | Path, rows: list[SummaryRow]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["record", "category", "regions", "bp_masked", "record_length", "pct_masked"])
        for r in rows:
            w.writerow([r.record, r.category, r.regions, r.bp, r.length, f"{r.pct:.3f}"])


@dataclass(frozen=True)
class SensitivityRow:
    min_len: int
    min_id: float
    regions: int
    bp: int
    pct: float
    chosen: bool


def repeat_sensitivity(
    hits, genome_len: int, chosen_len: int, chosen_id: float, lengths, identities
) -> list[SensitivityRow]:
    from refmasker.detectors.repeats import hits_to_regions

    grid = sorted({(L, float(i)) for L in lengths for i in identities} | {(chosen_len, chosen_id)})
    rows = []
    for L, i in grid:
        rs = hits_to_regions(hits, L, i)
        bp = union_bp(rs)
        rows.append(
            SensitivityRow(
                L, i, len(merge(rs)), bp, 100.0 * bp / genome_len, (L, i) == (chosen_len, chosen_id)
            )
        )
    return rows
