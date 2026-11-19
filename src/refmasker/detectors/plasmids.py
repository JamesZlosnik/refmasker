"""Plasmids: whole records flagged as plasmids."""

from __future__ import annotations

from refmasker.detectors import Context, Result
from refmasker.regions import Region
from refmasker.seqio import Record


def plasmid_ids(records: list[Record], named: list[str] | None = None) -> list[str]:
    """Records whose FASTA description mentions 'plasmid', plus any named explicitly."""
    named_set = set(named or [])
    unknown = named_set - {r.id for r in records}
    if unknown:
        raise ValueError(
            f"--plasmid-ids names records not in the FASTA: {', '.join(sorted(unknown))}"
        )
    return [r.id for r in records if r.id in named_set or "plasmid" in r.description.lower()]


def detect(ctx: Context) -> Result:
    ids = plasmid_ids(ctx.records, ctx.options.plasmid_ids)
    lengths = ctx.lengths
    regions = [Region(i, 0, lengths[i], "plasmid", "fasta", "whole plasmid") for i in ids]
    notes = (
        [] if ids else ["No plasmid records found (header containing 'plasmid' or --plasmid-ids)."]
    )
    return Result(regions=regions, notes=notes)
