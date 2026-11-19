"""User-supplied regions from --extra-bed."""

from __future__ import annotations

from refmasker.detectors import Context, Result
from refmasker.regions import read_bed


def detect(ctx: Context) -> Result:
    result = Result()
    lengths = ctx.lengths
    for path in ctx.options.extra_bed or []:
        for r in read_bed(path, category="user", source="extra-bed"):
            if r.chrom not in lengths:
                raise ValueError(
                    f"{path}: sequence '{r.chrom}' is not in the reference FASTA "
                    f"(have: {', '.join(lengths)})"
                )
            if r.end > lengths[r.chrom]:
                result.notes.append(
                    f"{path}: region {r.chrom}:{r.start}-{r.end} runs past the end of the sequence; clipped."
                )
            result.regions.append(r)
    return result
