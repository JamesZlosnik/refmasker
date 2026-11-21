"""User-supplied regions from --extra-bed."""

from __future__ import annotations

from refmasker.detectors import Context, Result
from refmasker.regions import read_bed


def detect(ctx: Context) -> Result:
    """Load the BED files given with --extra-bed.

    A sequence name not in the FASTA is an error. Regions running past the end of a sequence
    are kept with a note; the pipeline clips them to the sequence length.
    """
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
