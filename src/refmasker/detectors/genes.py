"""Genes chosen by name: annotated genes matching --mask-genes-matching."""

from __future__ import annotations

import re

from refmasker.detectors import Context, Result
from refmasker.regions import Region


def detect(ctx: Context) -> Result:
    pattern = ctx.options.mask_genes_matching
    if ctx.features is None:
        raise ValueError("--mask-genes-matching needs an --annotation file")
    try:
        rx = re.compile(pattern, re.IGNORECASE)
    except re.error as exc:
        raise ValueError(f"--mask-genes-matching: invalid regular expression: {exc}") from exc
    regions = []
    for f in ctx.features:
        if f.ftype not in {"CDS", "pseudogene", "rRNA", "tRNA", "ncRNA"}:
            continue
        if any(rx.search(x) for x in (f.product, f.gene, f.locus_tag) if x):
            name = " ".join(x for x in (f.locus_tag, f.product) if x)
            regions.append(Region(f.seqid, f.start, f.end, "genes", "annotation", name))
    notes = [] if regions else [f"No annotated genes matched '{pattern}'."]
    return Result(regions=regions, notes=notes)
