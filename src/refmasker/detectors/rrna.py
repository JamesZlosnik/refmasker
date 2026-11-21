"""rRNA operons: barrnap hits, merged into operons."""

from __future__ import annotations

import logging
from collections import defaultdict

from refmasker import external
from refmasker.detectors import Context, Result
from refmasker.regions import Region

log = logging.getLogger(__name__)


def parse_barrnap_gff(text: str) -> list[Region]:
    """Parse barrnap's GFF: one 'rRNA' line per gene.

    Names like '16S_rRNA' become '16S', and partial genes are marked '(partial)'.
    """
    regions = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) != 9 or cols[2] != "rRNA":
            continue
        attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
        name = attrs.get("Name", "rRNA").replace("_rRNA", "")
        if "partial" in attrs.get("product", ""):
            name += "(partial)"
        regions.append(Region(cols[0], int(cols[3]) - 1, int(cols[4]), "rrna", "barrnap", name))
    return regions


def group_operons(regions: list[Region], max_gap: int) -> list[Region]:
    """Join rRNA genes no more than `max_gap` bp apart into one operon region.

    The spacer between genes (often holding tRNAs) is covered too. max_gap=0 keeps
    the individual genes.
    """
    if max_gap <= 0:
        return list(regions)
    by_chrom: dict[str, list[Region]] = defaultdict(list)
    for r in regions:
        by_chrom[r.chrom].append(r)
    out = []
    for chrom, rs in by_chrom.items():
        rs.sort(key=lambda r: r.start)
        group = [rs[0]]
        for r in rs[1:]:
            # Gap between this gene and the furthest end reached by the current operon.
            if r.start - max(g.end for g in group) <= max_gap:
                group.append(r)
            else:
                out.append(_operon(chrom, group))
                group = [r]
        out.append(_operon(chrom, group))
    return out


def _operon(chrom: str, group: list[Region]) -> Region:
    """One region spanning a group of rRNA genes, e.g. 'operon 16S,23S,5S' (or a lone gene)."""
    names = ",".join(r.name for r in group)
    label = f"operon {names}" if len(group) > 1 else names
    return Region(
        chrom, min(r.start for r in group), max(r.end for r in group), "rrna", "barrnap", label
    )


def detect(ctx: Context) -> Result:
    """Run barrnap with its bacterial models, then join nearby rRNA genes into operons."""
    wd = ctx.workdir / "barrnap"
    wd.mkdir(parents=True, exist_ok=True)
    gff = wd / "barrnap.gff"
    external.run(
        ["barrnap", "--kingdom", "bac", "--threads", str(ctx.threads), "--quiet", str(ctx.fasta)],
        wd,
        "barrnap",
        stdout_path=gff,
    )
    genes = parse_barrnap_gff(gff.read_text())
    regions = group_operons(genes, ctx.options.rrna_operon_gap)
    log.info("rrna: %d rRNA genes in %d regions", len(genes), len(regions))
    return Result(regions=regions)
