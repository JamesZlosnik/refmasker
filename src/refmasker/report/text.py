"""Plain-language explanations shared by the HTML and PDF reports."""

# One-paragraph reason per category, shown in the "Why each category is masked" sections.
WHY = {
    "repeats": (
        "Sequence present in two or more copies at or above the length and identity thresholds. "
        "Reads from one copy can map to another, so differences between copies show up as false "
        "SNPs. Every copy is masked because snippy cannot tell which copy a read came from."
    ),
    "is": (
        "Insertion sequences and transposons. They are multi-copy and mobile, so their number and "
        "position differ between isolates, and reads from them map ambiguously."
    ),
    "prophage": (
        "Integrated phage. Prophages are gained and lost frequently and can differ strongly between "
        "isolates of the same ST, adding SNPs that do not reflect vertical descent."
    ),
    "rrna": (
        "Ribosomal RNA operons. Bacteria carry several near-identical copies, so reads collapse "
        "across copies. Genes within the operon gap are joined, which also masks the spacer."
    ),
    "lowcomplexity": (
        "Low-complexity sequence (dustmasker) and homopolymer runs, where alignment and "
        "sequencing errors produce spurious variants and indels."
    ),
    "crispr": (
        "CRISPR repeat-spacer arrays. Spacer content changes quickly between isolates and the "
        "repeats are, by definition, repetitive."
    ),
    "plasmid": (
        "Whole plasmid records. Plasmids are often absent from some isolates or moved horizontally, "
        "and can distort a chromosomal core-genome phylogeny."
    ),
    "genes": (
        "Annotated genes whose name or product matched --mask-genes-matching, for gene families "
        "known to be unreliable (for example PE/PPE genes in M. tuberculosis)."
    ),
    "user": "Regions supplied with --extra-bed.",
}

# Example command shown in the reports for using the mask with snippy-core.
SNIPPY_USAGE = (
    "snippy-core --ref reference.fna --mask mask.bed --prefix core snippy_dir1 snippy_dir2 ..."
)


def fmt_bp(n: int) -> str:
    """Format a base-pair count with thousands separators, e.g. 4411532 -> '4,411,532'."""
    return f"{n:,}"
