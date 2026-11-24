# Mask categories

Each category is a kind of region that makes SNP calls against a reference unreliable. Four are on
by default (`repeats`, `is`, `prophage`, `rrna`). Choose others with `--categories`, e.g.
`--categories default,crispr`. `all` selects every category that needs no extra input.

Every region in `mask.annotated.bed` and the report records its category, the tool that found it,
and a short name. A base can be in more than one category; `mask.bed` is the union.

## repeats (default)

**What:** sequence present in two or more copies in the reference.

**Why:** short reads from one copy map equally well to the other copies. Real differences between
the copies then show up as SNPs, and true SNPs in one copy can be missed. snippy cannot tell which
copy a read came from, so **every copy** is masked, not only the extra ones.

**How:** the genome is aligned to itself with `nucmer --maxmatch --nosimplify` (MUMmer4), and each
sequence's alignment to itself is dropped. An alignment is kept if both copies are at least
`--repeat-min-len` bp and share at least `--repeat-min-id`% identity. Both copies become regions.

**Thresholds:** `--repeat-min-len` defaults to max(200 bp, `--read-length`). Copies shorter than a
read are spanned by reads anchored in unique flanking sequence, so they map correctly. The 200 bp
floor stays because repeats only slightly shorter than a read still attract mismapped reads.
`--repeat-min-id` defaults to 90%: below that, mismapped reads carry so many mismatches that
aligners usually reject them. The report includes a sensitivity table showing the % masked at
100–500 bp and 80–95% identity, so you can see how much the choice matters for your genome.

**Limitations:** gene families diverged below the identity threshold are not caught (for example,
many *M. tuberculosis* PE/PPE genes). Mask those by name with `--mask-genes-matching`. Tandem
repeats show up as overlapping self-alignments and are masked.

## is (default)

**What:** insertion sequences and transposons.

**Why:** IS elements are multi-copy and move, so their number and position differ between isolates
of the same ST. Reads from them map ambiguously.

**How:** [ISEScan](https://github.com/xiezhq/ISEScan) predicts IS elements from sequence. By
default only calls ISEScan marks *complete* are used; add `--isescan-partial` to include partial
calls, which catches more IS remnants but also more false positives (on H37Rv, one partial call
covers *dnaA*). When `--annotation` is given, genes whose product mentions a transposase or an IS
family, and features typed `mobile_element` / `mobile_genetic_element`, are added (source
`annotation`). This catches transposases ISEScan misses, and RefSeq (PGAP) annotation is usually
good at labelling them.

**Limitations:** larger composite transposons and integrative conjugative elements are only
partly covered: their IS ends are masked, but the cargo between them is not.

## prophage (default)

**What:** integrated phage genomes.

**Why:** prophages are gained and lost frequently, even between closely related isolates, and
carry SNPs that don't reflect the isolates' shared ancestry.

**How:** [geNomad](https://github.com/apcamargo/genomad) `end-to-end` by default, using its
`*_virus_summary.tsv`. Proviruses are masked at geNomad's coordinates. A whole record that geNomad
calls viral is masked entirely. geNomad needs its database: `--genomad-db DIR` or
`REFMASKER_GENOMAD_DB`. Download it once with `refmasker download-db DIR`; see
[the README](../README.md#the-genomad-database) for what that does and how to verify it.

`--phage-tool phispy` uses [PhiSpy](https://github.com/linsalrob/PhiSpy) instead. PhiSpy works
from an annotated GenBank file, so it needs `--annotation file.gbk`, and its sequence IDs must
match the FASTA.

**Limitations:** prophage boundaries are approximate. Use `--pad` if you want a margin. Degraded
prophage remnants may be missed.

## rrna (default)

**What:** ribosomal RNA operons.

**Why:** most bacteria carry several near-identical rRNA operons. Reads collapse across the copies,
so variant calls inside them are unreliable.

**How:** [barrnap](https://github.com/tseemann/barrnap) finds 16S, 23S and 5S genes. Genes no more
than `--rrna-operon-gap` bp apart (default 1000) are joined into one operon region, which also
masks the spacer between them (usually holding tRNA genes). `--rrna-operon-gap 0` masks only the
genes.

**Note:** a genome with a single rRNA operon (like *M. tuberculosis*) has nothing for reads to
collapse onto. It is still masked by default for consistency; leave out `rrna` if you prefer.

## lowcomplexity (opt-in)

**What:** low-complexity sequence and homopolymer runs.

**Why:** aligners and Illumina chemistry both struggle in these regions, giving spurious SNPs and
indels.

**How:** `dustmasker` (from BLAST+) with default settings, plus homopolymer runs of at least
`--homopolymer-min-len` bases (default 8).

**Limitations:** in GC-rich genomes, dustmasker flags many short intervals. On H37Rv it adds 771
intervals and about 1.8% of the genome, so look at the report before keeping it.

## crispr (opt-in)

**What:** CRISPR repeat-spacer arrays, found with [MinCED](https://github.com/ctSkennerton/minced).

**Why:** spacers are gained and lost quickly, and the direct repeats are repetitive by definition.

**Limitations:** MinCED can call short tandem repeats in GC-rich genes as CRISPR arrays (two of the
four calls on H37Rv are inside PE_PGRS genes).

## genes (opt-in)

**What:** annotated genes whose product, gene name or locus tag matches a regular expression
(case-insensitive), e.g. `--mask-genes-matching '\bPE\b|PPE|PE_PGRS'`. Using the option turns the
category on, and it needs `--annotation`.

**Why:** some gene families are known to be unreliable for SNP calling even where they aren't
detected as repeats, for example *M. tuberculosis* PE/PPE genes, or surface-protein families with
internal repeats.

## plasmid (opt-in)

**What:** whole plasmid records. `--mask-plasmids` turns it on. A record counts as a plasmid if its
FASTA header contains "plasmid" (as RefSeq headers do), or if you name it with
`--plasmid-ids id1,id2`.

**Why:** plasmids may be missing from some isolates or exchanged horizontally, which can distort a
chromosomal core-genome tree.

## user (opt-in)

**What:** any BED file you supply with `--extra-bed` (repeatable), for example regions from a
published mask or ones you've found problematic. Sequence names must match the FASTA IDs. Using
the option turns the category on.

## Not included

- **Recombination.** It depends on your set of isolates, not on the reference, so it can't be
  masked in advance. Run Gubbins or ClonalFrameML on snippy-core's alignment.
- **Genomic islands.** Island predictors such as IslandPath-DIMOB depend on annotation and are not
  actively maintained. The mobile parts of most islands (IS elements, prophages) are already
  covered. Use `--extra-bed` if you have island coordinates you trust.
- **Low-mappability regions from simulated reads.** The self-alignment approximates this. It may
  be added later.
