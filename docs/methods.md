# Methods text and citations

Edit the numbers to match your run. `run_info.json` and the report's run details list the version
of every tool.

## Methods paragraph

> Regions of the reference genome likely to produce unreliable SNP calls were masked with
> refmasker v0.1.0. Repeats were identified by aligning the reference to itself with nucmer
> (MUMmer v4.0.1; `--maxmatch --nosimplify`) [1]. All copies of repeats of at least 200 bp with at
> least 90% identity between copies were masked. Insertion sequences were predicted with ISEScan
> v1.7.3 [2] (complete elements only), supplemented with transposase genes from the RefSeq
> annotation. Prophages were predicted with geNomad v1.x [3]. rRNA genes were found with barrnap
> v0.9 [4], and genes within 1 kb of each other were joined into operons. In total, X bp (Y%) of the
> reference was masked. The mask was supplied to snippy-core v4.x (`--mask`) [5] when building the
> core-genome SNP alignment.

Add sentences for any opt-in categories you used:

- Low-complexity regions were identified with dustmasker (BLAST+ v2.x) [6], together with
  homopolymer runs of at least 8 bp.
- CRISPR arrays were identified with MinCED v0.4.2 [7, 8].
- Genes annotated as \<family\> were masked.
- Plasmid sequences were masked in full.

## References

1. Marçais G, Delcher AL, Phillippy AM, Coston R, Salzberg SL, Zimin A. MUMmer4: a fast and
   versatile genome alignment system. *PLoS Comput Biol.* 2018;14(1):e1005944.
2. Xie Z, Tang H. ISEScan: automated identification of insertion sequence elements in prokaryotic
   genomes. *Bioinformatics.* 2017;33(21):3340–3347.
3. Camargo AP, Roux S, Schulz F, et al. Identification of mobile genetic elements with geNomad.
   *Nat Biotechnol.* 2024;42:1303–1312.
4. Seemann T. barrnap: Bacterial ribosomal RNA predictor. https://github.com/tseemann/barrnap
5. Seemann T. snippy: fast bacterial variant calling from NGS reads.
   https://github.com/tseemann/snippy
6. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. A fast and symmetric DUST implementation to mask
   low-complexity DNA sequences. *J Comput Biol.* 2006;13(5):1028–1040.
7. Skennerton CT. MinCED: Mining CRISPRs in Environmental Datasets.
   https://github.com/ctSkennerton/minced
8. Bland C, Ramsey TL, Sabree F, et al. CRISPR Recognition Tool (CRT): a tool for automatic
   detection of clustered regularly interspaced palindromic repeats. *BMC Bioinformatics.*
   2007;8:209.

If you use PhiSpy instead of geNomad:

- Akhter S, Aziz RK, Edwards RA. PhiSpy: a novel algorithm for finding prophages in bacterial
  genomes that combines similarity- and composition-based strategies. *Nucleic Acids Res.*
  2012;40(16):e126.
