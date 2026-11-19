# Changelog

## 0.1.0 (unreleased)

First version.

- `refmasker run`: BED mask for snippy from a reference FASTA, with optional GFF3/GenBank annotation.
- Default categories: repeats (nucmer self-alignment), IS elements (ISEScan + annotation),
  prophages (geNomad or PhiSpy), rRNA operons (barrnap).
- Opt-in categories: low complexity (dustmasker + homopolymers), CRISPR (MinCED), named genes
  (`--mask-genes-matching`), plasmids, user BED.
- Repeat length threshold follows `--read-length`.
- Reports: HTML, PDF, PNG/SVG figures (linear and circular), `summary.tsv`, `run_info.json`.
- `refmasker compare`, `refmasker check`, `refmasker download-db`.
- Thread use defaults to 1, or the SLURM/SGE allocation.
