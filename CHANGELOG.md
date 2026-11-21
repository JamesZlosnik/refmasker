# Changelog

## 0.1.0 (unreleased)

First version.

Security and HPC hardening:

- `--force` refuses to clear a non-empty folder that is not a previous refmasker output, so it can
  no longer delete an unrelated `work/` (e.g. Nextflow's).
- New `--redact`: reports show file names instead of full paths and omit FASTA header descriptions.
- `run_info.json` and the reports record the reference's SHA-256 alongside MD5, and the geNomad
  database version; `refmasker check` shows the database version.
- Python 3.11 or later (3.10 reached end of life in October 2026); the conda environment uses 3.12.
- `environment.yml` pins every package and no longer installs refmasker in editable mode; install it
  with `pip install --no-deps .`.
- CI: GitHub Actions pinned to commit SHAs, read-only token, ruff pinned; tests on 3.11 and 3.13.
- `scripts/fetch_test_refs.sh` fetches snippy's TB mask from a pinned commit.
- The output folder must belong to you and not be writable by every user; `--force` also refuses a
  folder containing symbolic links. Output files are written to a temporary file and moved into
  place, so a planted link is replaced rather than followed.
- New output folders get no access for other users; `work/` (which holds a copy of the reference)
  is private and is deleted after failed runs too, keeping the tool logs in `logs/`.
- Thread detection also reads `SLURM_CPUS_ON_NODE` and PBS Pro's `NCPUS`, warns when `--threads`
  is more than the job was given, and caps TensorFlow and Java thread pools.
- New `--max-genome-mb` (default 50) refuses oversized references while reading them; nucmer's
  alignments are streamed instead of read into memory at once.
- A `$` in a record ID no longer breaks the figures (matplotlib maths parsing is off).
- `run_info.json` and `refmasker.log` record each tool's full path; refmasker warns if a tool or
  the geNomad database can be modified by any user.
- CI installs dependencies from `.github/constraints.txt`, matching `environment.yml`.

- `refmasker run`: BED mask for snippy from a reference FASTA, with optional GFF3/GenBank annotation.
- Default categories: repeats (nucmer self-alignment), IS elements (ISEScan + annotation),
  prophages (geNomad or PhiSpy), rRNA operons (barrnap).
- Opt-in categories: low complexity (dustmasker + homopolymers), CRISPR (MinCED), named genes
  (`--mask-genes-matching`), plasmids, user BED.
- Repeat length threshold follows `--read-length`.
- Reports: HTML, PDF, PNG/SVG figures (linear and circular), `summary.tsv`, `run_info.json`.
- `refmasker compare`, `refmasker check`, `refmasker download-db`.
- Thread use defaults to 1, or the SLURM/SGE allocation.
