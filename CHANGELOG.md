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
- `refmasker download-db` now says what it is about to do (destination, download size, free
  space, which geNomad it found) and shows geNomad's downloader output live instead of hiding it
  in a log file, so the source URL and its progress are visible. Afterwards it records a SHA-256
  for every file of the database in a `sha256sum`-format manifest beside the folder, with a
  single fingerprint over the whole set, and prints it. It refuses to overwrite an existing
  database.
- New `refmasker download-db DIR --verify` re-checksums a database against that manifest and
  reports anything changed, missing or added; `refmasker check` shows the recorded fingerprint.

- `refmasker run`: BED mask for snippy from a reference FASTA, with optional GFF3/GenBank annotation.
- Default categories: repeats (nucmer self-alignment), IS elements (ISEScan + annotation),
  prophages (geNomad or PhiSpy), rRNA operons (barrnap).
- Opt-in categories: low complexity (dustmasker + homopolymers), CRISPR (MinCED), named genes
  (`--mask-genes-matching`), plasmids, user BED.
- Repeat length threshold follows `--read-length`.
- Reports: HTML, PDF, PNG/SVG figures (linear and circular), `summary.tsv`, `run_info.json`.
- `refmasker compare`, `refmasker check`, `refmasker download-db`.
- `refmasker compare` writes a self-contained HTML report by default: headline numbers, a bar of
  shared against unique bases, a track per mask along each sequence, a per-sequence table, and
  sortable tables of the regions each mask has that the other misses entirely. Named after the
  two inputs, or after `--out-prefix`; `--formats none` turns it off and `--formats html,png`
  also saves the figures. New `--ref` supplies the sequence lengths the BEDs do not carry, for
  percentages of the genome and a true-to-scale figure. The comparison also reports per-sequence
  numbers and warns when the two masks share no sequence names.
- The comparison report also says what kind of masking each difference is, when the categories
  are available: a table and a figure per mask category (bp, bp the other mask also covers, bp
  only this mask has, and the percentage the other misses), and the category and name on every
  row of the unique-region tables, filterable by category. A three-column `mask.bed` carries no
  categories, so refmasker reads them from a `mask.annotated.bed`: the mask itself when it is
  one, the file beside it, or `--annotated-a` / `--annotated-b` (`none` switches it off). The
  file used is named in the log and in the report.
- Thread use defaults to 1, or the SLURM/SGE allocation.
