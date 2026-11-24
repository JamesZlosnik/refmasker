# How refmasker works

A guide to what each part of refmasker does and how the parts fit together. For how to *use* it,
see the [README](../README.md) and [parameters](parameters.md); for what each output file contains,
see [outputs](outputs.md).

refmasker version 0.1.0.

## 1. What it does

refmasker takes a bacterial reference genome (FASTA, optionally with a GFF3 or GenBank annotation)
and writes a BED file of regions where SNP calls are unreliable: repeats, insertion sequences,
prophages, rRNA operons and, if asked, low-complexity sequence, CRISPR arrays, named genes,
plasmids and your own regions. You pass that BED to snippy (`snippy-core --mask mask.bed`) so
those regions are excluded from the core-genome alignment.

Each kind of region is a **category**. Each category has a **detector**: a small module that runs
one published tool (or reads the annotation) and turns its output into regions. refmasker does
not reimplement any detection algorithm; it runs the tools, parses their output, combines the
results and explains them.

## 2. The `refmasker run` workflow

Everything below happens in `pipeline.run()` (`src/refmasker/pipeline.py`).

```text
reference.fna (+ annotation)
        |
        v
 1. Prepare the output folder ......... refuse non-empty folders; --force only clears a
                                        previous refmasker output (work/, figures/)
 2. Decide threads .................... --threads, else SLURM/SGE allocation, else 1
 3. Resolve categories ................ "default", "all", or a list; check tools are on PATH
 4. Read the reference ................ validate FASTA; find plasmid records; warn on drafts and Ns
 5. Read the annotation (optional) .... GFF3 or GenBank -> features + overlap index
 6. Write work/reference.fna .......... same sequences, ID-only headers, so every tool
                                        reports the same sequence names
 7. Run each detector ................. one per category -> list of Regions
 8. Pad, clip, label genes, sort ...... --pad bp each side, clipped to sequence ends;
                                        each region lists the genes it overlaps
 9. Merge ............................. overlapping regions -> mask.bed
10. Summarise ......................... bp and % masked per category and record -> summary.tsv;
                                        repeat-threshold sensitivity grid
11. Record provenance ................. tool versions, geNomad DB version, checksums, parameters
12. Write reports ..................... figures, report.html, report.pdf, run_info.json
13. Clean up .......................... delete work/ unless --keep-work
```

Key points:

- **Coordinates.** Internally every region is BED-style: 0-based start, exclusive end. Each parser
  converts its tool's convention (most are 1-based inclusive; dustmasker is 0-based inclusive)
  as it reads. The HTML and PDF tables show 1-based starts for people; the BED files stay
  0-based.
- **Nothing is dropped silently.** Regions on sequences that are not in the FASTA are removed,
  but each removal is reported as a warning in the log and reports.
- **Categories overlap.** A region can be both a repeat and an IS element. Per-category figures
  count each base once within that category; the total counts each base once across all
  categories.

## 3. Commands

| Command | What it does | Code |
|---|---|---|
| `refmasker run REF -o DIR` | Builds the mask and reports (section 2) | `cli.cmd_run` -> `pipeline.run` |
| `refmasker compare A.bed B.bed` | Compares two masks: bp in both, Jaccard index, regions only in one, per-sequence breakdown. Prints a summary and writes a self-contained HTML report with two figures | `cli.cmd_compare` -> `compare.py`, `report/compare_html.py` |
| `refmasker check` | Lists each external tool, its version, and the geNomad database with the fingerprint recorded for it | `cli.cmd_check` |
| `refmasker download-db DIR` | Explains what it is about to do, runs `genomad download-database DIR` with geNomad's own output shown live, then SHA-256s every file of the database and writes `genomad_db.sha256sums` beside it | `cli.cmd_download_db` -> `database.py` |
| `refmasker download-db DIR --verify` | Re-checksums a downloaded database and reports files that changed, went missing or appeared (exit code 2 if it no longer matches) | `cli._verify_db` -> `database.verify` |

Exit codes: 0 success, 1 `check` found missing tools (or no command given), 2 an error such as
bad input, a missing tool or a refused output folder.

## 4. Repository layout

| Path | Purpose |
|---|---|
| `src/refmasker/` | The Python package |
| `src/refmasker/detectors/` | One module per mask category |
| `src/refmasker/report/` | Figures, HTML report, PDF report |
| `tests/` | pytest tests (fast tests use synthetic data; `slow` tests run the real tools) |
| `scripts/` | Shell scripts to download validation data and run the validation |
| `docs/` | User documentation, including this guide |
| `pyproject.toml` | Package metadata, dependencies, entry point, ruff and pytest settings |
| `environment.yml` | Pinned conda environment: Python, libraries and all external tools |
| `.github/workflows/` | CI: lint and fast tests on every push and pull request |

## 5. The Python package, file by file

### Entry points

| File | What it does |
|---|---|
| `__init__.py` | Holds the version number (`__version__`). |
| `__main__.py` | Lets you run `python -m refmasker`; calls `cli.main`. |
| `cli.py` | Builds the command-line parser (all options and their help text), sets up logging, and dispatches to `run`, `compare`, `check` or `download-db`, and prints what those last three report. Converts expected errors into a logged message and exit code 2 instead of a traceback. Also sets the default `--repeat-min-len` from `--read-length` (the larger of 200 bp and the read length). |

### The workflow and its inputs

| File | What it does |
|---|---|
| `pipeline.py` | The `run` workflow in section 2. Also: `resolve_categories` (turns `default,crispr` into an ordered list), `prepare_outdir` (the `--force` guard: only clears folders containing `run_info.json` or `refmasker.log`), `redactor` (for `--redact`: replaces input paths with file names in reports), `annotate_genes`, and `write_run_info` (writes `run_info.json`). |
| `seqio.py` | Reads FASTA (plain or gzipped) and checks it: no empty or duplicate IDs, no empty records. Writes the normalised FASTA used by the tools. `checksums` computes MD5 and SHA-256 of the reference in one pass. |
| `annotation.py` | Reads GFF3 (NCBI/PGAP, Prokka, Bakta, Ensembl) or GenBank (via Biopython) into `Feature` records. Works out gene names through GFF parent links, flags mobile elements (feature types such as `mobile_genetic_element`, or CDS products mentioning transposase/IS), and provides `FeatureIndex` for fast "which genes does this region overlap?" lookups. |
| `resources.py` | Decides the thread count without hogging a shared machine: an explicit `--threads` wins, then `SLURM_CPUS_PER_TASK`, `NSLOTS` or `PBS_NUM_PPN`, otherwise 1. Also caps OpenMP/BLAS thread pools to the same number. |
| `database.py` | The geNomad database: runs geNomad's downloader (`database.download`), checksums every file of the downloaded database and writes a `sha256sum`-format manifest beside the folder, and verifies a database against that manifest later. `fingerprint` reduces the whole set of checksums to one SHA-256, so two installs can be compared with a single line. The progress bar used while reading 1.4 GB lives here too. |
| `external.py` | Runs external tools. Commands are always argument lists (never through a shell). Each tool's stderr goes to `work/<tool>/<tool>.log`; on failure the last 20 lines are shown. `run_attached` instead lets a tool write straight to the terminal, so the geNomad database download shows its own progress. Also finds missing tools and reads tool versions. |

### Regions and numbers

| File | What it does |
|---|---|
| `regions.py` | The `Region` record (sequence, start, end, category, source tool, name, genes overlapped) and the interval maths: `merge`, `pad_and_clip`, `union_bp`, and a two-pointer sweep (`overlaps`) shared by `intersect_bp`, `intersect` and `subtract`. Reads BED (`--extra-bed`, `compare`) and writes `mask.bed` (3 columns, no header, as snippy needs) and `mask.annotated.bed`. |
| `summary.py` | Builds `summary.tsv`: regions, bp and % masked per category and per record, plus an overlap-aware total. Also the repeat sensitivity grid: how much would be masked at minimum lengths 100 to 500 bp and identities 80 to 95%. |
| `compare.py` | `refmasker compare`: merges each BED, then reports bp in each, bp in both, fraction of one covered by the other, Jaccard index, a row per sequence, and the regions found in only one file (optionally written as BED files). Keeps the per-base difference (`diff_a`/`diff_b`) apart from whole regions the other mask never touches (`only_a`/`only_b`). `--rename` fixes sequence-name mismatches; `--ref` adds the record lengths, which the BEDs do not carry, for percentages and a true-to-scale figure. Warns when the two masks share no sequence names, or when regions run past the end of the reference. |

### Detectors (`src/refmasker/detectors/`)

`detectors/__init__.py` defines the shared pieces: `Context` (the normalised FASTA, records,
work folder, threads, options and annotation that every detector receives), `Result` (the regions
plus notes for the report), and the **category registry**, which lists each category's name,
label, whether it is on by default, the tools it needs, and its `detect` function.

Each detector keeps the tool's output **parser as a separate pure function**, so parsing is
unit-tested on captured output without the tool installed.

| Category | Default | Module | Tool and command | What becomes a region |
|---|---|---|---|---|
| `repeats` | yes | `repeats.py` | `nucmer --maxmatch --nosimplify` (genome against itself), then `show-coords` | Both copies of every alignment at or above `--repeat-min-len` and `--repeat-min-id`; the trivial self-hit is dropped. Overlapping copies are merged into one region per locus, named by how many other copies it matches. |
| `is` | yes | `insertion_sequences.py` | `isescan.py` | ISEScan's complete IS calls (partial calls with `--isescan-partial`), plus transposase and mobile-element features from the annotation. |
| `prophage` | yes | `prophages.py` | `genomad end-to-end` (default) or `PhiSpy.py` | geNomad proviruses (`contig\|provirus_START_END`), or whole sequences geNomad calls viral; or PhiSpy's prophage coordinates. geNomad needs its database; PhiSpy needs a GenBank annotation. |
| `rrna` | yes | `rrna.py` | `barrnap --kingdom bac` | rRNA genes, joined into one operon region when within `--rrna-operon-gap` bp (default 1000), which also covers the spacer between them. |
| `lowcomplexity` | no | `lowcomplexity.py` | `dustmasker -outfmt interval` | Low-complexity intervals, plus homopolymer runs of at least `--homopolymer-min-len` (found in Python with a regular expression). |
| `crispr` | no | `crispr.py` | `minced` | CRISPR repeat-spacer arrays. |
| `genes` | no | `genes.py` | none (annotation) | Annotated genes whose product, gene name or locus tag matches `--mask-genes-matching` (a case-insensitive regular expression), e.g. PE/PPE genes in *M. tuberculosis*. |
| `plasmid` | no | `plasmids.py` | none (FASTA) | Whole records whose header mentions "plasmid", or named with `--plasmid-ids`. |
| `user` | no | `user.py` | none (BED) | Regions from `--extra-bed` files; rejects sequence names not in the FASTA and notes regions that run past the end. |

### Reports (`src/refmasker/report/`)

| File | What it does |
|---|---|
| `__init__.py` | Parses `--formats` and calls the writers below in order. |
| `plots.py` | Linear and circular genome figures with matplotlib (one track per category, plus the merged mask), using the Okabe-Ito colour-blind-safe palette. Saved as PNG and SVG under `figures/`. Also the two `compare` figures (a shared-against-unique bar, and a track per mask along each sequence), which are rendered to a string and inlined rather than written to disk. |
| `compare_html.py` | The `refmasker compare` report: tiles, the two comparison figures inlined as SVG, a per-sequence table and sortable tables of the regions each mask has that the other misses. Reuses the page styles and table sorter from `html.py`. |
| `html.py` | A self-contained `report.html` (figures inlined as SVG, works offline): summary tiles, warnings, figures, tables, the sensitivity grid, a sortable and searchable table of every region, how to use the mask with snippy, and run details. All text from inputs is HTML-escaped. |
| `pdf.py` | `report.pdf` with the same content, drawn with matplotlib. |
| `text.py` | Plain-language explanations of why each category is masked, and the snippy usage text, shared by both reports. |

## 6. Scripts (`scripts/`)

Both scripts use `set -euo pipefail` and quote their variables.

| Script | What it does |
|---|---|
| `fetch_test_refs.sh [OUTDIR]` | Downloads three public validation references (FASTA + GFF3) from NCBI into `test_data/`: *M. tuberculosis* H37Rv (NC_000962.3), *S. aureus* USA300_FPR3757 (NC_007793.1 + three plasmids) and *P. aeruginosa* PAO1 (NC_002516.2). Also downloads snippy's published H37Rv mask, from a pinned snippy commit, as the benchmark. Pauses between NCBI requests to stay under its rate limit. |
| `validate.sh [DATA] [OUT] [THREADS]` | Runs refmasker on the three references, plus a second H37Rv run that also masks PE/PPE genes by name. Compares both H37Rv masks with snippy's using `refmasker compare` and prints the fraction masked for each run. Needs the conda environment and the geNomad database. |

## 7. Configuration and CI

| File | What it does |
|---|---|
| `pyproject.toml` | Package metadata; Python 3.11 or later; runtime dependencies (biopython, matplotlib); the `refmasker` command; ruff lint rules; the pytest `slow` marker. |
| `environment.yml` | The conda environment, with every package pinned to an exact version (Python 3.12, libraries, nucmer, ISEScan, geNomad, PhiSpy, barrnap, BLAST, MinCED, pytest, ruff). It does not install refmasker itself: run `pip install --no-deps .` afterwards. |
| `.github/workflows/lint.yml` | `ruff check` and `ruff format --check` on every push to main and every pull request. |
| `.github/workflows/tests.yml` | Fast tests (`pytest -m "not slow"`) on Python 3.11 and 3.13. Actions are pinned to commit SHAs and the job token is read-only. |
| `.gitignore` | Keeps build files, caches, `test_data/`, `runs/` and mask output folders out of git. |

## 8. Tests (`tests/`)

| File | What it covers |
|---|---|
| `conftest.py` | The synthetic test genome: a 30 kb chromosome with a 1 kb repeat planted twice and a homopolymer, plus a 5 kb plasmid (fixed random seeds, so it is the same every run). Also a gzipped copy. |
| `test_seqio.py` | FASTA reading (plain and gzipped), error cases, and that written FASTA has ID-only headers. |
| `test_annotation.py` | NCBI and Ensembl GFF3, GenBank, and the overlap index. |
| `test_regions.py` | Merging, padding and clipping, overlap and intersection maths, BED reading and writing. |
| `test_summary.py` | Overlapping categories counted once in the total; the sensitivity grid marks the chosen setting. |
| `test_resources.py` | Thread count: default 1, scheduler allocation, explicit value wins, zero rejected. |
| `test_parsers.py` | Every tool's output parser, on captured output: show-coords, ISEScan, geNomad (and its database version), PhiSpy, barrnap and operon grouping, dustmasker, homopolymers, MinCED, plasmid detection. |
| `test_cli.py` | End-to-end runs with the categories that need no external tools: outputs exist and are correct, the output-folder and `--force` rules, `--redact`, missing tools, `--formats none`, and `compare`. |
| `test_tools_slow.py` | Marked `slow`: runs the real tools on the synthetic genome (the planted repeat must be found; the other tools must run). Skips tools that are not installed. |

Run them with `pytest -m "not slow"` (seconds, no tools needed) and `pytest -m slow` (needs the
conda environment).

## 9. Safety on shared systems

- **No shell.** External tools are run with argument lists, so file names and options can never
  be interpreted as shell commands.
- **Bounded resources.** One thread unless told otherwise or given an allocation by the scheduler;
  numerical libraries are capped to the same count.
- **Output folder guard.** refmasker never writes into a non-empty folder without `--force`, and
  `--force` only clears folders that a previous refmasker run created.
- **Temporary files** stay inside the output folder (`work/`), never in a shared `/tmp`, and are
  deleted at the end unless `--keep-work`.
- **Reports and identifiers.** `--redact` keeps full paths and FASTA header descriptions out of
  the HTML, PDF and `run_info.json`. Record IDs are always kept (the mask needs them), and
  `refmasker.log` keeps full paths for debugging.
- **Provenance.** `run_info.json` records the refmasker and tool versions, the geNomad database
  version, every parameter, and the reference's MD5 and SHA-256.
- **Permissions** follow your umask; see the README section on shared HPC systems.

## 10. Adding a new category

1. Create `src/refmasker/detectors/<name>.py` with a pure `parse_...` function for the tool's
   output and a `detect(ctx) -> Result` function that runs the tool with `external.run` and
   returns `Region`s tagged with the new category.
2. Register it in `_registry()` in `detectors/__init__.py` (name, label, default on or off,
   tools needed, `detect`).
3. Add a colour in `report/plots.py` (`COLOURS`) and an explanation in `report/text.py` (`WHY`).
4. Add the tool, pinned, to `environment.yml` and its version command to
   `external.VERSION_COMMANDS`.
5. Add a parser test with captured output to `tests/test_parsers.py`, and a `slow` test if
   useful.
6. Document it in `docs/categories.md` and `CHANGELOG.md`.

## 11. Further reading

- [categories.md](categories.md): what each category masks and why.
- [parameters.md](parameters.md): every option.
- [outputs.md](outputs.md): every output file.
- [methods.md](methods.md): a methods paragraph and references for papers.
- [validation.md](validation.md): how the masks compare with snippy's H37Rv mask.
- [faq.md](faq.md): common questions.
