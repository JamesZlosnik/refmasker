# refmasker

**Build a BED mask of problematic regions in a bacterial reference genome, for core-genome SNP
analysis with [snippy](https://github.com/tseemann/snippy).**

snippy ships a hand-curated mask for *M. tuberculosis* H37Rv, but not for the per-ST references used
for most other species. refmasker builds one from a reference FASTA. It finds repeats, insertion
sequences, prophages and rRNA operons, and optionally other region types. It writes a BED file for
`snippy --mask` / `snippy-core --mask` and a report showing what was masked, why, and how much of
the genome that is.

```
refmasker run reference.fna --annotation reference.gff -o reference_mask
snippy-core --ref reference.fna --mask reference_mask/mask.bed --prefix core isolate1 isolate2 ...
```

## What gets masked

| Category | Default | Found by | Why |
| --- | --- | --- | --- |
| `repeats` | on | nucmer self-alignment | Reads map to the wrong copy, giving false SNPs |
| `is` | on | ISEScan, plus transposase genes from the annotation | Multi-copy, mobile |
| `prophage` | on | geNomad (or PhiSpy) | Gained and lost between isolates |
| `rrna` | on | barrnap, joined into operons | Several near-identical copies |
| `lowcomplexity` | off | dustmasker, homopolymer runs | Alignment and sequencing artefacts |
| `crispr` | off | MinCED | Spacers change between isolates |
| `genes` | off | `--mask-genes-matching REGEX` on the annotation | Gene families you know are unreliable (e.g. PE/PPE) |
| `plasmid` | off | `--mask-plasmids`, `--plasmid-ids` | Plasmids are often absent or mobile |
| `user` | off | `--extra-bed` | Your own regions |

See [docs/categories.md](docs/categories.md) for how each one is detected and its limitations.
Recombination is **not** handled here because it depends on your isolates, not the reference: run
Gubbins or ClonalFrameML on the core alignment as usual.

## Install

refmasker is a Python package; the detection tools come from bioconda.

```bash
git clone https://github.com/<you>/refmasker.git
cd refmasker
conda env create -f environment.yml      # or mamba: installs the pinned tools, not refmasker
conda activate refmasker
pip install --no-deps .                  # installs a copy of refmasker into the environment
refmasker check                          # lists each tool and its version
```

`pip install --no-deps .` copies refmasker into the environment, so the environment keeps working if
the checkout is moved, edited or deleted, and `--no-deps` stops pip pulling anything from PyPI (conda
has already installed the dependencies). Avoid editable installs (`pip install -e`) for shared
environments: anyone who can edit the checkout would change the code every user of the environment
runs.

### Shared HPC systems

- refmasker uses 1 thread unless told otherwise, or the allocation from SLURM (`SLURM_CPUS_PER_TASK`)
  or SGE (`NSLOTS`). Run it inside a job, not on a login node.
- Output files get your umask's permissions; on many clusters that makes them readable by other
  users. For sensitive references, run with a private umask:

  ```bash
  umask 077      # only you; or umask 027 to share with your group
  refmasker run ref.fna -o ref_mask
  ```

- To share reports without revealing your folder layout or FASTA header text, add `--redact`.
- If your home folder is read-only or over quota on compute nodes, matplotlib may warn that it cannot
  write its cache; set `export MPLCONFIGDIR="$TMPDIR/mpl"` in the job script.

The default prophage caller, geNomad, needs its database (about 1.4 GB, downloaded once):

```bash
refmasker download-db /path/to/dbs
export REFMASKER_GENOMAD_DB=/path/to/dbs/genomad_db     # add to ~/.bashrc
```

If you don't want the download, use `--phage-tool phispy` (needs a GenBank `--annotation`) or leave
prophages out with `--categories repeats,is,rrna`.

The conda environment targets Linux (x86-64). ISEScan and geNomad may not have macOS ARM builds;
run refmasker on your Linux server.

## Quickstart

1. Get the reference you use for the ST, plus its annotation. For a RefSeq genome, download the
   FASTA and GFF3 together, for example with NCBI Datasets:

   ```bash
   datasets download genome accession <GCF_accession> --include genome,gff3
   ```

2. Build the mask (here with 4 threads on a compute node):

   ```bash
   refmasker run ST8_ref.fna --annotation ST8_ref.gff -o ST8_mask --threads 4
   ```

3. Read `ST8_mask/report.html`. Check how much is masked and look over the biggest regions.

4. Use the mask with snippy:

   ```bash
   snippy-core --ref ST8_ref.fna --mask ST8_mask/mask.bed --prefix core snps/*
   ```

   The sequence names in `mask.bed` are the FASTA IDs, so give snippy the **same** FASTA.
   For nullarbor, pass the mask through to snippy-core (see [docs/faq.md](docs/faq.md)).

## Outputs

| File | What it is |
| --- | --- |
| `mask.bed` | Merged regions, 3 columns, no header: what snippy needs |
| `mask.annotated.bed` | Every region before merging, with category, tool, name and overlapping genes |
| `summary.tsv` | bp and % masked per category and per record, plus an overlap-aware total |
| `report.html` | Self-contained report: plots, tables, sortable region list, run details |
| `report.pdf` | The same report as a PDF |
| `figures/linear.{svg,png}`, `figures/circular.{svg,png}` | Genome plots, one track per category |
| `run_info.json` | Parameters, tool versions, input checksum, headline numbers |
| `refmasker.log` | Full log, including every command run |

Choose report formats with `--formats html,pdf,png` (default: all three) or `--formats none`.
More detail in [docs/outputs.md](docs/outputs.md).

## Key options

```
--read-length 150        read length the mask is for; the minimum repeat length is
                         max(200, read length), so use --read-length 300 for 2x300 data
--repeat-min-id 90       minimum % identity between repeat copies
--categories default     e.g. default,crispr,lowcomplexity  or  repeats,is,rrna
--pad 0                  extend every region by N bp on each side
--threads                default 1, or the SLURM/SGE allocation
```

Full list: `refmasker run --help` and [docs/parameters.md](docs/parameters.md).

## Running on a shared server

refmasker uses **one thread by default**. Inside a SLURM or SGE job it uses the cores the job was
given (`SLURM_CPUS_PER_TASK` / `NSLOTS`). It never grabs every core on the machine, so it is safe to
run on a head node. Ask for more with `--threads N`. ISEScan is the slowest step: about
7 minutes for H37Rv on 4 threads, 12 minutes for a 5.2 Mb genome on 2.

## Validation

`scripts/fetch_test_refs.sh` downloads H37Rv, MRSA USA300_FPR3757 and *P. aeruginosa* PAO1.
`scripts/validate.sh` runs refmasker on them and compares the H37Rv mask with snippy's curated
TB mask. Results so far are in [docs/validation.md](docs/validation.md).

## Other commands

```
refmasker compare a.bed b.bed      bp in common, Jaccard, regions unique to each mask
refmasker check                    which external tools are installed
refmasker download-db DIR          download the geNomad database
```

## Citing

If you use refmasker in a paper, please cite the tools it runs. There is a ready-to-paste
methods paragraph with references in [docs/methods.md](docs/methods.md).

## Development

```bash
conda activate refmasker                 # environment.yml includes pytest and ruff
pip install --no-deps -e .               # editable for development, in your own environment only
pytest -m "not slow"       # fast tests, no external tools needed (what CI runs)
pytest -m slow             # runs the real tools; needs the conda env
ruff check src tests && ruff format --check src tests
```

## Publishing this repository

This repository was created locally. To publish it, create an empty repository on GitHub with no
README or licence, then:

```bash
git remote add origin git@github.com:<you>/refmasker.git
git push -u origin main
```

## Licence

MIT. See [LICENSE](LICENSE).
