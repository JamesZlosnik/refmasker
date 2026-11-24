# Parameters

`refmasker run REFERENCE -o OUTDIR [options]`

| Option | Default | Description |
| --- | --- | --- |
| `REFERENCE` | required | Reference FASTA, single or multi-record, plain or gzipped |
| `-o, --outdir` | required | Output folder. Must be empty or absent unless `--force` |
| `--annotation FILE` | none | GFF3 or GenBank for the same reference (format detected automatically). Adds overlapping genes to the region table, transposase genes to `is`, and is needed for `genes` and `--phage-tool phispy` |
| `--categories LIST` | `repeats,is,prophage,rrna` | Comma-separated categories. Also accepts `default` and `all` (everything needing no extra input), e.g. `default,crispr` |
| `--read-length N` | 150 | Read length the mask is built for. Sets the default `--repeat-min-len` |
| `--repeat-min-len N` | max(200, read length) | Minimum repeat length (bp); both copies must reach it |
| `--repeat-min-id PCT` | 90 | Minimum % identity between repeat copies |
| `--phage-tool` | `genomad` | `genomad` or `phispy` |
| `--genomad-db DIR` | `$REFMASKER_GENOMAD_DB` | geNomad database folder (the `genomad_db` directory), downloaded once with `refmasker download-db` ([README](../README.md#the-genomad-database)) |
| `--isescan-partial` | off | Also mask IS elements ISEScan calls partial |
| `--mask-genes-matching REGEX` | none | Mask annotated genes whose product, gene name or locus tag matches (case-insensitive). Turns on `genes` |
| `--rrna-operon-gap N` | 1000 | Join rRNA genes this close into one operon region; 0 keeps them separate |
| `--homopolymer-min-len N` | 8 | Shortest homopolymer masked by `lowcomplexity` |
| `--extra-bed BED` | none | Your own regions (category `user`); repeat for several files |
| `--mask-plasmids` | off | Mask whole plasmid records (category `plasmid`) |
| `--plasmid-ids LIST` | none | Record IDs to treat as plasmids, besides headers mentioning "plasmid" |
| `--pad N` | 0 | Extend every region by N bp each side (clipped at sequence ends) |
| `--formats LIST` | `html,pdf,png` | Report formats; `none` writes only the BED and TSV files |
| `--threads N` | 1 or scheduler allocation | Threads for external tools. Without `--threads`, uses `SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, `NSLOTS`, `NCPUS` or `PBS_NUM_PPN` (first one set), otherwise 1. Warns if it is more than the job's allocation |
| `--max-genome-mb N` | 50 | Refuse a reference with more than N Mb of sequence; 0 = no limit |
| `--keep-work` | off | Keep raw tool outputs in `OUTDIR/work`. Otherwise `work/` is deleted, even after a failed run (failed runs keep the tool logs in `OUTDIR/logs`) |
| `--force` | off | Overwrite a previous refmasker output folder. Refuses any non-empty folder without `run_info.json` or `refmasker.log`, so it cannot clear an unrelated folder (such as a Nextflow `work/`), and any folder containing symbolic links |
| `--redact` | off | For sharing reports: HTML, PDF and `run_info.json` show input file names instead of full paths and leave out FASTA header descriptions. Record IDs stay (the mask needs them); `refmasker.log` keeps full paths |
| `-q, --quiet` | off | Print only warnings and errors (the log file still gets everything) |

## Choosing the repeat thresholds

The defaults suit 2x150 bp Illumina data. For 2x250 or 2x300 data, set `--read-length` and the
minimum repeat length follows. Look at the sensitivity table in the report: if % masked changes a
lot between 90% and 95% identity, your genome has many moderately diverged paralogues, and you
should decide deliberately which side of that line you want.

## Thread use

All external tools get the same thread count. refmasker also limits OpenMP/BLAS thread pools to
that number, so numerical libraries don't start one thread per core, and sets
`TF_NUM_INTRAOP_THREADS` (TensorFlow, used by geNomad) and `-XX:ActiveProcessorCount` in
`JAVA_TOOL_OPTIONS` (Java, used by MinCED) to match. Whether geNomad and MinCED stay within the
limit has not yet been measured on a cluster; check with `seff` or your scheduler's job report. ISEScan is by far the
slowest step: 12 minutes for a 5.2 Mb *B. fragilis* genome on 2 threads, and 7 minutes for H37Rv
on 4. nucmer, barrnap, dustmasker and MinCED each take seconds. geNomad's run time is not yet
measured (see [validation.md](validation.md)).

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | Success |
| 1 | No command given; or `refmasker check` found a default tool missing |
| 2 | Input, option or tool error, or a database that failed `download-db --verify` (message printed) |
