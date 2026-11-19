# Outputs

Everything goes into the `-o` folder.

## mask.bed

The file for snippy. Regions from all selected categories are padded (`--pad`), clipped to sequence
ends, and merged, so overlapping or touching regions become one line. It has three tab-separated
columns (sequence ID, 0-based start, exclusive end) and **no header**: snippy-core reads every line
as a region and rejects lines whose first column is not a reference sequence. Lines follow the
order of records in the FASTA.

## mask.annotated.bed

Every region before merging across categories, one line per region, with a `#` header line. Within the `repeats` category, the pairwise alignments of a multi-copy element are combined into one region per copy.

| Column | Meaning |
| --- | --- |
| chrom, start, end | BED coordinates (0-based, end exclusive) |
| category | `repeats`, `is`, `prophage`, ... |
| source | tool or input that found it: `nucmer`, `ISEScan`, `annotation`, `geNomad`, `barrnap`, ... |
| length | bp |
| name | detail from the tool, e.g. `matches 3 other copies, up to 99.6% identity`, `IS6/IS6_1 (complete)`, `operon 16S,23S,5S` |
| genes | genes the region overlaps, from `--annotation` (`locus_tag(gene)`), at most 50 listed |

Load it in IGV or the UCSC browser next to the annotation, or filter it with `awk` to build a mask
from a subset of categories:

```bash
awk -F'\t' '!/^#/ && ($4=="repeats" || $4=="is")' mask.annotated.bed | cut -f1-3 \
  | sort -k1,1 -k2,2n | bedtools merge > repeats_is.bed
```

## summary.tsv

| record | category | regions | bp_masked | record_length | pct_masked |
| --- | --- | --- | --- | --- | --- |

One block for the whole genome (`record` = `ALL`), then one per FASTA record. Within a block there
is one row per category plus `total`. A category's bp counts each base once, even where its own
regions overlap. `total` counts each base once across all categories, so it is less than the sum of
the category rows when categories overlap. `regions` is the number of merged intervals.

To compare masks across many ST references:

```bash
for d in */; do awk -v ref="${d%/}" '$1=="ALL" {print ref"\t"$0}' "$d/summary.tsv"; done
```

## report.html

A single file with no external dependencies: open it offline or attach it to a run. It has:

- headline numbers (genome size, bp and % masked, number of regions)
- warnings (draft assembly, Ns in the reference, annotation names not matching the FASTA)
- a linear and a circular genome plot, with one track or ring per category plus an "All masked"
  track
- bp and % masked per category and per record
- the reason for each category and which tools found its regions
- the repeat threshold sensitivity table
- a sortable, searchable table of every region and the genes it overlaps
- run details: command, parameters, tool versions, reference checksum

## report.pdf

The same content as the HTML, for printing or archiving. The region table is truncated to fit the
page; the full list is in `mask.annotated.bed`.

## figures/

`linear.svg`, `circular.svg`, and PNG versions at 200 dpi when `png` is in `--formats`. SVG text is
kept as text, so you can edit the figures in Inkscape or Illustrator. References with more than
six records (draft assemblies) are drawn with contigs placed end to end.

## run_info.json

Machine-readable run record: refmasker version, command, reference path with MD5 and SHA-256,
records (ID, length and FASTA header description), plasmid records, parameters, tool versions
(and the geNomad database version when prophages are called), total bp and % masked, warnings.

With `--redact`, paths are reduced to file names, header descriptions are left out and
`"redacted": true` is recorded. Record IDs are always kept: if your FASTA IDs are themselves
sample identifiers, rename them before running.

## refmasker.log

The full log, including every external command, with full paths even under `--redact`: don't
share it if paths are sensitive. With `--keep-work`, each tool's stderr is in
`work/<tool>/<tool>.log` next to its raw output.
