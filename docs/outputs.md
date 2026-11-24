# Outputs

Everything goes into the `-o` folder. A new folder gets no access for other users (group access
follows your umask), and refmasker refuses an existing folder that belongs to someone else, that
any user can write to, or (with `--force`) that contains symbolic links.

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
(and the geNomad database version when prophages are called), full paths of the tools used
(`tool_paths`), total bp and % masked, warnings.

With `--redact`, paths are reduced to file names, header descriptions are left out and
`"redacted": true` is recorded; `tool_paths` is `null`. Record IDs are always kept: if your FASTA IDs are themselves
sample identifiers, rename them before running.

## refmasker.log

The full log, including every external command, with full paths even under `--redact`: don't
share it if paths are sensitive. With `--keep-work`, each tool's stderr is in
`work/<tool>/<tool>.log` next to its raw output.

## logs/

Only after a failed run (without `--keep-work`): the stderr log of each tool that ran, moved out of
`work/` before it is deleted.

## The comparison report (`refmasker compare`)

`refmasker compare A.bed B.bed` writes one self-contained HTML page, named after the two inputs
(`compare_<A>_vs_<B>.html`) or after `--out-prefix` (`PREFIX.compare.html`). Both figures are
inlined, so the page needs no other file; `--formats html,png` also saves them as
`<base>.bar.png` and `<base>.tracks.png`, and `--formats none` writes nothing at all.

The page holds:

| Section | What it shows |
| --- | --- |
| Tiles | bp masked by each mask, bp masked by both, Jaccard index, and what fraction of each mask the other covers. Percentages of the genome only with `--ref` |
| How much the masks share | One bar per mask, split into bases the other mask also covers and bases unique to it |
| Where they agree and differ | One row per mask along each sequence: green where both masks cover a stretch, the mask's own colour where only it does |
| What kind of masking the differences are | A bar per mask category: grey for the part the other mask also covers, the category's own colour for what it misses. Only with an annotated BED |
| By category | Per category: regions, bp, bp the other mask also covers, bp only this mask has, and the percentage the other misses. Only with an annotated BED |
| By sequence | Per-sequence bp for each mask, shared bp, bp unique to each, and the Jaccard index for that sequence |
| Regions the other mask misses entirely | Every merged region of one mask that shares no base with the other, sortable and searchable, with the category and name behind it when an annotated BED supplied them |
| Run details | The two paths, the reference, any `--rename`, and when the page was written |

Two differences are reported separately, because they answer different questions: *bases only in
one mask* is a per-base difference (a region the other mask half covers counts half), while
*regions the other mask misses entirely* counts whole regions with no overlap at all.

`--out-prefix` additionally writes `PREFIX.only_a.bed` and `PREFIX.only_b.bed`, the merged
regions each mask has that the other does not touch.

### Without matplotlib

The figures need matplotlib; the rest of the page does not. If it is missing or broken, the
report is still written, with every table and number, and a note at the top saying the figures
were left out and why. `refmasker run` needs it for `report.pdf` and the figures.

### Where the categories come from

`mask.bed` has three columns and no categories. They live in `mask.annotated.bed`, so the
category sections appear only when one is available. refmasker looks, in order, at
`--annotated-a` / `--annotated-b` (the literal `none` switches it off), the mask itself when it
is already an annotated BED, and a `mask.annotated.bed` beside a plain one. It reports which
file it used, in the log and in the report's run details, because the last of those reads a file
you did not name. Comparing two refmasker runs gives categories for both masks; comparing
against a published three-column mask gives them for yours only.
