# Validation

## Status (v0.1.0)

| Reference | Run | Result |
| --- | --- | --- |
| *M. tuberculosis* H37Rv | repeats, IS, rRNA (+ low complexity, CRISPR) | Done; compared with snippy's curated mask below |
| *B. fragilis* FDAARGOS_1225 (RefSeq, chromosome + plasmid) | repeats, IS, rRNA, low complexity, CRISPR, plasmid | Done; 4.4% masked |
| snippy's `test/example` (3-contig draft, GenBank) | repeats, prophage (PhiSpy), rRNA | Done; PhiSpy path and draft-assembly warning work |
| MRSA USA300_FPR3757 | defaults | **To do:** run `scripts/validate.sh` |
| *P. aeruginosa* PAO1 | defaults | **To do:** run `scripts/validate.sh` |
| geNomad prophage calls | | **Not yet run end to end:** the build environment could not download the geNomad database. The output parser is unit-tested on geNomad's documented format |

The first runs were made in an environment without access to NCBI, so H37Rv came from the
[tbdb](https://github.com/jodyphelan/tbdb) copy of GCA_000195955.2 (identical sequence to
NC_000962.3, renamed) with its Ensembl GFF. *B. fragilis* came from the nf-core test datasets.
`scripts/fetch_test_refs.sh` downloads the RefSeq versions of H37Rv, USA300 and PAO1 directly.

## H37Rv vs snippy's TB mask

snippy ships `etc/Mtb_NC_000962.3_mask.bed`: 731 regions, 311 after merging, covering 209,178 bp
(4.7%). It is a curated list built mainly from PE/PPE genes, IS elements and phage regions.
refmasker builds its mask from sequence, so the two are not expected to match exactly. The useful
question is how much of the curated mask refmasker recovers, and what it adds.

`refmasker compare A.bed snippy_mask.bed`:

| refmasker run | bp masked | % genome | Covers this much of snippy's mask | Jaccard |
| --- | --- | --- | --- | --- |
| repeats, IS, rRNA (defaults without prophage) | 157,222 | 3.6% | 59.2% | 0.51 |
| + low complexity, CRISPR | 230,479 | 5.2% | 63.3% | 0.43 |
| repeats, IS, rRNA + `--mask-genes-matching '\bPE\b\|PPE\|PE_PGRS\|PE-PGRS'` | 410,560 | 9.3% | 74.7% | 0.34 |

Prophage calling was not part of these runs (see Status). 79% of the bases refmasker masks with
the defaults are also in snippy's mask, so what it adds is small.

What this shows:

- **Repeats, IS elements and rRNA recover about 60% of the curated mask.** Most of the rest is PE/PPE
  genes that are TB-specific and diverged below 90% identity, so self-alignment can't see them as
  repeats. The sensitivity table shows this: dropping to 80% identity raises the repeat fraction
  from 2.7% to 4.2% of the genome.
- **Masking PE/PPE genes by name brings coverage to 75%**, but masks 6.5% of the genome. snippy's
  curated mask covers only some PE/PPE genes, or parts of them. This is the trade-off the
  `genes` category is for: it's explicit and easy to justify, but blunt.
- **Low complexity adds about 1.7% of H37Rv** in 771 short intervals, because the genome is 65% GC.
  It contributes little to agreement with the curated mask, which is why it's opt-in.
- **ISEScan partial calls are excluded by default.** On H37Rv, one partial "IS21" call covered
  *dnaA* at the start of the genome; leaving partial calls out removes this and 8 others.

For the species refmasker is aimed at (MRSA, *P. aeruginosa*), no curated mask exists to compare
against. The H37Rv comparison checks that the detection approach recovers the bulk of what an
expert mask contains, and shows where it doesn't.

## B. fragilis FDAARGOS_1225 (NZ_CP069563.1 + plasmid NZ_CP069564.1)

| Category | Regions | bp | % genome |
| --- | --- | --- | --- |
| Repeats | 79 | 92,837 | 1.77 |
| IS / transposons | 14 | 12,361 | 0.24 |
| rRNA operons | 6 | 30,926 | 0.59 |
| Low complexity | 2,584 | 39,295 | 0.75 |
| CRISPR arrays | 2 | 1,851 | 0.04 |
| Plasmid (whole record) | 1 | 98,235 | 1.87 |
| **Total** | | **230,566** | **4.39** |

The six rRNA operons are each about 5.2 kb once the 16S–23S spacer is included. The plasmid is
recognised from its RefSeq header and masked whole because `--mask-plasmids` was used.

## Run times

On 4 threads of a cloud VM: H37Rv took 7 minutes, almost all of it ISEScan. On 2 threads,
*B. fragilis* (5.2 Mb) took 13 minutes. nucmer self-alignment takes under 10 seconds for a 5 Mb
genome.

## Reproduce

```bash
scripts/fetch_test_refs.sh test_data
scripts/validate.sh test_data runs 4
```
