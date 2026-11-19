# FAQ

**How much of the genome should be masked?**
There is no fixed target. With the default categories, expect a few percent for most species:
3.6% of H37Rv with repeats, IS and rRNA, and 4.4% of a *B. fragilis* genome including its plasmid (see [validation.md](validation.md)). Genomes rich in IS
elements or prophages will be higher. A mask over about 10% deserves a look at the report:
usually one category (often `lowcomplexity` or a loose repeat threshold) is driving it.

**Does it matter that the reference has plasmids?**
Plasmid records are kept and searched like the chromosome, but not masked whole unless you ask for
it with `--mask-plasmids`. Repeats between chromosome and plasmid (shared IS elements, for example)
are found either way.

**My reference is a draft assembly.**
It works, but the report warns you. Repeats and IS elements tend to be collapsed or broken at
contig ends in short-read assemblies, so they are under-detected. A closed RefSeq genome is a
better reference.

**Do I need the annotation?**
No, but it helps. With `--annotation` the region table lists the genes each region removes, `is`
also uses annotated transposases, and you can use `--mask-genes-matching`. Use the GFF3 or GenBank
file that came with the same RefSeq assembly, so the sequence IDs match the FASTA. The report warns
if they don't.

**How do I use the mask with nullarbor?**
nullarbor runs snippy and snippy-core for you. If your nullarbor version has no option to pass a
mask through, re-run snippy-core yourself on the per-isolate snippy folders nullarbor made:

```bash
snippy-core --ref ref.fna --mask ref_mask/mask.bed --prefix core_masked nullarbor_out/*/
```

Check that each folder you pass is a snippy output (it contains `snps.aligned.fa`).

**snippy-core says "BED ... has CHR that is not in the reference".**
The mask was built from a different FASTA, or one whose IDs differ (e.g. `NC_007793` vs
`NC_007793.1`). Use exactly the same FASTA for refmasker and snippy.

**Can I combine refmasker's mask with snippy's own `--mask auto`?**
snippy-core's `--mask auto` masks repeats found by a minimap2 self-alignment of the reference.
That overlaps with refmasker's `repeats` category but covers nothing else, so use refmasker's
`mask.bed` instead.

**Can I run it on many references at once?**
Yes, with a loop or job array; each run is independent:

```bash
for ref in refs/*.fna; do
  st=$(basename "$ref" .fna)
  refmasker run "$ref" --annotation "refs/$st.gff" -o "masks/$st" --threads 4
done
```

A Nextflow wrapper is planned.

**Why is recombination not masked?**
Recombination depends on your isolates, not the reference. Run Gubbins or ClonalFrameML on the
masked core alignment.
