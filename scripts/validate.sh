#!/usr/bin/env bash
# Run refmasker on the validation references and compare H37Rv with snippy's TB mask.
#
#   scripts/fetch_test_refs.sh test_data
#   scripts/validate.sh test_data runs [THREADS]
#
# Needs the conda environment and REFMASKER_GENOMAD_DB (or edit CATS to drop prophage).
set -euo pipefail

# Arguments: folder made by fetch_test_refs.sh, output folder, number of threads.
data=${1:-test_data}
out=${2:-runs}
threads=${3:-4}
mkdir -p "$out"

# 1. Build a mask with the default categories for each reference. The GFF3 annotation labels
#    each masked region with the genes it overlaps. --force reruns over previous output.
for ref in H37Rv USA300_FPR3757 PAO1; do
    refmasker run "$data/$ref.fna" --annotation "$data/$ref.gff" -o "$out/$ref" \
        --threads "$threads" --force
done

# 2. The TB mask in snippy is mostly PE/PPE genes, so also try masking those by name.
#    The regex matches PE (as a whole word), PPE, PE_PGRS and PE-PGRS in gene names,
#    products or locus tags.
refmasker run "$data/H37Rv.fna" --annotation "$data/H37Rv.gff" -o "$out/H37Rv_pe_ppe" \
    --mask-genes-matching '\bPE\b|PPE|PE_PGRS|PE-PGRS' --threads "$threads" --force

# 3. Compare both H37Rv masks with snippy's. Prints the overlap statistics and writes the
#    regions found by only one of the masks to vs_snippy.only_a.bed / vs_snippy.only_b.bed.
for run in H37Rv H37Rv_pe_ppe; do
    echo "== $run vs snippy's Mtb mask"
    refmasker compare "$out/$run/mask.bed" "$data/snippy_Mtb_NC_000962.3_mask.bed" \
        --out-prefix "$out/$run/vs_snippy"
done

# 4. Print the genome-wide total from each run's summary.tsv: the row with record "ALL" and
#    category "total"; column 6 is pct_masked.
echo
echo "Fraction masked:"
for ref in H37Rv H37Rv_pe_ppe USA300_FPR3757 PAO1; do
    awk -v r="$ref" '$1=="ALL" && $2=="total" {printf "  %-16s %s%%\n", r, $6}' "$out/$ref/summary.tsv"
done
