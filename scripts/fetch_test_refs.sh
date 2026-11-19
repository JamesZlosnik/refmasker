#!/usr/bin/env bash
# Download the validation references (FASTA + GFF3) from NCBI, and snippy's TB mask.
#
#   scripts/fetch_test_refs.sh [OUTDIR]     (default: test_data/)
#
# References:
#   H37Rv           NC_000962.3                         M. tuberculosis (benchmark vs snippy's mask)
#   USA300_FPR3757  NC_007793.1 + plasmids NC_007790-2  MRSA ST8
#   PAO1            NC_002516.2                         P. aeruginosa
set -euo pipefail

out=${1:-test_data}
mkdir -p "$out"
eutils="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

fetch() {  # fetch NAME ACC [ACC...]  -> NAME.fna, NAME.gff
    local name=$1; shift
    local ids
    ids=$(IFS=,; echo "$*")
    echo "Downloading $name ($ids)"
    curl -fsSL "$eutils?db=nuccore&id=$ids&rettype=fasta&retmode=text" > "$out/$name.fna"
    : > "$out/$name.gff"
    for acc in "$@"; do
        curl -fsSL "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$acc" \
            | { if [ -s "$out/$name.gff" ]; then grep -v '^##gff-version'; else cat; fi; } >> "$out/$name.gff"
        sleep 0.5  # stay under NCBI's request limit
    done
}

fetch H37Rv NC_000962.3
fetch USA300_FPR3757 NC_007793.1 NC_007790.1 NC_007791.1 NC_007792.1
fetch PAO1 NC_002516.2

curl -fsSL https://raw.githubusercontent.com/tseemann/snippy/master/etc/Mtb_NC_000962.3_mask.bed \
    > "$out/snippy_Mtb_NC_000962.3_mask.bed"

grep -c '>' "$out"/*.fna
echo "Done: $out"
