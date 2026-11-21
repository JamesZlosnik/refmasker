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

# Output folder: the first argument, or test_data/ if none is given.
out=${1:-test_data}
mkdir -p "$out"
# NCBI E-utilities endpoint that returns sequences as FASTA.
eutils="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# Pinned so the benchmark mask cannot change under us.
snippy_commit=feebb229c5e448af65167fbad9e454e2366c47b8  # tseemann/snippy master, 2026-11-19

# fetch NAME ACC [ACC...]
#   Download one genome (one or more accessions, e.g. chromosome + plasmids) into
#   $out/NAME.fna (all records in one multi-FASTA) and $out/NAME.gff (one combined GFF3).
fetch() {
    local name=$1; shift       # first argument is the output name; the rest are accessions
    local ids
    ids=$(IFS=,; echo "$*")    # join the accessions with commas, as efetch expects
    echo "Downloading $name ($ids)"
    # One efetch request returns every record as a single multi-FASTA.
    curl -fsSL "$eutils?db=nuccore&id=$ids&rettype=fasta&retmode=text" > "$out/$name.fna"
    # GFF3 comes one accession at a time from NCBI's sequence viewer. Start from an empty
    # file, then append each accession's GFF3, keeping only the first '##gff-version' line
    # so the combined file has a single header.
    : > "$out/$name.gff"
    for acc in "$@"; do
        curl -fsSL "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$acc" \
            | { if [ -s "$out/$name.gff" ]; then grep -v '^##gff-version'; else cat; fi; } >> "$out/$name.gff"
        sleep 0.5  # stay under NCBI's request limit
    done
}

# The three validation genomes (see the list at the top of this file).
fetch H37Rv NC_000962.3
fetch USA300_FPR3757 NC_007793.1 NC_007790.1 NC_007791.1 NC_007792.1
fetch PAO1 NC_002516.2

# snippy's published M. tuberculosis mask: the benchmark that validate.sh compares against.
curl -fsSL "https://raw.githubusercontent.com/tseemann/snippy/${snippy_commit}/etc/Mtb_NC_000962.3_mask.bed" \
    > "$out/snippy_Mtb_NC_000962.3_mask.bed"

# Sanity check: print the number of FASTA records in each downloaded genome.
grep -c '>' "$out"/*.fna
echo "Done: $out"
