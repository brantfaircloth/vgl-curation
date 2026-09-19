#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 [-k] <input.fa> <output.fa>" >&2
    echo "  -k    keep intermediate files (default: delete them)" >&2
    exit 1
}

KEEP=0
while getopts ":k" opt; do
    case $opt in
        k) KEEP=1 ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -ne 2 ]; then
    usage
fi

IN=$1
OUT=$2

seqkit grep -nrp '^SUPER_' "$IN" | seqkit sort -N > super.fa
seqkit grep -nrp '^Scaffold_' "$IN" | seqkit sort -l -r -2 > scaf.sorted.fa
paste \
  <(seqkit seq -n -i scaf.sorted.fa) \
  <(seqkit seq -n -i scaf.sorted.fa | awk '{printf "Scaffold_%d\n", NR}') \
  > scaffold_rename.tsv
seqkit replace -p '^.+$' -r 'Scaffold_{nr}' scaf.sorted.fa > scaf.renamed.fa

cat super.fa scaf.renamed.fa \
    <(seqkit grep -vnrp '^(SUPER|Scaffold)_' "$IN") > "$OUT"

if [ "$KEEP" -eq 0 ]; then
    rm -f super.fa scaf.sorted.fa scaffold_rename.tsv scaf.renamed.fa
fi
