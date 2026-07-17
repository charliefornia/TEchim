#!/bin/bash

# Usage: ./filter_gene_te_by_reads.sh input.tsv output.tsv 20

input="$1"
output="$2"
min_reads="$3"

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 input.tsv output.tsv min_total_reads"
    exit 1
fi

awk -F'\t' -v OFS='\t' -v MIN="$min_reads" '
{
    gene = $1
    te = $5
    read_count = $9 + 0  # force numeric

    # Normalise TE name by stripping _LTR
    te_core = te
    sub(/_LTR$/, "", te_core)

    group_key = gene "|" te_core
    total_reads[group_key] += read_count

    # Save original line and its normalised group key
    lines[NR] = $0
    keys[NR] = group_key
}
END {
    for (i = 1; i <= NR; i++) {
        if (total_reads[keys[i]] >= MIN) {
            print lines[i]
        }
    }
}
' "$input" | sort -t $'\t' -k2,2 -k4,4n > "$output"