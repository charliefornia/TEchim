#!/bin/bash

# Usage: ./TEchim_postproc_SumUpnonTEsplicerate_v0.1.sh input.tsv output.tsv

input="$1"
output="$2"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.tsv output.tsv"
    exit 1
fi

gawk -F'\t' -v OFS='\t' '
{
    if (NF < 16) {
        print $0
        next
    }

    sum = 0
    text = $16
    while (match(text, /\(([0-9]+)\)/, m)) {
        sum += m[1]
        text = substr(text, RSTART + RLENGTH)
    }
    print $0, sum
}
' "$input" > "$output"