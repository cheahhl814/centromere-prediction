#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 4 ]; then
  echo "Usage: $0 <density_bed> [scaffolds_file] [output_file] [top_n]" >&2
  exit 1
fi

density_bed="$1"
scaffolds_file="${2:-scaffolds.txt}"
output_file="${3:-candidate_regions.tsv}"
top_n="${4:-10}"

if [ ! -f "$density_bed" ]; then
  echo "Error: density BED file not found: $density_bed" >&2
  exit 1
fi

if [ ! -f "$scaffolds_file" ]; then
  echo "Error: scaffolds file not found: $scaffolds_file" >&2
  exit 1
fi

: > "$output_file"

while read -r scaffold; do
  [ -z "$scaffold" ] && continue

  awk -v sc="$scaffold" 'BEGIN{FS="\t"; OFS="\t"} $1 == sc' "$density_bed" | \
    sort -k7,7rn | sed -n "1,${top_n}p" | sort -k2,2n | \
    awk 'BEGIN{FS="\t"; OFS="\t"} NR==1{s=$2; chr=$1} {e=$3} END{if(chr) print chr, s, e}' \
    >> "$output_file"
done < "$scaffolds_file"
