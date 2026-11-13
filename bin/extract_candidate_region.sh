#!/usr/bin/env bash
# Pick one candidate centromere region per scaffold from a bedtools-coverage density BED.
#
#   window_bp > 0 (default 1000000): the contiguous run of windows with the highest summed
#       density, at most window_bp long. Bounded, so stringDecomposer memory stays predictable.
#   window_bp = 0 (legacy): the span of the top_n densest windows merged into one region.
#       Scattered hits give a near-whole-scaffold span (the v1.0.0 behaviour that OOM-killed
#       stringDecomposer), so only use it for comparison.
#
# Output: tab-separated  scaffold  start(0-based)  end

set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 5 ]; then
  echo "Usage: $0 <density_bed> [scaffolds_file] [output_file] [top_n] [window_bp]" >&2
  exit 1
fi

density_bed="$1"
scaffolds_file="${2:-scaffolds.txt}"
output_file="${3:-candidate_regions.tsv}"
top_n="${4:-10}"
window_bp="${5:-1000000}"

if [ ! -f "$density_bed" ]; then
  echo "Error: density BED file not found: $density_bed" >&2
  exit 1
fi

if [ ! -f "$scaffolds_file" ]; then
  echo "Error: scaffolds file not found: $scaffolds_file" >&2
  exit 1
fi

: > "$output_file"

if [ "$window_bp" -gt 0 ]; then
  # Sliding window over each scaffold's bins (density = column 7). Ties keep the first window.
  while read -r scaffold; do
    [ -z "$scaffold" ] && continue
    awk -v sc="$scaffold" -v wbp="$window_bp" 'BEGIN{FS="\t"; OFS="\t"}
      $1 == sc { n++; s[n]=$2; e[n]=$3; d[n]=$7 }
      END {
        if (n == 0) exit
        bin = e[1] - s[1]; k = int(wbp / bin); if (k < 1) k = 1; if (k > n) k = n
        for (i = 1; i <= k; i++) cur += d[i]
        best = cur; bi = 1
        for (i = 2; i + k - 1 <= n; i++) {
          cur += d[i + k - 1] - d[i - 1]
          if (cur > best) { best = cur; bi = i }
        }
        print sc, s[bi], e[bi + k - 1]
      }' "$density_bed" >> "$output_file"
  done < "$scaffolds_file"
else
  while read -r scaffold; do
    [ -z "$scaffold" ] && continue

    awk -v sc="$scaffold" 'BEGIN{FS="\t"; OFS="\t"} $1 == sc' "$density_bed" | \
      sort -k7,7rn | sed -n "1,${top_n}p" | sort -k2,2n | \
      awk 'BEGIN{FS="\t"; OFS="\t"} NR==1{s=$2; chr=$1} {e=$3} END{if(chr) print chr, s, e}' \
      >> "$output_file"
  done < "$scaffolds_file"
fi
