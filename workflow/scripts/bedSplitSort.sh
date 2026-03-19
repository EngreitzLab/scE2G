#!/bin/bash
## Memory-efficient replacement for `bedtools sort -g sizes` that sorts BED records
## in the chromosome order defined by a genome sizes file.
##
## Strategy: split input by chromosome, sort each file by position, then concatenate
## in the order specified by the sizes file.
##
## Uses `sort-bed` from bedops if available (~2x faster), otherwise falls back to
## GNU sort. Supports reading from stdin or from a file (optionally gzipped).

set -euo pipefail

USAGE="Usage: $(basename "$0") -i input.bed(.gz)|/dev/stdin -g chromosome_sizes [-t split_dir] [-p threads] [-s buffer_size] [-T sort_tmp_dir]"

INPUT=
SIZES=
SPLITDIR=
THREADS=1
BUFFER_SIZE=
SORT_TMPDIR=

while getopts ":i:g:t:p:s:T:h" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    g) SIZES="$OPTARG" ;;
    t) SPLITDIR="$OPTARG" ;;
    p) THREADS="$OPTARG" ;;
    s) BUFFER_SIZE="$OPTARG" ;;
    T) SORT_TMPDIR="$OPTARG" ;;
    h) echo "$USAGE"; exit 0 ;;
    *) echo "Unknown option: -$OPTARG" >&2; echo "$USAGE" >&2; exit 1 ;;
  esac
done

if [ -z "$INPUT" ] || [ -z "$SIZES" ]; then
    echo "$USAGE" >&2
    exit 1
fi

if [ -z "$SPLITDIR" ]; then
    SPLITDIR=$(mktemp -d "${INPUT}.split.XXXXXX")
fi

mkdir -p "$SPLITDIR"

## Determine sort command: prefer sort-bed if available, fall back to GNU sort
if command -v sort-bed &>/dev/null; then
    sort_chr() {
        sort-bed "$1"
    }
else
    sort_chr() {
        local sort_args=(-k2,2n -k3,3n)
        if [ "$THREADS" -gt 1 ]; then
            sort_args+=(--parallel "$THREADS")
        fi
        if [ -n "$BUFFER_SIZE" ]; then
            sort_args+=(-S "$BUFFER_SIZE")
        fi
        if [ -n "$SORT_TMPDIR" ]; then
            sort_args+=(-T "$SORT_TMPDIR")
        fi
        LC_ALL=C sort "${sort_args[@]}" "$1"
    }
fi

## Split input by chromosome
if [ "$INPUT" = "-" ] || [ "$INPUT" = "/dev/stdin" ]; then
    awk -v tmpdir="$SPLITDIR" '{ print $0 > tmpdir "/" $1 ".bed" }'
else
    zcat --force "$INPUT" | awk -v tmpdir="$SPLITDIR" '{ print $0 > tmpdir "/" $1 ".bed" }'
fi

## Warn about chromosomes in the input that are not in the sizes file
for file in "$SPLITDIR"/*.bed; do
    [ -f "$file" ] || continue
    chr=$(basename "$file" .bed)
    if ! awk -v c="$chr" '$1 == c { found=1; exit } END { exit !found }' "$SIZES"; then
        >&2 echo "Warning: Found entries in BED file for ${chr}, but this chromosome is not in the genome sizes file ($SIZES)"
    fi
done

## Output sorted records in the chromosome order from the sizes file
while read -r chr size; do
    if [ -f "$SPLITDIR/${chr}.bed" ]; then
        sort_chr "$SPLITDIR/${chr}.bed"
    fi
done < "$SIZES"

## Clean up
rm -r "$SPLITDIR"
