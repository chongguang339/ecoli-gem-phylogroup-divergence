#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  build_gapseq_compatible_ecoli_fastdb.sh -g <gapseq_bacteria_seq_dir> -o <output_faa> -r <regex>

Example:
  bash build_gapseq_compatible_ecoli_fastdb.sh \
    -g /path/to/gapseq/dat/seq/Bacteria \
    -o gapseq_assets/enterobacterales_gapseq_fastdb.faa \
    -r 'Escherichia|Shigella|Salmonella|Klebsiella|Enterobacter|Citrobacter'

This script builds a compact gapseq-compatible FASTA database by selecting
sequence records whose FASTA headers match the supplied taxonomic regex.
EOF
}

GENOME_DB=""
OUTPUT=""
REGEX=""

while getopts ":g:o:r:h" opt; do
  case "$opt" in
    g) GENOME_DB="$OPTARG" ;;
    o) OUTPUT="$OPTARG" ;;
    r) REGEX="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

if [[ -z "$GENOME_DB" || -z "$OUTPUT" || -z "$REGEX" ]]; then
  usage
  exit 1
fi

if [[ ! -d "$GENOME_DB" ]]; then
  echo "Input gapseq sequence directory does not exist: $GENOME_DB" >&2
  exit 2
fi

mkdir -p "$(dirname "$OUTPUT")"

awk -v regex="$REGEX" '
  /^>/ {
    keep = ($0 ~ regex)
  }
  keep {
    print
  }
' "$GENOME_DB"/*.faa "$GENOME_DB"/*.fasta 2>/dev/null > "$OUTPUT"

if [[ ! -s "$OUTPUT" ]]; then
  echo "No records were written. Check input file extensions and regex." >&2
  exit 3
fi

echo "Wrote custom gapseq FASTA database: $OUTPUT"
