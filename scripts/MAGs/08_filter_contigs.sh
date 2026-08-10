#!/usr/bin/env bash

###############################################################################
# 08_filter_contigs.sh
#
# Retain assembled contigs with a minimum length of 2,000 bp.
#
# Usage:
#   bash 08_filter_contigs.sh INPUT_CONTIGS OUTPUT_CONTIGS
#
# Example:
#   bash 08_filter_contigs.sh \
#       data/megahit_assembly/mergedRh/final.contigs.fa \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa
###############################################################################

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 INPUT_CONTIGS OUTPUT_CONTIGS" >&2
    exit 1
fi

input_contigs="$1"
output_contigs="$2"

if [[ ! -f "$input_contigs" ]]; then
    echo "ERROR: Input contig file not found: $input_contigs" >&2
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    echo "ERROR: SeqKit is not available in the active environment." >&2
    exit 1
fi

mkdir -p "$(dirname "$output_contigs")"

seqkit seq \
    -m 2000 \
    "$input_contigs" \
    -o "$output_contigs"

echo "Filtered contigs written to: $output_contigs"