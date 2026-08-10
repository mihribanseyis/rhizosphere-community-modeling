#!/usr/bin/env bash

###############################################################################
# 09_bwa_index.sh
#
# Build a BWA index for the filtered assembly contigs.
#
# Usage:
#   bash 09_bwa_index.sh CONTIGS
#
# Example:
#   bash 09_bwa_index.sh \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa
###############################################################################

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash $0 CONTIGS" >&2
    exit 1
fi

contigs="$1"

if [[ ! -f "$contigs" ]]; then
    echo "ERROR: Contig file not found: $contigs" >&2
    exit 1
fi

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: BWA is not available." >&2
    exit 1
fi

bwa index "$contigs"

echo "BWA index created for: $contigs"