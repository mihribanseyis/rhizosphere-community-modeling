#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=extract_RhizoMAGs
#SBATCH --mail-type=END,FAIL

###############################################################################
# 16_extract_hq_mags.sh
#
# Copy the DAS Tool-refined MAGs that passed the CheckM completeness and
# contamination thresholds into a separate output directory.
#
# Usage:
#   sbatch 16_extract_hq_mags.sh \
#       BIN_DIR SELECTED_BINS_TABLE OUTPUT_DIR
#
# The selected-bins table must contain:
#   bin_id  completeness  contamination
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 BIN_DIR SELECTED_BINS_TABLE OUTPUT_DIR" >&2
    exit 1
fi

bin_dir="$1"
selected_bins_table="$2"
output_dir="$3"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$bin_dir" ]]; then
    echo "ERROR: Bin directory not found: $bin_dir" >&2
    exit 1
fi

if [[ ! -s "$selected_bins_table" ]]; then
    echo "ERROR: Selected-bin table not found or empty: $selected_bins_table" >&2
    exit 1
fi

mkdir -p "$output_dir"


# ---------------------------------------------------------------------------
# Extract selected MAGs
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting selected RhizoMAG extraction"
echo "Bin directory:       $bin_dir"
echo "Selected-bin table:  $selected_bins_table"
echo "Output directory:    $output_dir"
echo

copied_count=0
missing_count=0

while IFS=$'\t' read -r bin_id completeness contamination _; do
    [[ -z "${bin_id:-}" ]] && continue
    [[ "$bin_id" == "Bin Id" ]] && continue

    source_file="${bin_dir}/${bin_id}.fa"
    destination_file="${output_dir}/${bin_id}.fa"

    if [[ ! -s "$source_file" ]]; then
        echo "ERROR: FASTA file not found or empty for ${bin_id}: ${source_file}" >&2
        missing_count=$((missing_count + 1))
        continue
    fi

    cp "$source_file" "$destination_file"

    echo "Copied ${bin_id}: completeness=${completeness}, contamination=${contamination}"

    copied_count=$((copied_count + 1))

done < "$selected_bins_table"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_seconds=$((end_time - start_time))

echo
echo "[$(date)] RhizoMAG extraction completed"
echo "MAGs copied:  $copied_count"
echo "MAGs missing: $missing_count"
echo "Output:       $output_dir"
echo "Runtime:      ${runtime_seconds} seconds"

if [[ "$copied_count" -eq 0 ]]; then
    echo "ERROR: No MAG files were copied." >&2
    exit 1
fi

if [[ "$missing_count" -gt 0 ]]; then
    echo "ERROR: Some selected MAG files could not be found." >&2
    exit 1
fi