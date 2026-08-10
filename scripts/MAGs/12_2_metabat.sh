#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --job-name=metabat_global
#SBATCH --mail-type=END,FAIL

###############################################################################
# 12_2_metabat_global.sh
#
# Perform global MetaBAT2 binning using the shared rhizosphere assembly and
# the combined depth file generated from all rhizosphere samples.
#
# Usage:
#   sbatch 12_2_metabat_global.sh \
#       ASSEMBLY DEPTH_FILE OUTPUT_DIR [SUMMARY_LOG]
#
# Example:
#   sbatch 12_2_metabat_global.sh \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa \
#       data/mappingRh/all_samples_combined_depth.txt \
#       data/binningRh/all_samples/metabat \
#       logs/bins/metabat_summary.tsv
#
# Requirements:
#   MetaBAT2 2.17
#
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

CONDA_PROFILE="${CONDA_PROFILE:-$HOME/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="binning-env"

if [[ ! -f "$CONDA_PROFILE" ]]; then
    echo "ERROR: Cannot find conda.sh: $CONDA_PROFILE" >&2
    exit 1
fi

source "$CONDA_PROFILE"
conda activate "$CONDA_ENV"

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: sbatch $0 ASSEMBLY DEPTH_FILE OUTPUT_DIR [SUMMARY_LOG]" >&2
    exit 1
fi

assembly="$1"
depth_file="$2"
output_dir="$3"
summary_log="${4:-}"

threads="${SLURM_CPUS_PER_TASK:-8}"
output_prefix="${output_dir}/allRh_metabat_bin"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -s "$assembly" ]]; then
    echo "ERROR: Assembly file not found or empty: $assembly" >&2
    exit 1
fi

if [[ ! -s "$depth_file" ]]; then
    echo "ERROR: Combined depth file not found or empty: $depth_file" >&2
    exit 1
fi

if ! command -v metabat2 >/dev/null 2>&1; then
    echo "ERROR: metabat2 is not available in the active environment." >&2
    exit 1
fi

mkdir -p "$output_dir"

if [[ -n "$summary_log" ]]; then
    mkdir -p "$(dirname "$summary_log")"
fi


# ---------------------------------------------------------------------------
# Run MetaBAT2
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting global MetaBAT2 binning"
echo "Assembly:       $assembly"
echo "Depth file:     $depth_file"
echo "Output prefix:  $output_prefix"
echo "Threads:        $threads"
echo "Unbinned files: enabled"
echo

metabat2 \
    -i "$assembly" \
    -a "$depth_file" \
    -o "$output_prefix" \
    -t "$threads" \
    --unbinned


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Global MetaBAT2 binning completed."
echo "Output directory: $output_dir"
echo "Runtime:          ${runtime_minutes} minutes"


# ---------------------------------------------------------------------------
# Optional summary table
# ---------------------------------------------------------------------------

if [[ -n "$summary_log" ]]; then
    if [[ ! -f "$summary_log" ]]; then
        printf "Run\tJob_ID\tDuration_min\tEnd_Time\tOutput_Dir\n" > "$summary_log"
    fi

    printf "all_samples\t%s\t%s\t%s\t%s\n" \
        "${SLURM_JOB_ID:-NA}" \
        "$runtime_minutes" \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$output_dir" \
        >> "$summary_log"

    echo "Summary log:      $summary_log"
fi