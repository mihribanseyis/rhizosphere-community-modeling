#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --job-name=gtdbtk_RhizoMAGs
#SBATCH --mail-type=END,FAIL

###############################################################################
# 17_gtdbtk_classification.sh
#
# Perform taxonomic classification of the selected RhizoMAGs using GTDB-Tk.
#
# Usage:
#   sbatch 17_gtdbtk_classification.sh GENOME_DIR OUTPUT_DIR
#
# Requirements:
#   GTDB-Tk 2.1.1
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 GENOME_DIR OUTPUT_DIR" >&2
    exit 1
fi

genome_dir="$1"
output_dir="$2"

threads="${SLURM_CPUS_PER_TASK:-8}"


# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

conda_profile="${CONDA_PROFILE:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
conda_environment="${GTDBTK_CONDA_ENV:-gtdbtk-2.1.1}"

if [[ ! -f "$conda_profile" ]]; then
    echo "ERROR: Conda initialization script not found: $conda_profile" >&2
    exit 1
fi

source "$conda_profile"
conda activate "$conda_environment"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$genome_dir" ]]; then
    echo "ERROR: Genome directory not found: $genome_dir" >&2
    exit 1
fi

if ! command -v gtdbtk >/dev/null 2>&1; then
    echo "ERROR: gtdbtk is not available in the active environment." >&2
    exit 1
fi

mag_count=$(
    find "$genome_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.fa' \
        | wc -l
)

if [[ "$mag_count" -eq 0 ]]; then
    echo "ERROR: No MAG files found in: $genome_dir" >&2
    exit 1
fi

mkdir -p "$output_dir"


# ---------------------------------------------------------------------------
# Run GTDB-Tk
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting GTDB-Tk classification"
echo "Genome directory: $genome_dir"
echo "Number of MAGs:   $mag_count"
echo "Output directory: $output_dir"
echo "Threads:          $threads"
echo

gtdbtk classify_wf \
    --genome_dir "$genome_dir" \
    --out_dir "$output_dir" \
    --extension fa \
    --skip_ani_screen \
    --cpus "$threads"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] GTDB-Tk classification completed"
echo "MAGs classified: $mag_count"
echo "Output:          $output_dir"
echo "Runtime:         ${runtime_minutes} minutes"