#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=drep_RhizoMAGs
#SBATCH --mail-type=END,FAIL

###############################################################################
# 18_dereplication.sh
#
# Dereplicate the selected RhizoMAGs using dRep.
#
# Usage:
#   sbatch 18_dereplication.sh INPUT_MAG_DIR OUTPUT_DIR
#
# Requirements:
#   dRep
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 INPUT_MAG_DIR OUTPUT_DIR" >&2
    exit 1
fi

input_dir="$1"
output_dir="$2"

threads="${SLURM_CPUS_PER_TASK:-8}"


# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

conda_profile="${CONDA_PROFILE:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
conda_environment="${DREP_CONDA_ENV:-binning-env}"

if [[ ! -f "$conda_profile" ]]; then
    echo "ERROR: Conda initialization script not found: $conda_profile" >&2
    exit 1
fi

source "$conda_profile"
conda activate "$conda_environment"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$input_dir" ]]; then
    echo "ERROR: Input MAG directory not found: $input_dir" >&2
    exit 1
fi

if ! command -v dRep >/dev/null 2>&1; then
    echo "ERROR: dRep is not available in the active environment." >&2
    exit 1
fi

mapfile -t genome_files < <(
    find "$input_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.fa' \
        -print \
        | sort
)

mag_count="${#genome_files[@]}"

if [[ "$mag_count" -eq 0 ]]; then
    echo "ERROR: No MAG files with the .fa extension found in: $input_dir" >&2
    exit 1
fi

if [[ -e "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "ERROR: Output directory already exists and is not empty: $output_dir" >&2
    echo "dRep should be run with a new or empty output directory." >&2
    exit 1
fi

mkdir -p "$output_dir"


# ---------------------------------------------------------------------------
# Run dRep
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting dRep dereplication"
echo "Input MAG directory: $input_dir"
echo "Number of MAGs:      $mag_count"
echo "Output directory:    $output_dir"
echo "Threads:             $threads"
echo "Minimum completeness: 70%"
echo "Maximum contamination: 10%"
echo "Secondary ANI:         0.99"
echo

dRep dereplicate "$output_dir" \
    -g "${genome_files[@]}" \
    -p "$threads" \
    -comp 70 \
    -con 10 \
    -sa 0.99


# ---------------------------------------------------------------------------
# Validate outputs
# ---------------------------------------------------------------------------

representative_dir="${output_dir}/dereplicated_genomes"
data_tables_dir="${output_dir}/data_tables"

if [[ ! -d "$representative_dir" ]]; then
    echo "ERROR: Expected dRep representative-genome directory was not created:" >&2
    echo "       $representative_dir" >&2
    exit 1
fi

representative_count=$(
    find "$representative_dir" \
        -maxdepth 1 \
        -type f \
        | wc -l
)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] dRep dereplication completed"
echo "Input MAGs:             $mag_count"
echo "Representative MAGs:    $representative_count"
echo "Representatives:        $representative_dir"
echo "Data tables:            $data_tables_dir"
echo "Runtime:                ${runtime_minutes} minutes"