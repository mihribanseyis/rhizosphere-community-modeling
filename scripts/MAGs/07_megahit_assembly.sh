#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000G
#SBATCH --time=120:00:00
#SBATCH --job-name=megahit_assembly
#SBATCH --mail-type=END,FAIL

###############################################################################
# 07_megahit_assembly.sh
#
# Perform paired-end metagenomic co-assembly using MEGAHIT.
#
# Usage:
#   sbatch 07_megahit_assembly.sh READ1 READ2 OUTPUT_DIR
#
# Example:
#   sbatch \
#       --output=logs/assembly/megahit_%j.log \
#       --error=logs/assembly/megahit_%j.err \
#       07_megahit_assembly.sh \
#       data/filtered_reads/mergedRhizo/mergedRh_non_host_R1.fastq.gz \
#       data/filtered_reads/mergedRhizo/mergedRh_non_host_R2.fastq.gz \
#       data/megahit_assembly/mergedRh
#
# Requirements:
#   MEGAHIT 1.2.9
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate megahit_env


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 READ1 READ2 OUTPUT_DIR" >&2
    exit 1
fi

read1="$1"
read2="$2"
output_dir="$3"

threads="${SLURM_CPUS_PER_TASK:-32}"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$read1" ]]; then
    echo "ERROR: Forward-read file not found: $read1" >&2
    exit 1
fi

if [[ ! -f "$read2" ]]; then
    echo "ERROR: Reverse-read file not found: $read2" >&2
    exit 1
fi

if ! command -v megahit >/dev/null 2>&1; then
    echo "ERROR: MEGAHIT is not available in the active environment." >&2
    exit 1
fi

mkdir -p "$(dirname "$output_dir")"


# ---------------------------------------------------------------------------
# Run MEGAHIT
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting MEGAHIT assembly"
echo "Forward reads:    $read1"
echo "Reverse reads:    $read2"
echo "Output directory: $output_dir"
echo "Preset:           meta-large"
echo "Minimum contig:   200 bp"
echo "Threads:          $threads"
echo

megahit \
    -1 "$read1" \
    -2 "$read2" \
    -o "$output_dir" \
    --presets meta-large \
    --min-contig-len 200 \
    --num-cpu-threads "$threads" \
    --continue


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] MEGAHIT assembly completed."
echo "Runtime: $runtime_minutes minutes"
echo "Output:  $output_dir"

conda deactivate