#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --job-name=prodigal_drep_RhizoMAGs
#SBATCH --mail-type=END,FAIL

###############################################################################
# 19_prodigal_drep.sh
#
# Predict protein-coding genes in one dereplicated MAG per Slurm array task
# using Prodigal in metagenomic mode.
#
# Usage:
#   sbatch --array=1-N%50 19_prodigal_drep.sh \
#       GENOME_MANIFEST OUTPUT_DIR
#
# The manifest must contain one genome FASTA path per line.
#
# Requirements:
#   Prodigal 2.6.3
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch --array=1-N%50 $0 GENOME_MANIFEST OUTPUT_DIR" >&2
    exit 1
fi

genome_manifest="$1"
output_dir="$2"

task_id="${SLURM_ARRAY_TASK_ID:-}"


# ---------------------------------------------------------------------------
# Validate Slurm array task and inputs
# ---------------------------------------------------------------------------

if [[ -z "$task_id" ]]; then
    echo "ERROR: This script must be submitted as a Slurm array job." >&2
    exit 1
fi

if [[ ! -s "$genome_manifest" ]]; then
    echo "ERROR: Genome manifest not found or empty: $genome_manifest" >&2
    exit 1
fi

mkdir -p "$output_dir"


# ---------------------------------------------------------------------------
# Load Prodigal
# ---------------------------------------------------------------------------

module load bio/prodigal/2.6.3-GCCcore-10.2.0

if ! command -v prodigal >/dev/null 2>&1; then
    echo "ERROR: Prodigal is not available after loading the module." >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Select genome for this array task
# ---------------------------------------------------------------------------

genome_file="$(
    sed -n "${task_id}p" "$genome_manifest"
)"

if [[ -z "$genome_file" ]]; then
    echo "ERROR: No genome is assigned to array task ${task_id}." >&2
    exit 1
fi

if [[ ! -s "$genome_file" ]]; then
    echo "ERROR: Genome file not found or empty: $genome_file" >&2
    exit 1
fi

genome_filename="$(basename "$genome_file")"
genome_id="${genome_filename%.*}"

faa_file="${output_dir}/${genome_id}.faa"
ffn_file="${output_dir}/${genome_id}.ffn"
gff_file="${output_dir}/${genome_id}.gff"


# ---------------------------------------------------------------------------
# Run Prodigal
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting Prodigal"
echo "Array task:   $task_id"
echo "Genome ID:    $genome_id"
echo "Input genome: $genome_file"
echo "FAA output:   $faa_file"
echo "FFN output:   $ffn_file"
echo "GFF output:   $gff_file"
echo

if [[ -s "$faa_file" && -s "$ffn_file" && -s "$gff_file" ]]; then
    echo "[$(date)] All Prodigal outputs already exist; skipping ${genome_id}"
    exit 0
fi

prodigal \
    -i "$genome_file" \
    -a "$faa_file" \
    -d "$ffn_file" \
    -o "$gff_file" \
    -p meta \
    -q


# ---------------------------------------------------------------------------
# Validate outputs
# ---------------------------------------------------------------------------

for output_file in "$faa_file" "$ffn_file" "$gff_file"; do
    if [[ ! -s "$output_file" ]]; then
        echo "ERROR: Expected Prodigal output was not created or is empty:" >&2
        echo "       $output_file" >&2
        exit 1
    fi
done


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Prodigal completed for ${genome_id}"
echo "Runtime: ${runtime_minutes} minutes"