#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=4-00:00:00
#SBATCH --job-name=concoct_global
#SBATCH --mail-type=END,FAIL

###############################################################################
# 12_3_concoct.sh
#
# Perform global CONCOCT binning using the shared rhizosphere assembly and
# BAM files from all rhizosphere samples.
#
# Usage:
#   sbatch 12_3_concoct.sh \
#       SAMPLE_LIST CONTIGS MAPPING_DIR OUTPUT_DIR
#
# Requirements:
#   CONCOCT
#   cut_up_fasta.py
#   concoct_coverage_table.py
#   merge_cutup_clustering.py
#   extract_fasta_bins.py
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 SAMPLE_LIST CONTIGS MAPPING_DIR OUTPUT_DIR" >&2
    exit 1
fi

sample_list="$1"
contigs="$2"
mapping_dir="$3"
output_dir="$4"

work_dir="${output_dir}/tmp"
run_dir="${output_dir}/run"
bins_dir="${run_dir}/bins"
threads="${SLURM_CPUS_PER_TASK:-8}"


# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

conda_profile="${CONDA_PROFILE:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

if [[ ! -f "$conda_profile" ]]; then
    echo "ERROR: Conda initialization script not found: $conda_profile" >&2
    exit 1
fi

source "$conda_profile"
conda activate binning-env


# ---------------------------------------------------------------------------
# Validate inputs and software
# ---------------------------------------------------------------------------

if [[ ! -f "$sample_list" ]]; then
    echo "ERROR: Sample list not found: $sample_list" >&2
    exit 1
fi

if [[ ! -s "$contigs" ]]; then
    echo "ERROR: Contig file not found or empty: $contigs" >&2
    exit 1
fi

if [[ ! -d "$mapping_dir" ]]; then
    echo "ERROR: Mapping directory not found: $mapping_dir" >&2
    exit 1
fi

required_programs=(
    cut_up_fasta.py
    concoct_coverage_table.py
    concoct
    merge_cutup_clustering.py
    extract_fasta_bins.py
)

for program in "${required_programs[@]}"; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

mkdir -p "$work_dir" "$run_dir" "$bins_dir"


# ---------------------------------------------------------------------------
# Read sample IDs and collect BAM files
# ---------------------------------------------------------------------------

mapfile -t samples < <(
    grep -v '^[[:space:]]*$' "$sample_list" |
    grep -v '^[[:space:]]*#'
)

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "ERROR: No sample IDs found in: $sample_list" >&2
    exit 1
fi

bam_paths=()

for sample in "${samples[@]}"; do
    bam="${mapping_dir}/${sample}/${sample}.sorted.bam"

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM file not found or empty: $bam" >&2
        exit 1
    fi

    bam_paths+=("$bam")
done


# ---------------------------------------------------------------------------
# Run CONCOCT workflow
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting global CONCOCT binning"
echo "Sample list:       $sample_list"
echo "Number of samples: ${#samples[@]}"
echo "Contigs:           $contigs"
echo "Mapping directory: $mapping_dir"
echo "Output directory:  $output_dir"
echo "Threads:           $threads"
echo
echo "BAM files:"
printf '  %s\n' "${bam_paths[@]}"
echo

echo "[$(date)] Cutting contigs into 10 kb fragments"

cut_up_fasta.py "$contigs" \
    -c 10000 \
    -o 0 \
    --merge_last \
    -b "${work_dir}/contigs_10k.bed" \
    > "${work_dir}/contigs_10k.fa"

echo "[$(date)] Computing CONCOCT coverage table"

concoct_coverage_table.py \
    "${work_dir}/contigs_10k.bed" \
    "${bam_paths[@]}" \
    > "${work_dir}/coverage.tsv"

echo "[$(date)] Running CONCOCT"

concoct \
    --composition_file "${work_dir}/contigs_10k.fa" \
    --coverage_file "${work_dir}/coverage.tsv" \
    -b "${run_dir}/" \
    -t "$threads"

echo "[$(date)] Merging cut-up contig clusters"

merge_cutup_clustering.py \
    "${run_dir}/clustering_gt1000.csv" \
    > "${run_dir}/clustering_merged.csv"

echo "[$(date)] Extracting CONCOCT bins"

extract_fasta_bins.py \
    "$contigs" \
    "${run_dir}/clustering_merged.csv" \
    --output_path "$bins_dir"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Global CONCOCT binning completed"
echo "Output bins: $bins_dir"
echo "Runtime:     ${runtime_minutes} minutes"