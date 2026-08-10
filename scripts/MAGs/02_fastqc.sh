#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=fastqc_%j.log
#SBATCH --error=fastqc_%j.err

###############################################################################
# 02_fastqc.sh
#
# Run FastQC on one stage of the metagenomic workflow.
#
# Usage:
#   sbatch 02_fastqc.sh STAGE INPUT_DIR OUTPUT_DIR
#
# Supported stages:
#   raw       Raw paired-end reads stored in sample directories
#   trimmed   Trimmed reads stored in sample directories
#   filtered  Host-filtered reads stored in sample directories
#   merged    Merged host-filtered reads stored directly in INPUT_DIR
#
# Examples:
#   sbatch 02_fastqc.sh raw data/rawdata results/fastqc/raw
#
#   sbatch 02_fastqc.sh trimmed data/trimmed results/fastqc/trimmed
#
#   sbatch 02_fastqc.sh filtered \
#       data/filtered_reads \
#       results/fastqc/filtered
#
#   sbatch 02_fastqc.sh merged \
#       data/filtered_reads/mergedRhizo \
#       results/fastqc/mergedRhizo
#
# Requirements:
#   FastQC 0.11.9
###############################################################################

set -u
set -o pipefail
shopt -s nullglob


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

module load bio/FastQC/0.11.9-Java-1.8


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    echo "Usage: sbatch $0 STAGE INPUT_DIR OUTPUT_DIR" >&2
    echo "Supported stages: raw, trimmed, filtered, merged" >&2
    exit 1
fi

stage="$1"
input_dir="$2"
output_dir="$3"

threads="${SLURM_CPUS_PER_TASK:-4}"


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

case "$stage" in
    raw|trimmed|filtered|merged)
        ;;
    *)
        echo "ERROR: Unsupported stage: $stage" >&2
        echo "Supported stages: raw, trimmed, filtered, merged" >&2
        exit 1
        ;;
esac

if [[ ! -d "$input_dir" ]]; then
    echo "ERROR: Input directory not found: $input_dir" >&2
    exit 1
fi

if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not available after loading the module." >&2
    exit 1
fi

mkdir -p "$output_dir"


# ---------------------------------------------------------------------------
# Run information
# ---------------------------------------------------------------------------

start_time=$(date +%s)
processed_files=0

echo "FastQC stage:     $stage"
echo "Input directory:  $input_dir"
echo "Output directory: $output_dir"
echo "Threads:          $threads"
echo "Started:          $(date)"
echo


# ---------------------------------------------------------------------------
# Raw paired-end reads
#
# Expected input:
#   INPUT_DIR/SRR*/SRR*_1.fastq.gz
#   INPUT_DIR/SRR*/SRR*_2.fastq.gz
# ---------------------------------------------------------------------------

if [[ "$stage" == "raw" ]]; then

    sample_dirs=("$input_dir"/SRR*)

    if [[ ${#sample_dirs[@]} -eq 0 ]]; then
        echo "ERROR: No SRR sample directories found in: $input_dir" >&2
        exit 1
    fi

    for sample_dir in "${sample_dirs[@]}"; do
        [[ -d "$sample_dir" ]] || continue

        sample=$(basename "$sample_dir")
        sample_output_dir="$output_dir/$sample"

        fq1="$sample_dir/${sample}_1.fastq.gz"
        fq2="$sample_dir/${sample}_2.fastq.gz"

        if [[ -f "$fq1" && -f "$fq2" ]]; then
            mkdir -p "$sample_output_dir"

            echo "Running FastQC on raw sample: $sample"

            if fastqc \
                --threads "$threads" \
                --outdir "$sample_output_dir" \
                "$fq1" "$fq2"
            then
                processed_files=$((processed_files + 2))
            else
                echo "WARNING: FastQC failed for raw sample: $sample" >&2
            fi
        else
            echo "WARNING: Paired FASTQ files not found for: $sample" >&2
        fi
    done
fi


# ---------------------------------------------------------------------------
# Trimmed reads
#
# Expected input:
#   INPUT_DIR/SRR*/trimmed_*.fastq.gz
# ---------------------------------------------------------------------------

if [[ "$stage" == "trimmed" ]]; then

    sample_dirs=("$input_dir"/SRR*)

    if [[ ${#sample_dirs[@]} -eq 0 ]]; then
        echo "ERROR: No SRR sample directories found in: $input_dir" >&2
        exit 1
    fi

    for sample_dir in "${sample_dirs[@]}"; do
        [[ -d "$sample_dir" ]] || continue

        sample=$(basename "$sample_dir")
        sample_output_dir="$output_dir/$sample"
        fastq_files=("$sample_dir"/trimmed_*.fastq.gz)

        if [[ ${#fastq_files[@]} -eq 0 ]]; then
            echo "WARNING: No trimmed FASTQ files found for: $sample" >&2
            continue
        fi

        mkdir -p "$sample_output_dir"

        echo "Running FastQC on trimmed sample: $sample"

        if fastqc \
            --threads "$threads" \
            --outdir "$sample_output_dir" \
            "${fastq_files[@]}"
        then
            processed_files=$((processed_files + ${#fastq_files[@]}))
        else
            echo "WARNING: FastQC failed for trimmed sample: $sample" >&2
        fi
    done
fi


# ---------------------------------------------------------------------------
# Host-filtered reads
#
# Expected input:
#   INPUT_DIR/SRR*/non_host_R*_*.fastq.gz
# ---------------------------------------------------------------------------

if [[ "$stage" == "filtered" ]]; then

    sample_dirs=("$input_dir"/SRR*)

    if [[ ${#sample_dirs[@]} -eq 0 ]]; then
        echo "ERROR: No SRR sample directories found in: $input_dir" >&2
        exit 1
    fi

    for sample_dir in "${sample_dirs[@]}"; do
        [[ -d "$sample_dir" ]] || continue

        sample=$(basename "$sample_dir")
        sample_output_dir="$output_dir/$sample"
        fastq_files=("$sample_dir"/non_host_R*_*.fastq.gz)

        if [[ ${#fastq_files[@]} -eq 0 ]]; then
            echo "WARNING: No host-filtered FASTQ files found for: $sample" >&2
            continue
        fi

        mkdir -p "$sample_output_dir"

        echo "Running FastQC on filtered sample: $sample"

        if fastqc \
            --threads "$threads" \
            --outdir "$sample_output_dir" \
            "${fastq_files[@]}"
        then
            processed_files=$((processed_files + ${#fastq_files[@]}))
        else
            echo "WARNING: FastQC failed for filtered sample: $sample" >&2
        fi
    done
fi


# ---------------------------------------------------------------------------
# Merged host-filtered reads
#
# Expected input:
#   INPUT_DIR/merged_non_host_R*.fastq.gz
# ---------------------------------------------------------------------------

if [[ "$stage" == "merged" ]]; then

    fastq_files=("$input_dir"/merged_non_host_R*.fastq.gz)

    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        echo "ERROR: No merged host-filtered FASTQ files found in:" >&2
        echo "       $input_dir" >&2
        exit 1
    fi

    echo "Running FastQC on merged host-filtered reads."

    if fastqc \
        --threads "$threads" \
        --outdir "$output_dir" \
        "${fastq_files[@]}"
    then
        processed_files=${#fastq_files[@]}
    else
        echo "ERROR: FastQC failed for merged host-filtered reads." >&2
        exit 1
    fi
fi


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
duration_seconds=$((end_time - start_time))
duration_minutes=$((duration_seconds / 60))

echo
echo "FastQC completed."
echo "Stage:           $stage"
echo "Files processed: $processed_files"
echo "Runtime:         $duration_minutes minutes"
echo "Finished:        $(date)"

if [[ "$processed_files" -eq 0 ]]; then
    echo "ERROR: No FASTQ files were processed." >&2
    exit 1
fi