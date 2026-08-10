#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=96:00:00
#SBATCH --job-name=mapping
#SBATCH --mail-type=END,FAIL

###############################################################################
# 10_mapping.sh
#
# Map paired-end host-filtered reads to the filtered metagenomic assembly,
# sort the alignments, and index the resulting BAM file.
#
# Usage:
#   sbatch --array=0-N 10_mapping.sh \
#       SAMPLE_LIST FILTERED_READS_DIR CONTIGS OUTPUT_DIR
#
# Example:
#   mkdir -p logs/mapping
#
#   sbatch \
#       --array=0-11 \
#       --output=logs/mapping/mapping_%A_%a.log \
#       --error=logs/mapping/mapping_%A_%a.err \
#       10_mapping.sh \
#       config/rhizosphere_samples.txt \
#       data/filtered_reads \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa \
#       data/mappingRh
#
# Sample-list format:
#   One sample accession per line, without a header.
#
# Requirements:
#   BWA 0.7.17
#   SAMtools 1.19.2
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Software
# ---------------------------------------------------------------------------

module load bio/BWA/0.7.17-GCC-10.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch --array=0-N $0 SAMPLE_LIST FILTERED_READS_DIR CONTIGS OUTPUT_DIR" >&2
    exit 1
fi

sample_list="$1"
filtered_reads_dir="$2"
contigs="$3"
output_dir="$4"

threads="${SLURM_CPUS_PER_TASK:-8}"
task_id="${SLURM_ARRAY_TASK_ID:-}"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ -z "$task_id" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set." >&2
    echo "Submit this script using sbatch --array=0-N." >&2
    exit 1
fi

if [[ ! -f "$sample_list" ]]; then
    echo "ERROR: Sample list not found: $sample_list" >&2
    exit 1
fi

if [[ ! -f "$contigs" ]]; then
    echo "ERROR: Contig file not found: $contigs" >&2
    exit 1
fi

if ! command -v bwa >/dev/null 2>&1; then
    echo "ERROR: BWA is not available." >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: SAMtools is not available." >&2
    exit 1
fi

if ! command -v zcat >/dev/null 2>&1; then
    echo "ERROR: zcat is not available." >&2
    exit 1
fi

if [[ ! -f "${contigs}.bwt" ]]; then
    echo "ERROR: BWA index is missing for: $contigs" >&2
    echo "Run 09_bwa_index.sh before starting the mapping step." >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Select sample
# ---------------------------------------------------------------------------

mapfile -t samples < <(
    grep -v '^[[:space:]]*$' "$sample_list" |
    grep -v '^[[:space:]]*#'
)

if (( task_id >= ${#samples[@]} )); then
    echo "ERROR: Array task ID $task_id exceeds the sample-list range." >&2
    echo "Number of samples: ${#samples[@]}" >&2
    exit 1
fi

sample="${samples[$task_id]}"

read1="${filtered_reads_dir}/${sample}/non_host_R1_${sample}.fastq.gz"
read2="${filtered_reads_dir}/${sample}/non_host_R2_${sample}.fastq.gz"

sample_output_dir="${output_dir}/${sample}"
bam="${sample_output_dir}/${sample}.sorted.bam"

mkdir -p "$sample_output_dir"


# ---------------------------------------------------------------------------
# Validate reads
# ---------------------------------------------------------------------------

if [[ ! -f "$read1" ]]; then
    echo "ERROR: Forward-read file not found: $read1" >&2
    exit 1
fi

if [[ ! -f "$read2" ]]; then
    echo "ERROR: Reverse-read file not found: $read2" >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Mapping
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting mapping"
echo "Sample:           $sample"
echo "Forward reads:    $read1"
echo "Reverse reads:    $read2"
echo "Reference:        $contigs"
echo "Output BAM:       $bam"
echo "Threads:          $threads"
echo

bwa mem \
    -t "$threads" \
    "$contigs" \
    <(zcat "$read1") \
    <(zcat "$read2") |
samtools view \
    -bS \
    - |
samtools sort \
    -@ "$threads" \
    -o "$bam" \
    -

samtools index "$bam"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

properly_paired=$(samtools view -c -f 0x2 "$bam")

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Mapping completed."
echo "Sample:                         $sample"
echo "Properly paired mapped reads:   $properly_paired"
echo "Runtime:                        $runtime_minutes minutes"
echo "Sorted BAM:                     $bam"
echo "BAM index:                      ${bam}.bai"