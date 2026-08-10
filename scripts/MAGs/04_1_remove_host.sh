#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=120:00:00

###############################################################################
# 04_1_remove_host.sh
#
# Align paired trimmed reads against the maize reference genome and retain
# read pairs for which both mates are unmapped.
#
# Usage:
#   sbatch 04_1_remove_host.sh \
#       SAMPLE TRIMMED_BASE REFERENCE OUTPUT_BASE [SUMMARY_LOG]
#
# Example:
#   sbatch 04_1_remove_host.sh \
#       SRR16095324 \
#       data/trimmed \
#       data/reference/maizeREF.fna \
#       data/filtered_reads \
#       logs/remove_host/remove_host_summary.tsv
#
# Expected input:
#   TRIMMED_BASE/SAMPLE/trimmed_R1_paired.fastq.gz
#   TRIMMED_BASE/SAMPLE/trimmed_R2_paired.fastq.gz
#
# Output:
#   OUTPUT_BASE/SAMPLE/sorted_SAMPLE.bam
#   OUTPUT_BASE/SAMPLE/non_host.bam
#   OUTPUT_BASE/SAMPLE/non_host_R1_SAMPLE.fastq.gz
#   OUTPUT_BASE/SAMPLE/non_host_R2_SAMPLE.fastq.gz
#
# Requirements:
#   BWA 0.7.17
#   SAMtools 1.19.2
#
# The reference genome must already be indexed with:
#   bwa index REFERENCE
###############################################################################

set -u
set -o pipefail


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

module load bio/BWA/0.7.17-GCC-10.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: sbatch $0 SAMPLE TRIMMED_BASE REFERENCE OUTPUT_BASE [SUMMARY_LOG]" >&2
    exit 1
fi

sample="$1"
trimmed_base="$2"
reference="$3"
output_base="$4"
summary_log="${5:-}"

threads="${SLURM_CPUS_PER_TASK:-4}"
job_id="${SLURM_JOB_ID:-not_slurm}"


# ---------------------------------------------------------------------------
# Input and output paths
# ---------------------------------------------------------------------------

trimmed_dir="${trimmed_base}/${sample}"
output_dir="${output_base}/${sample}"

input_r1="${trimmed_dir}/trimmed_R1_paired.fastq.gz"
input_r2="${trimmed_dir}/trimmed_R2_paired.fastq.gz"

sorted_bam="${output_dir}/sorted_${sample}.bam"
unmapped_bam="${output_dir}/non_host.bam"

fq1_out="${output_dir}/non_host_R1_${sample}.fastq.gz"
fq2_out="${output_dir}/non_host_R2_${sample}.fastq.gz"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$input_r1" ]]; then
    echo "ERROR: Forward trimmed reads not found: $input_r1" >&2
    exit 1
fi

if [[ ! -f "$input_r2" ]]; then
    echo "ERROR: Reverse trimmed reads not found: $input_r2" >&2
    exit 1
fi

if [[ ! -f "$reference" ]]; then
    echo "ERROR: Reference genome not found: $reference" >&2
    exit 1
fi

# BWA generates several index files. Checking one confirms that the expected
# index basename is present.
if [[ ! -f "${reference}.bwt" ]]; then
    echo "ERROR: BWA index not found for reference: $reference" >&2
    echo "Run: bwa index \"$reference\"" >&2
    exit 1
fi

for command in bwa samtools; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: Required command not available: $command" >&2
        exit 1
    fi
done

mkdir -p "$output_dir"

if [[ -n "$summary_log" ]]; then
    mkdir -p "$(dirname "$summary_log")"
fi


# ---------------------------------------------------------------------------
# Host alignment and removal
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "Starting host-read removal"
echo "Sample:           $sample"
echo "Reference:        $reference"
echo "Forward reads:    $input_r1"
echo "Reverse reads:    $input_r2"
echo "Output directory: $output_dir"
echo "Threads:          $threads"
echo "Started:          $(date)"
echo


# Align reads to the host genome, convert SAM output to BAM, and sort it.
if ! bwa mem \
    -t "$threads" \
    "$reference" \
    "$input_r1" \
    "$input_r2" \
    | samtools view -bS - \
    | samtools sort -o "$sorted_bam" -
then
    echo "ERROR: Alignment or BAM processing failed for $sample." >&2
    exit 1
fi


# Retain read pairs for which both the read and its mate are unmapped.
#
# SAM flag 12 combines:
#   4 = read unmapped
#   8 = mate unmapped
if ! samtools view \
    -b \
    -f 12 \
    "$sorted_bam" \
    -o "$unmapped_bam"
then
    echo "ERROR: Extraction of unmapped reads failed for $sample." >&2
    exit 1
fi


# Convert the unmapped paired reads back to FASTQ.
if ! samtools fastq \
    -1 "$fq1_out" \
    -2 "$fq2_out" \
    -0 /dev/null \
    -s /dev/null \
    -n \
    "$unmapped_bam"
then
    echo "ERROR: BAM-to-FASTQ conversion failed for $sample." >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Runtime summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
duration_minutes=$(( (end_time - start_time) / 60 ))
end_timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

echo
echo "Host-read removal completed."
echo "Sample:   $sample"
echo "Runtime:  $duration_minutes minutes"
echo "Job ID:   $job_id"
echo "Finished: $end_timestamp"

if [[ -n "$summary_log" ]]; then
    printf "%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$job_id" \
        "$duration_minutes" \
        "$end_timestamp" >> "$summary_log"
fi