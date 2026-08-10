#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=04:00:00

###############################################################################
# 03_1_trimmomatic.sh
#
# Perform paired-end quality trimming for one sample using Trimmomatic.
#
# Usage:
#   sbatch 03_1_trimmomatic.sh SAMPLE RAW_DIR OUTPUT_BASE [TIMING_LOG]
#
# Example:
#   sbatch 03_1_trimmomatic.sh \
#       SRR16095324 \
#       data/rawdata \
#       data/trimmed \
#       logs/trim/timing_summary.tsv
#
# Expected input:
#   RAW_DIR/SAMPLE/SAMPLE_1.fastq.gz
#   RAW_DIR/SAMPLE/SAMPLE_2.fastq.gz
#
# Output:
#   OUTPUT_BASE/SAMPLE/trimmed_R1_paired.fastq.gz
#   OUTPUT_BASE/SAMPLE/trimmed_R1_unpaired.fastq.gz
#   OUTPUT_BASE/SAMPLE/trimmed_R2_paired.fastq.gz
#   OUTPUT_BASE/SAMPLE/trimmed_R2_unpaired.fastq.gz
###############################################################################

set -u
set -o pipefail


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

module load bio/Trimmomatic/0.39-Java-1.8


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: sbatch $0 SAMPLE RAW_DIR OUTPUT_BASE [TIMING_LOG]" >&2
    exit 1
fi

sample="$1"
raw_dir="$2"
output_base="$3"
timing_log="${4:-}"

threads="${SLURM_CPUS_PER_TASK:-4}"


# ---------------------------------------------------------------------------
# Input and output paths
# ---------------------------------------------------------------------------

input_r1="${raw_dir}/${sample}/${sample}_1.fastq.gz"
input_r2="${raw_dir}/${sample}/${sample}_2.fastq.gz"

output_dir="${output_base}/${sample}"

output_r1_paired="${output_dir}/trimmed_R1_paired.fastq.gz"
output_r1_unpaired="${output_dir}/trimmed_R1_unpaired.fastq.gz"
output_r2_paired="${output_dir}/trimmed_R2_paired.fastq.gz"
output_r2_unpaired="${output_dir}/trimmed_R2_unpaired.fastq.gz"

trimmomatic_jar="${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar"
adapter_file="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE-2.fa"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$input_r1" ]]; then
    echo "ERROR: Forward read file not found: $input_r1" >&2
    exit 1
fi

if [[ ! -f "$input_r2" ]]; then
    echo "ERROR: Reverse read file not found: $input_r2" >&2
    exit 1
fi

if [[ ! -f "$trimmomatic_jar" ]]; then
    echo "ERROR: Trimmomatic JAR not found: $trimmomatic_jar" >&2
    exit 1
fi

if [[ ! -f "$adapter_file" ]]; then
    echo "ERROR: Adapter file not found: $adapter_file" >&2
    exit 1
fi

mkdir -p "$output_dir"

if [[ -n "$timing_log" ]]; then
    mkdir -p "$(dirname "$timing_log")"
fi


# ---------------------------------------------------------------------------
# Run Trimmomatic
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "Starting trimming for sample: $sample"
echo "Forward reads: $input_r1"
echo "Reverse reads: $input_r2"
echo "Output directory: $output_dir"
echo "Threads: $threads"
echo "Started: $(date)"

java -Xmx4G \
    -jar "$trimmomatic_jar" \
    PE \
    -threads "$threads" \
    "$input_r1" \
    "$input_r2" \
    "$output_r1_paired" \
    "$output_r1_unpaired" \
    "$output_r2_paired" \
    "$output_r2_unpaired" \
    "ILLUMINACLIP:${adapter_file}:2:30:10" \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50

exit_status=$?

if [[ $exit_status -ne 0 ]]; then
    echo "ERROR: Trimmomatic failed for sample: $sample" >&2
    exit "$exit_status"
fi


# ---------------------------------------------------------------------------
# Runtime summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))
job_id="${SLURM_JOB_ID:-not_slurm}"

if [[ -n "$timing_log" ]]; then
    if [[ ! -f "$timing_log" ]]; then
        printf "sample\truntime_minutes\tjob_id\n" > "$timing_log"
    fi

    printf "%s\t%s\t%s\n" \
        "$sample" \
        "$runtime_minutes" \
        "$job_id" >> "$timing_log"
fi

echo "Trimming completed for $sample in $runtime_minutes minutes."
echo "Job ID: $job_id"