#!/usr/bin/env bash

###############################################################################
# 03_2_trimmomatic_batch.sh
#
# Submit one Trimmomatic Slurm job per SRR accession.
#
# Usage:
#   bash 03_2_trimmomatic_batch.sh \
#       SRR_LIST RAW_DIR OUTPUT_BASE LOG_DIR [JOB_SCRIPT]
#
# Example:
#   bash 03_2_trimmomatic_batch.sh \
#       config/srr_ids.txt \
#       data/rawdata \
#       data/trimmed \
#       logs/trim
#
###############################################################################

set -u
set -o pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: bash $0 SRR_LIST RAW_DIR OUTPUT_BASE LOG_DIR [JOB_SCRIPT]" >&2
    exit 1
fi

srr_list="$1"
raw_dir="$2"
output_base="$3"
log_dir="$4"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
job_script="${5:-${script_dir}/03_trimmomatic.sh}"

timing_log="${log_dir}/timing_summary.tsv"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$srr_list" ]]; then
    echo "ERROR: SRR list not found: $srr_list" >&2
    exit 1
fi

if [[ ! -f "$job_script" ]]; then
    echo "ERROR: Trimmomatic job script not found: $job_script" >&2
    exit 1
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available." >&2
    exit 1
fi

mkdir -p "$output_base"
mkdir -p "$log_dir"


# ---------------------------------------------------------------------------
# Submit one job per sample
# ---------------------------------------------------------------------------

submitted_jobs=0

while IFS= read -r sample || [[ -n "$sample" ]]; do

    sample="${sample//$'\r'/}"

    # Ignore empty lines and comments.
    [[ -z "$sample" || "$sample" == \#* ]] && continue

    input_r1="${raw_dir}/${sample}/${sample}_1.fastq.gz"
    input_r2="${raw_dir}/${sample}/${sample}_2.fastq.gz"

    if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
        echo "WARNING: Paired input files not found for $sample; job not submitted." >&2
        continue
    fi

    echo "Submitting trimming job for $sample..."

    sbatch \
        --job-name="trim_${sample}" \
        --output="${log_dir}/trimmomatic_${sample}_%j.log" \
        --error="${log_dir}/trimmomatic_${sample}_%j.err" \
        "$job_script" \
        "$sample" \
        "$raw_dir" \
        "$output_base" \
        "$timing_log"

    if [[ $? -eq 0 ]]; then
        submitted_jobs=$((submitted_jobs + 1))
    else
        echo "WARNING: Job submission failed for $sample." >&2
    fi

done < "$srr_list"

echo "Submitted $submitted_jobs trimming jobs."