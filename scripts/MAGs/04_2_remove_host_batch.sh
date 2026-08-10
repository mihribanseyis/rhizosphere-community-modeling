#!/usr/bin/env bash

###############################################################################
# 04_2_remove_host_batch.sh
#
# Submit one host-removal Slurm job for each sample listed in an SRR file.
#
# Usage:
#   bash 04_2_remove_host_batch.sh \
#       SRR_LIST TRIMMED_BASE REFERENCE OUTPUT_BASE LOG_DIR [JOB_SCRIPT]
#
# Example:
#   bash 04_2_remove_host_batch.sh \
#       config/srr_ids.txt \
#       data/trimmed \
#       data/reference/maizeREF.fna \
#       data/filtered_reads \
#       logs/remove_host
#
###############################################################################

set -u
set -o pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: bash $0 SRR_LIST TRIMMED_BASE REFERENCE OUTPUT_BASE LOG_DIR [JOB_SCRIPT]" >&2
    exit 1
fi

srr_list="$1"
trimmed_base="$2"
reference="$3"
output_base="$4"
log_dir="$5"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
job_script="${6:-${script_dir}/04_remove_host.sh}"

summary_log="${log_dir}/remove_host_summary.tsv"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$srr_list" ]]; then
    echo "ERROR: SRR list not found: $srr_list" >&2
    exit 1
fi

if [[ ! -f "$reference" ]]; then
    echo "ERROR: Reference genome not found: $reference" >&2
    exit 1
fi

if [[ ! -f "$job_script" ]]; then
    echo "ERROR: Host-removal job script not found: $job_script" >&2
    exit 1
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available." >&2
    exit 1
fi

mkdir -p "$output_base"
mkdir -p "$log_dir"


# Create the shared summary file before submitting jobs. This avoids multiple
# jobs attempting to create the header at the same time.
if [[ ! -f "$summary_log" ]]; then
    printf "SRR_ID\tJob_ID\tDuration_min\tEnd_Time\n" > "$summary_log"
fi


# ---------------------------------------------------------------------------
# Submit one job per sample
# ---------------------------------------------------------------------------

submitted_jobs=0

while IFS= read -r sample || [[ -n "$sample" ]]; do

    # Remove carriage returns from Windows-formatted files.
    sample="${sample//$'\r'/}"

    # Ignore empty lines and comment lines.
    [[ -z "$sample" || "$sample" == \#* ]] && continue

    input_r1="${trimmed_base}/${sample}/trimmed_R1_paired.fastq.gz"
    input_r2="${trimmed_base}/${sample}/trimmed_R2_paired.fastq.gz"

    if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
        echo "WARNING: Paired trimmed reads not found for $sample; skipping." >&2
        continue
    fi

    echo "Submitting host-removal job for $sample..."

    if sbatch \
        --job-name="hostrm_${sample}" \
        --output="${log_dir}/bwa_euk_rm_${sample}_%j.log" \
        --error="${log_dir}/bwa_euk_rm_${sample}_%j.err" \
        "$job_script" \
        "$sample" \
        "$trimmed_base" \
        "$reference" \
        "$output_base" \
        "$summary_log"
    then
        submitted_jobs=$((submitted_jobs + 1))
    else
        echo "WARNING: Job submission failed for $sample." >&2
    fi

done < "$srr_list"

echo "Submitted $submitted_jobs host-removal jobs."