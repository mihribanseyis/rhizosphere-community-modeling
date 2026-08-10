#!/usr/bin/env bash

###############################################################################
# 01_download.sh
#
# Purpose:
#   Download sequencing runs from the NCBI Sequence Read Archive using
#   SRA Toolkit, convert them to paired FASTQ files, compress the FASTQ files,
#   and remove the downloaded SRA files to reduce storage usage.
#
# Requirements:
#   - SRA Toolkit: prefetch and fasterq-dump
#   - gzip
#
# Usage:
#   bash 01_download.sh SRR_list.txt OUTPUT_DIR [LOG_DIR] [THREADS]
#
# Example:
#   bash 01_download.sh \
#       config/srr_batches/batch_01.txt \
#       data/rawdata \
#       logs \
#       4
#
# Arguments:
#   1. Text file containing one SRR accession per line
#   2. Base output directory for downloaded reads
#   3. Log directory (optional; default: logs)
#   4. Threads per fasterq-dump process (optional; default: 4)
#
# Output structure:
#   OUTPUT_DIR/
#   ├── SRRXXXXXXXX/
#   │   ├── SRRXXXXXXXX_1.fastq.gz
#   │   └── SRRXXXXXXXX_2.fastq.gz
#   └── ...
###############################################################################

set -u
set -o pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: bash $0 SRR_LIST OUTPUT_DIR [LOG_DIR] [THREADS]" >&2
    exit 1
fi

srr_list="$1"
base_output_dir="$2"
log_dir="${3:-logs}"
threads="${4:-4}"


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

if [[ ! -f "$srr_list" ]]; then
    echo "ERROR: SRR list not found: $srr_list" >&2
    exit 1
fi

if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 1
fi

for command in prefetch fasterq-dump gzip; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $command" >&2
        exit 1
    fi
done


# ---------------------------------------------------------------------------
# Output and log directories
# ---------------------------------------------------------------------------

mkdir -p "$base_output_dir"
mkdir -p "$log_dir"

batch_name="$(basename "$srr_list")"
batch_name="${batch_name%.*}"

logfile="${log_dir}/download_log_${batch_name}.txt"

{
    echo "============================================================"
    echo "Download batch started: $(date)"
    echo "SRR list: $srr_list"
    echo "Output directory: $base_output_dir"
    echo "Threads: $threads"
    echo "============================================================"
} >> "$logfile"


# ---------------------------------------------------------------------------
# Download and convert each sequencing run
# ---------------------------------------------------------------------------

while IFS= read -r srr || [[ -n "$srr" ]]; do

    # Remove carriage returns that may occur in Windows-formatted text files.
    srr="${srr//$'\r'/}"

    # Ignore empty lines and comment lines.
    [[ -z "$srr" || "$srr" == \#* ]] && continue

    echo "========== Starting $srr =========="
    echo "$(date) - START $srr" >> "$logfile"

    start_time=$(date +%s)

    sample_dir="${base_output_dir}/${srr}"
    mkdir -p "$sample_dir"


    # Download the SRA record.
    if ! prefetch "$srr" -O "$sample_dir"; then
        echo "ERROR: prefetch failed for $srr" | tee -a "$logfile"
        continue
    fi


    # Convert the downloaded SRA record to paired FASTQ files.
    if ! fasterq-dump \
        --split-files \
        --threads "$threads" \
        --outdir "$sample_dir" \
        "$sample_dir/$srr"
    then
        echo "ERROR: fasterq-dump failed for $srr" | tee -a "$logfile"
        continue
    fi


    # Compress the generated FASTQ files.
    fastq_files=("$sample_dir"/*.fastq)

    if [[ ! -e "${fastq_files[0]}" ]]; then
        echo "ERROR: No FASTQ files were produced for $srr" | tee -a "$logfile"
        continue
    fi

    if ! gzip "${fastq_files[@]}"; then
        echo "ERROR: FASTQ compression failed for $srr" | tee -a "$logfile"
        continue
    fi


    # Remove the downloaded SRA directory only after successful conversion
    # and compression.
    rm -rf "$sample_dir/$srr"


    end_time=$(date +%s)
    runtime=$(( (end_time - start_time) / 60 ))

    echo "Sample $srr completed in $runtime minutes."
    echo "$(date) - DONE $srr in $runtime min" >> "$logfile"
    echo "========== Completed $srr =========="

done < "$srr_list"

echo "$(date) - Download batch completed" >> "$logfile"