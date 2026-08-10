#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --job-name=merge_filtered

###############################################################################
# 06_merge.sh
#
# Concatenate host-filtered paired-end FASTQ files across samples.
#
# Usage:
#   sbatch 06_merge.sh FILTERED_BASE MERGED_DIR
#
# Example:
#   sbatch \
#       --output=logs/merge/merge_filtered_%j.log \
#       --error=logs/merge/merge_filtered_%j.err \
#       06_merge.sh \
#       data/filtered_reads \
#       data/filtered_reads/mergedRhizo
#
# Expected input:
#   FILTERED_BASE/SRR*/non_host_R1_SRR*.fastq.gz
#   FILTERED_BASE/SRR*/non_host_R2_SRR*.fastq.gz
#
# Outputs:
#   MERGED_DIR/mergedRh_non_host_R1.fastq.gz
#   MERGED_DIR/mergedRh_non_host_R2.fastq.gz
###############################################################################

set -euo pipefail
shopt -s nullglob


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 FILTERED_BASE MERGED_DIR" >&2
    exit 1
fi

filtered_base="$1"
merged_dir="$2"

merged_r1="${merged_dir}/mergedRh_non_host_R1.fastq.gz"
merged_r2="${merged_dir}/mergedRh_non_host_R2.fastq.gz"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$filtered_base" ]]; then
    echo "ERROR: Filtered-read directory not found: $filtered_base" >&2
    exit 1
fi

for command in gzip wc cat; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $command" >&2
        exit 1
    fi
done

sample_dirs=("$filtered_base"/SRR*)

if [[ ${#sample_dirs[@]} -eq 0 ]]; then
    echo "ERROR: No SRR sample directories found in: $filtered_base" >&2
    exit 1
fi

mkdir -p "$merged_dir"


# ---------------------------------------------------------------------------
# Merge reads
# ---------------------------------------------------------------------------

echo "[$(date)] Starting merge of filtered paired-end reads."

start_time=$(date +%s)

# Truncate existing output files before appending.
: > "$merged_r1"
: > "$merged_r2"

total_lines_r1=0
total_lines_r2=0
sample_count=0

for sample_dir in "${sample_dirs[@]}"; do
    [[ -d "$sample_dir" ]] || continue

    sample="$(basename "$sample_dir")"

    fq1="${sample_dir}/non_host_R1_${sample}.fastq.gz"
    fq2="${sample_dir}/non_host_R2_${sample}.fastq.gz"

    if [[ ! -s "$fq1" || ! -s "$fq2" ]]; then
        echo "WARNING: Missing or empty paired FASTQ files for $sample; skipping." >&2
        continue
    fi

    echo "Processing $sample"

    lines_r1=$(gzip -cd "$fq1" | wc -l)
    lines_r2=$(gzip -cd "$fq2" | wc -l)

    if (( lines_r1 % 4 != 0 || lines_r2 % 4 != 0 )); then
        echo "ERROR: FASTQ line count is not divisible by four for $sample." >&2
        exit 1
    fi

    reads_r1=$((lines_r1 / 4))
    reads_r2=$((lines_r2 / 4))

    if [[ "$reads_r1" -ne "$reads_r2" ]]; then
        echo "ERROR: Paired read counts differ for $sample:" >&2
        echo "       R1=$reads_r1, R2=$reads_r2" >&2
        exit 1
    fi

    total_lines_r1=$((total_lines_r1 + lines_r1))
    total_lines_r2=$((total_lines_r2 + lines_r2))

    # Concatenating gzip files is valid: decompression tools read the
    # concatenated gzip members as one continuous stream.
    cat "$fq1" >> "$merged_r1"
    cat "$fq2" >> "$merged_r2"

    sample_count=$((sample_count + 1))
done


# ---------------------------------------------------------------------------
# Post-merge checks
# ---------------------------------------------------------------------------

if [[ "$sample_count" -eq 0 ]]; then
    echo "ERROR: No samples were merged." >&2
    rm -f "$merged_r1" "$merged_r2"
    exit 1
fi

merged_lines_r1=$(gzip -cd "$merged_r1" | wc -l)
merged_lines_r2=$(gzip -cd "$merged_r2" | wc -l)

if (( merged_lines_r1 % 4 != 0 || merged_lines_r2 % 4 != 0 )); then
    echo "ERROR: Merged FASTQ line count is not divisible by four." >&2
    exit 1
fi

merged_reads_r1=$((merged_lines_r1 / 4))
merged_reads_r2=$((merged_lines_r2 / 4))

if [[ "$merged_lines_r1" -ne "$total_lines_r1" ]]; then
    echo "ERROR: Merged R1 line count does not match the input total." >&2
    exit 1
fi

if [[ "$merged_lines_r2" -ne "$total_lines_r2" ]]; then
    echo "ERROR: Merged R2 line count does not match the input total." >&2
    exit 1
fi

if [[ "$merged_reads_r1" -ne "$merged_reads_r2" ]]; then
    echo "ERROR: Merged R1 and R2 read counts differ." >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
duration_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "Merge completed."
echo "Samples merged: $sample_count"
echo "Merged R1 reads: $merged_reads_r1"
echo "Merged R2 reads: $merged_reads_r2"
echo "Output R1:       $merged_r1"
echo "Output R2:       $merged_r2"
echo
ls -lh "$merged_r1" "$merged_r2"
echo
echo "Runtime: $duration_minutes minutes"
echo "Finished: $(date)"