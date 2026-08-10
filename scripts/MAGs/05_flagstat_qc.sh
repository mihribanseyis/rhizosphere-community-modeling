#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --job-name=flagstat_qc

###############################################################################
# 05_flagstat_qc.sh
#
# Run SAMtools flagstat on sorted BAM files and generate a combined mapping
# quality summary.
#
# Usage:
#   sbatch 05_flagstat_qc.sh MAPPED_DIR [SUMMARY_FILE]
#
# Example:
#   sbatch \
#       --output=logs/mapping/flagstat_%j.log \
#       --error=logs/mapping/flagstat_%j.err \
#       05_flagstat_qc.sh \
#       data/mappingRh \
#       results/mapping_qc_summary.tsv
#
# Expected input structure:
#   MAPPED_DIR/
#   ├── SRRXXXXXXXX/
#   │   └── SRRXXXXXXXX.sorted.bam
#   └── ...
#
# Outputs:
#   - One flagstat report beside each BAM file:
#       SRRXXXXXXXX.sorted.bam.flagstat.txt
#   - One combined tab-separated summary file
###############################################################################

set -euo pipefail
shopt -s nullglob


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

module load bio/SAMtools/1.19.2-GCC-13.2.0


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch $0 MAPPED_DIR [SUMMARY_FILE]" >&2
    exit 1
fi

mapped_dir="$1"
summary_file="${2:-${mapped_dir}/mapping_qc_summary.tsv}"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$mapped_dir" ]]; then
    echo "ERROR: Mapping directory not found: $mapped_dir" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: SAMtools is not available after loading the module." >&2
    exit 1
fi

mkdir -p "$(dirname "$summary_file")"

bam_files=("$mapped_dir"/SRR*/*.sorted.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "ERROR: No sorted BAM files found under:" >&2
    echo "       $mapped_dir/SRR*/*.sorted.bam" >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Run flagstat
# ---------------------------------------------------------------------------

start_time=$(date +%s)

printf "Sample\tTotal_Reads\tMapped_Reads\tMapped_%%\tProperly_Paired\tProperly_Paired_%%\n" \
    > "$summary_file"

processed_bams=0

for bam_file in "${bam_files[@]}"; do

    sample="$(basename "$bam_file" .sorted.bam)"
    flagstat_file="${bam_file}.flagstat.txt"

    echo "[$(date)] Running flagstat for: $sample"

    samtools flagstat \
        --threads "${SLURM_CPUS_PER_TASK:-2}" \
        "$bam_file" > "$flagstat_file"

    # Extract values from the standard SAMtools flagstat output.
    total="$(
        awk '/in total/ {
            print $1
            exit
        }' "$flagstat_file"
    )"

    mapped="$(
        awk '$4 == "mapped" {
            print $1
            exit
        }' "$flagstat_file"
    )"

    mapped_pct="$(
        awk '$4 == "mapped" {
            gsub(/[()%]/, "", $5)
            print $5
            exit
        }' "$flagstat_file"
    )"

    properly_paired="$(
        awk '$4 == "properly" && $5 == "paired" {
            print $1
            exit
        }' "$flagstat_file"
    )"

    properly_paired_pct="$(
        awk '$4 == "properly" && $5 == "paired" {
            gsub(/[()%]/, "", $6)
            print $6
            exit
        }' "$flagstat_file"
    )"

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$total" \
        "$mapped" \
        "$mapped_pct" \
        "$properly_paired" \
        "$properly_paired_pct" \
        >> "$summary_file"

    processed_bams=$((processed_bams + 1))

done


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "Mapping QC completed."
echo "BAM files processed: $processed_bams"
echo "Summary file:        $summary_file"
echo "Runtime:             $runtime_minutes minutes"