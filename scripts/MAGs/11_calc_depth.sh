#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=depth_generation
#SBATCH --mail-type=END,FAIL

###############################################################################
# 11_calc_depth.sh
#
# Generate:
#   1. One MetaBAT2-compatible depth file for each sample
#   2. One combined depth file containing all samples
#
# Usage:
#   sbatch 11_calc_depth.sh \
#       SAMPLE_LIST MAPPED_DIR CONTIGS OUTPUT_DEPTH_FILE
#
# Example:
#   sbatch \
#       --output=logs/depth_%j.log \
#       --error=logs/depth_%j.err \
#       11_calc_depth.sh \
#       config/rhizosphere_samples.txt \
#       data/mappingRh \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa \
#       data/mappingRh/all_samples_combined_depth.txt
#
# Requirements:
#   jgi_summarize_bam_contig_depths
#   gzip
#   awk
#
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

CONDA_PROFILE="${CONDA_PROFILE:-$HOME/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="binning-env"

if [[ ! -f "$CONDA_PROFILE" ]]; then
    echo "ERROR: Cannot find conda.sh: $CONDA_PROFILE" >&2
    exit 1
fi

source "$CONDA_PROFILE"
conda activate "$CONDA_ENV"


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 SAMPLE_LIST MAPPED_DIR CONTIGS OUTPUT_DEPTH_FILE" >&2
    exit 1
fi

sample_list="$1"
mapped_dir="$2"
contigs="$3"
combined_depth_file="$4"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -f "$sample_list" ]]; then
    echo "ERROR: Sample list not found: $sample_list" >&2
    exit 1
fi

if [[ ! -d "$mapped_dir" ]]; then
    echo "ERROR: Mapping directory not found: $mapped_dir" >&2
    exit 1
fi

if [[ ! -f "$contigs" ]]; then
    echo "ERROR: Contig file not found: $contigs" >&2
    exit 1
fi

for program in jgi_summarize_bam_contig_depths gzip awk; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "$combined_depth_file")"


# ---------------------------------------------------------------------------
# Read sample IDs
# ---------------------------------------------------------------------------

mapfile -t samples < <(
    grep -v '^[[:space:]]*$' "$sample_list" |
    grep -v '^[[:space:]]*#'
)

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "ERROR: No sample IDs found in: $sample_list" >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Check BAM files and collect paths
# ---------------------------------------------------------------------------

bam_files=()

for sample in "${samples[@]}"; do
    sample_dir="${mapped_dir}/${sample}"
    bam="${sample_dir}/${sample}.sorted.bam"

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM file missing or empty for $sample: $bam" >&2
        exit 1
    fi

    bam_files+=("$bam")
done


# ---------------------------------------------------------------------------
# Generate individual depth files
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting depth-file generation"
echo "Sample list:          $sample_list"
echo "Number of samples:    ${#samples[@]}"
echo "Mapped BAM directory: $mapped_dir"
echo "Reference contigs:    $contigs"
echo "Combined output:      $combined_depth_file"
echo

for sample in "${samples[@]}"; do
    sample_dir="${mapped_dir}/${sample}"
    bam="${sample_dir}/${sample}.sorted.bam"
    depth_file="${sample_dir}/${sample}_depth.txt"

    echo "[$(date)] Processing $sample"

    jgi_summarize_bam_contig_depths \
        --outputDepth "$depth_file" \
        --referenceFasta "$contigs" \
        "$bam"

    gzip -c "$depth_file" > "${depth_file}.gz"

    negative_count=$(
        awk '
            NR > 1 {
                for (i = 2; i <= NF; i++) {
                    if ($i < 0) {
                        count++
                        break
                    }
                }
            }
            END {
                print count + 0
            }
        ' "$depth_file"
    )

    if [[ "$negative_count" -gt 0 ]]; then
        echo "WARNING: $sample has $negative_count contigs with negative depth values." >&2
    else
        echo "$sample: no negative depth values found."
    fi
done


# ---------------------------------------------------------------------------
# Generate combined depth file
# ---------------------------------------------------------------------------

echo
echo "[$(date)] Generating combined depth file"

jgi_summarize_bam_contig_depths \
    --outputDepth "$combined_depth_file" \
    --referenceFasta "$contigs" \
    "${bam_files[@]}"

gzip -c "$combined_depth_file" > "${combined_depth_file}.gz"


# ---------------------------------------------------------------------------
# Check combined depth file
# ---------------------------------------------------------------------------

negative_count_all=$(
    awk '
        NR > 1 {
            for (i = 2; i <= NF; i++) {
                if ($i < 0) {
                    count++
                    break
                }
            }
        }
        END {
            print count + 0
        }
    ' "$combined_depth_file"
)

if [[ "$negative_count_all" -gt 0 ]]; then
    echo "WARNING: Combined depth file has $negative_count_all contigs with negative values." >&2
else
    echo "Combined depth file: no negative depth values found."
fi


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Depth-file generation completed"
echo "Combined depth file: ${combined_depth_file}"
echo "Compressed copy:     ${combined_depth_file}.gz"
echo "Runtime:             ${runtime_minutes} minutes"