#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=5-00:00:00
#SBATCH --job-name=maxbin_group
#SBATCH --mail-type=END,FAIL

###############################################################################
# 12_1_maxbin_group.sh
#
# Perform group-wise MaxBin2 binning using treatment-specific sample depth
# files generated from the shared rhizosphere co-assembly.
#
# Usage:
#   sbatch 12_1_maxbin_group.sh \
#       GROUP SAMPLE_LIST ASSEMBLY MAPPING_DIR OUTPUT_BASE
#
# Example:
#   sbatch \
#       12_1_maxbin_group.sh \
#       CK \
#       config/CK_samples.txt \
#       data/megahit_assembly/mergedRh/contigs_2kb.fa \
#       data/mappingRh \
#       data/binningRh/maxbin
#
# Requirements:
#   MaxBin2 2.2.7
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

if [[ $# -ne 5 ]]; then
    echo "Usage: sbatch $0 GROUP SAMPLE_LIST ASSEMBLY MAPPING_DIR OUTPUT_BASE" >&2
    exit 1
fi

group="$1"
sample_list="$2"
assembly="$3"
mapping_dir="$4"
output_base="$5"

threads="${SLURM_CPUS_PER_TASK:-8}"
output_dir="${output_base}/${group}"
output_prefix="${output_dir}/${group}_maxbin"


# ---------------------------------------------------------------------------
# Validate group name
# ---------------------------------------------------------------------------

case "$group" in
    CK|NP|NPM)
        ;;
    *)
        echo "ERROR: Unknown group: $group" >&2
        echo "Allowed groups: CK, NP, NPM" >&2
        exit 1
        ;;
esac


# ---------------------------------------------------------------------------
# Validate inputs and software
# ---------------------------------------------------------------------------

if [[ ! -f "$sample_list" ]]; then
    echo "ERROR: Sample list not found: $sample_list" >&2
    exit 1
fi

if [[ ! -s "$assembly" ]]; then
    echo "ERROR: Assembly file not found or empty: $assembly" >&2
    exit 1
fi

if [[ ! -d "$mapping_dir" ]]; then
    echo "ERROR: Mapping directory not found: $mapping_dir" >&2
    exit 1
fi

if ! command -v run_MaxBin.pl >/dev/null 2>&1; then
    echo "ERROR: run_MaxBin.pl is not available in the active environment." >&2
    exit 1
fi

mkdir -p "$output_dir"


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
# Build MaxBin2 abundance arguments
# ---------------------------------------------------------------------------

abundance_args=()

for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    depth_file="${mapping_dir}/${sample}/${sample}_depth.txt"

    if [[ ! -s "$depth_file" ]]; then
        echo "ERROR: Depth file not found or empty: $depth_file" >&2
        exit 1
    fi

    abundance_args+=("-abund${i}" "$depth_file")
done


# ---------------------------------------------------------------------------
# Run MaxBin2
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting MaxBin2 group-wise binning"
echo "Group:              $group"
echo "Sample list:        $sample_list"
echo "Number of samples:  ${#samples[@]}"
echo "Assembly:           $assembly"
echo "Mapping directory:  $mapping_dir"
echo "Output prefix:      $output_prefix"
echo "Minimum contig:     2000 bp"
echo "Threads:            $threads"
echo
echo "Depth files:"

for sample in "${samples[@]}"; do
    echo "  ${mapping_dir}/${sample}/${sample}_depth.txt"
done

echo

run_MaxBin.pl \
    -contig "$assembly" \
    "${abundance_args[@]}" \
    -out "$output_prefix" \
    -min_contig_length 2000 \
    -thread "$threads"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] MaxBin2 group-wise binning completed."
echo "Group:         $group"
echo "Output prefix: $output_prefix"
echo "Runtime:       ${runtime_minutes} minutes"