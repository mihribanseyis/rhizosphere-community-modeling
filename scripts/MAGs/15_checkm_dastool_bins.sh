#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --job-name=checkm_RhizoMAGs
#SBATCH --mail-type=END,FAIL

###############################################################################
# 15_checkm_dastool_bins.sh
#
# Rename DAS Tool-refined bins and evaluate their completeness and
# contamination using CheckM.
#
# Usage:
#   sbatch 15_checkm_dastool_bins.sh DASTOOL_DIR
#
# Requirements:
#   CheckM 1.2.3
###############################################################################

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 DASTOOL_DIR" >&2
    exit 1
fi

dastool_dir="$1"

dastool_bins="${dastool_dir}/allRh_DASTool_bins"
renamed_bins="${dastool_dir}/RhizoMAGs"
checkm_output="${dastool_dir}/checkm_RhizoMAGs"
mapping_file="${renamed_bins}/RhizoMAG_name_mapping.tsv"

threads="${SLURM_CPUS_PER_TASK:-8}"


# ---------------------------------------------------------------------------
# Conda environment
# ---------------------------------------------------------------------------

conda_profile="${CONDA_PROFILE:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
conda_environment="${CHECKM_CONDA_ENV:-base}"

if [[ ! -f "$conda_profile" ]]; then
    echo "ERROR: Conda initialization script not found: $conda_profile" >&2
    exit 1
fi

source "$conda_profile"
conda activate "$conda_environment"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$dastool_bins" ]]; then
    echo "ERROR: DAS Tool bin directory not found: $dastool_bins" >&2
    exit 1
fi

if ! command -v checkm >/dev/null 2>&1; then
    echo "ERROR: checkm is not available in the active environment." >&2
    exit 1
fi

if ! compgen -G "${dastool_bins}/*.fa" >/dev/null; then
    echo "ERROR: No DAS Tool bin files found in: $dastool_bins" >&2
    exit 1
fi

mkdir -p "$renamed_bins" "$checkm_output"


# ---------------------------------------------------------------------------
# Rename DAS Tool bins
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting CheckM evaluation of DAS Tool bins"
echo "DAS Tool bins:  $dastool_bins"
echo "Renamed bins:   $renamed_bins"
echo "CheckM output:  $checkm_output"
echo "Threads:        $threads"
echo

find "$renamed_bins" \
    -maxdepth 1 \
    -type f \
    -name 'RhizoMAG_*.fa' \
    -delete

printf "new_name\toriginal_name\n" > "$mapping_file"

bin_index=1

while IFS= read -r bin_file; do
    original_name="$(basename "$bin_file")"
    new_filename="$(printf 'RhizoMAG_%03d.fa' "$bin_index")"
    new_bin_id="${new_filename%.fa}"

    cp "$bin_file" "${renamed_bins}/${new_filename}"

    printf "%s\t%s\n" \
        "$new_bin_id" \
        "$original_name" \
        >> "$mapping_file"

    bin_index=$((bin_index + 1))
done < <(
    find "$dastool_bins" \
        -maxdepth 1 \
        -type f \
        -name '*.fa' \
        -print \
        | sort
)

bin_count=$((bin_index - 1))

if [[ "$bin_count" -eq 0 ]]; then
    echo "ERROR: No bins were renamed." >&2
    exit 1
fi

echo "[$(date)] Renamed $bin_count bins"
echo "[$(date)] Name mapping: $mapping_file"


# ---------------------------------------------------------------------------
# Run CheckM
# ---------------------------------------------------------------------------

echo "[$(date)] Running CheckM lineage workflow"

checkm lineage_wf \
    -t "$threads" \
    -x fa \
    "$renamed_bins" \
    "$checkm_output"


# ---------------------------------------------------------------------------
# Generate CheckM QA table
# ---------------------------------------------------------------------------

summary_file="${checkm_output}/checkm_summary.txt"

echo "[$(date)] Generating CheckM QA table"

checkm qa \
    "${checkm_output}/lineage.ms" \
    "$checkm_output" \
    -o 2 \
    --tab_table \
    > "$summary_file"


# ---------------------------------------------------------------------------
# Select bins using the thesis thresholds
# ---------------------------------------------------------------------------

minimum_completeness=70
maximum_contamination=10
selected_file="${checkm_output}/high_quality_bins.txt"

awk -F '\t' \
    -v minimum_completeness="$minimum_completeness" \
    -v maximum_contamination="$maximum_contamination" \
    '
    $1 == "Bin Id" {
        found_header = 1
        next
    }

    found_header && NF > 6 {
        completeness = $6 + 0
        contamination = $7 + 0

        if (
            completeness >= minimum_completeness &&
            contamination <= maximum_contamination
        ) {
            print $1 "\t" completeness "\t" contamination
        }
    }
    ' "$summary_file" \
    > "$selected_file"

selected_count=$(wc -l < "$selected_file")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] CheckM evaluation completed"
echo "Bins evaluated:      $bin_count"
echo "Bins selected:       $selected_count"
echo "Selection threshold: completeness >= ${minimum_completeness}%, contamination <= ${maximum_contamination}%"
echo "CheckM summary:      $summary_file"
echo "Selected bins:       $selected_file"
echo "Renamed bins:        $renamed_bins"
echo "Name mapping:        $mapping_file"
echo "Runtime:             ${runtime_minutes} minutes"

if [[ "$selected_count" -gt 0 ]]; then
    echo
    echo "Selected bins:"
    cat "$selected_file"
fi