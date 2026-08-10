#!/usr/bin/env bash
#
# ==============================================================================
# Script: runCarveMe.sh
#
# Purpose:
#   Reconstruct genome-scale metabolic models from Prokka protein FASTA files
#   using CarveMe. A Slurm array assigns one treatment-specific medium to each
#   task:
#
#       Array task 0: CK
#       Array task 1: NP
#       Array task 2: NPM
#
#   Only MAGs whose protein FASTA filenames begin with the selected medium name
#   are processed. Each model is gap-filled using its corresponding medium.
#
# Usage:
#   sbatch runCarveMe.sh \
#       PROKKA_DIR \
#       MEDIA_DATABASE \
#       OUTPUT_BASE_DIR
#
# Example:
#   sbatch runCarveMe.sh \
#       /path/to/prokka_HQ_MAGs \
#       /path/to/media_groups.tsv \
#       /path/to/carveme_models
#
# Optional environment variables:
#   CARVEME_ENV   Conda environment containing CarveMe.
#                 Default: carveme_env
#
#   CONDA_SH      Path to conda.sh.
#                 If unset, the script attempts to locate it automatically.
#
# Expected input structure:
#   PROKKA_DIR/
#   ├── MAG_directory_1/
#   │   └── CK_example.faa
#   ├── MAG_directory_2/
#   │   └── NP_example.faa
#   └── MAG_directory_3/
#       └── NPM_example.faa
#
# Outputs:
#   OUTPUT_BASE_DIR/
#   ├── CK/*.xml
#   ├── NP/*.xml
#   └── NPM/*.xml
#
# Reproducibility-critical CarveMe settings:
#   --fbc2
#   --solver scip
#   --gapfill <CK|NP|NPM>
#   --mediadb <MEDIA_DATABASE>
#
# Notes:
#   - The media database must be a tab-separated file whose first column
#     contains the medium identifiers CK, NP, and NPM.
#   - Existing XML files with identical names may be overwritten by CarveMe.
#   - Slurm log files are written to the directory from which sbatch is run,
#     unless --output and --error are supplied during submission.
# ==============================================================================

#SBATCH --job-name=carveme
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=carveme_%x_%A_%a.log
#SBATCH --error=carveme_%x_%A_%a.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    cat >&2 <<EOF
Usage:
  sbatch $0 PROKKA_DIR MEDIA_DATABASE OUTPUT_BASE_DIR

Arguments:
  PROKKA_DIR       Directory containing per-MAG Prokka output directories.
  MEDIA_DATABASE  CarveMe media database in TSV format.
  OUTPUT_BASE_DIR Base directory for medium-specific XML model directories.
EOF
    exit 1
fi

PROKKA_DIR="$1"
MEDIA_DB="$2"
OUTPUT_BASE_DIR="$3"

CARVEME_ENV="${CARVEME_ENV:-carveme_env}"


# ------------------------------------------------------------------------------
# Select the medium based on the Slurm array index
# ------------------------------------------------------------------------------

MEDIA_LIST=("CK" "NP" "NPM")

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set." >&2
    echo "Submit this script as a Slurm array job using sbatch." >&2
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 0 ||
      SLURM_ARRAY_TASK_ID >= ${#MEDIA_LIST[@]} )); then
    echo "ERROR: Invalid Slurm array index: ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

MED="${MEDIA_LIST[$SLURM_ARRAY_TASK_ID]}"
OUTDIR="${OUTPUT_BASE_DIR}/${MED}"

echo "[$(date)] Starting CarveMe reconstruction"
echo "Medium: ${MED}"
echo "Slurm job ID: ${SLURM_JOB_ID:-not_available}"
echo "Slurm array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Protein input directory: ${PROKKA_DIR}"
echo "Media database: ${MEDIA_DB}"
echo "Output directory: ${OUTDIR}"

start_time=$(date +%s)


# ------------------------------------------------------------------------------
# Locate and activate the Conda environment
# ------------------------------------------------------------------------------

if [[ -n "${CONDA_SH:-}" ]]; then
    if [[ ! -f "$CONDA_SH" ]]; then
        echo "ERROR: CONDA_SH does not exist: $CONDA_SH" >&2
        exit 1
    fi
elif [[ -n "${CONDA_EXE:-}" ]]; then
    conda_base="$(dirname "$(dirname "$CONDA_EXE")")"
    CONDA_SH="${conda_base}/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
    CONDA_SH="${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
    CONDA_SH="${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: Could not locate conda.sh." >&2
    echo "Set CONDA_SH to the full path of your conda.sh file." >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$CARVEME_ENV"


# ------------------------------------------------------------------------------
# Validate software availability
# ------------------------------------------------------------------------------

if ! command -v carve >/dev/null 2>&1; then
    echo "ERROR: The 'carve' command was not found after activating" >&2
    echo "       Conda environment '${CARVEME_ENV}'." >&2
    exit 1
fi

echo "Conda environment: ${CONDA_DEFAULT_ENV:-unknown}"
echo "CarveMe executable: $(command -v carve)"

if carve --version >/dev/null 2>&1; then
    echo "CarveMe version: $(carve --version 2>&1 | head -n 1)"
fi


# ------------------------------------------------------------------------------
# Validate input paths
# ------------------------------------------------------------------------------

if [[ ! -d "$PROKKA_DIR" ]]; then
    echo "ERROR: Prokka directory does not exist: $PROKKA_DIR" >&2
    exit 1
fi

if [[ ! -r "$PROKKA_DIR" ]]; then
    echo "ERROR: Prokka directory is not readable: $PROKKA_DIR" >&2
    exit 1
fi

if [[ ! -s "$MEDIA_DB" ]]; then
    echo "ERROR: Media database was not found or is empty: $MEDIA_DB" >&2
    exit 1
fi

if [[ ! -r "$MEDIA_DB" ]]; then
    echo "ERROR: Media database is not readable: $MEDIA_DB" >&2
    exit 1
fi

if ! tail -n +2 "$MEDIA_DB" |
     cut -f1 |
     grep -Fqx "$MED"; then
    echo "ERROR: Medium '${MED}' was not found in: $MEDIA_DB" >&2
    exit 1
fi

mkdir -p "$OUTDIR"


# ------------------------------------------------------------------------------
# Select only MAGs belonging to the current medium
#
# This preserves the original input-selection rule:
#   PROKKA_DIR/*/<MEDIUM>_*.faa
# ------------------------------------------------------------------------------

shopt -s nullglob
faa_files=("${PROKKA_DIR}"/*/"${MED}"_*.faa)
shopt -u nullglob

if [[ ${#faa_files[@]} -eq 0 ]]; then
    echo "No ${MED} MAG protein FASTA files were found."
    echo "Expected pattern: ${PROKKA_DIR}/*/${MED}_*.faa"
    exit 0
fi

echo "Number of ${MED} MAGs found: ${#faa_files[@]}"


# ------------------------------------------------------------------------------
# Run CarveMe
#
# Each MAG is reconstructed and gap-filled only with its corresponding medium.
# The CarveMe options below are intentionally unchanged from the original run.
# ------------------------------------------------------------------------------

models_completed=0

for faa in "${faa_files[@]}"; do
    id="$(basename "$faa" .faa)"
    out="${OUTDIR}/${id}.xml"

    if [[ ! -s "$faa" ]]; then
        echo "ERROR: Protein FASTA file is empty: $faa" >&2
        exit 1
    fi

    echo
    echo "[$(date)] Building GEM"
    echo "MAG: ${id}"
    echo "Input: ${faa}"
    echo "Medium: ${MED}"
    echo "Output: ${out}"

    carve "$faa" \
        -o "$out" \
        --fbc2 \
        --solver scip \
        --gapfill "$MED" \
        --mediadb "$MEDIA_DB"

    if [[ ! -s "$out" ]]; then
        echo "ERROR: CarveMe did not produce a non-empty model: $out" >&2
        exit 1
    fi

    ((models_completed += 1))
done


# ------------------------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------------------------

end_time=$(date +%s)
elapsed_seconds=$((end_time - start_time))
elapsed_minutes=$((elapsed_seconds / 60))

echo
echo "[$(date)] CarveMe reconstruction completed successfully"
echo "Medium: ${MED}"
echo "Models completed: ${models_completed}"
echo "Elapsed time: ${elapsed_minutes} minutes"
echo "Output directory: ${OUTDIR}/"
