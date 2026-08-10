#!/usr/bin/env bash
#
# ==============================================================================
# Script: run_memote.sh
#
# Purpose:
#   Generate an independent MEMOTE snapshot report for every CarveMe SBML model
#   in a specified directory.
#
# Usage:
#   bash run_memote.sh MODELS_DIR OUTPUT_DIR
#
# Example:
#   bash run_memote.sh \
#       /path/to/models \
#       /path/to/memote_reports
#
# Optional environment variables:
#   MEMOTE_ENV
#       Name of the Conda environment containing MEMOTE.
#       When set, the script activates this environment before running MEMOTE.
#
#   CONDA_SH
#       Full path to conda.sh. Required only when MEMOTE_ENV is set and
#       conda.sh cannot be found automatically.
#
# Outputs:
#   One HTML snapshot report per input model:
#
#       OUTPUT_DIR/<model_id>.html
#
# Reproducibility:
#   This script preserves the original command:
#
#       memote report snapshot MODEL --filename REPORT.html
#
#   Reports are placed in an explicit output directory rather than whichever
#   directory happens to be active when the loop is launched.
# ==============================================================================

set -euo pipefail


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    cat >&2 <<EOF
Usage:
  bash $0 MODELS_DIR OUTPUT_DIR

Arguments:
  MODELS_DIR  Directory containing SBML models in XML format.
  OUTPUT_DIR  Directory in which MEMOTE HTML reports will be written.
EOF
    exit 1
fi

MODELS_DIR="$1"
OUTPUT_DIR="$2"


# ------------------------------------------------------------------------------
# Validate input and create output directory
# ------------------------------------------------------------------------------

if [[ ! -d "$MODELS_DIR" ]]; then
    echo "ERROR: Models directory does not exist: $MODELS_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"


# ------------------------------------------------------------------------------
# Optionally activate a dedicated MEMOTE Conda environment
# ------------------------------------------------------------------------------

if [[ -n "${MEMOTE_ENV:-}" ]]; then
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
        echo "Set CONDA_SH to the full path of conda.sh." >&2
        exit 1
    fi

    # shellcheck source=/dev/null
    source "$CONDA_SH"
    conda activate "$MEMOTE_ENV"
fi


# ------------------------------------------------------------------------------
# Validate MEMOTE availability
# ------------------------------------------------------------------------------

if ! command -v memote >/dev/null 2>&1; then
    echo "ERROR: The 'memote' command was not found." >&2

    if [[ -n "${MEMOTE_ENV:-}" ]]; then
        echo "Activated Conda environment: $MEMOTE_ENV" >&2
    else
        echo "Activate the MEMOTE environment or set MEMOTE_ENV." >&2
    fi

    exit 1
fi

echo "=== MEMOTE snapshot reports ==="
echo "MEMOTE executable: $(command -v memote)"

if memote --version >/dev/null 2>&1; then
    echo "MEMOTE version:    $(memote --version 2>&1 | head -n 1)"
fi

echo "Models directory:  $MODELS_DIR"
echo "Output directory:  $OUTPUT_DIR"
echo


# ------------------------------------------------------------------------------
# Select models
#
# This intentionally preserves the original non-recursive *.xml search.
# ------------------------------------------------------------------------------

shopt -s nullglob
models=("${MODELS_DIR}"/*.xml)
shopt -u nullglob

if [[ ${#models[@]} -eq 0 ]]; then
    echo "ERROR: No XML models found under: ${MODELS_DIR}/*.xml" >&2
    exit 1
fi

echo "Models found: ${#models[@]}"
echo


# ------------------------------------------------------------------------------
# Generate one MEMOTE snapshot per model
# ------------------------------------------------------------------------------

reports_completed=0

for model in "${models[@]}"; do
    base="$(basename "$model" .xml)"
    report="${OUTPUT_DIR}/${base}.html"

    echo "[$(date)] Running MEMOTE"
    echo "Model:  $model"
    echo "Report: $report"

    memote report snapshot "$model" \
        --filename "$report"

    if [[ ! -s "$report" ]]; then
        echo "ERROR: MEMOTE did not produce a non-empty report: $report" >&2
        exit 1
    fi

    ((reports_completed += 1))

    echo "[$(date)] Completed: $base"
    echo
done


# ------------------------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------------------------

echo "=== MEMOTE completed successfully ==="
echo "Reports completed: $reports_completed"
echo "Reports directory: $OUTPUT_DIR/"
