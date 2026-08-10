#!/usr/bin/env python3
"""
Select community exchange reactions for downstream differential analysis
using MICOM FVA results.

This script summarizes reaction-wise flux variability across all sample-level
FVA result files and selects candidate reactions using three criteria:

1. Median absolute FVA range greater than a minimum threshold.
2. Median relative FVA range below a maximum threshold.
3. Reaction classified as variable in at least a specified fraction of samples.

The selected reactions are then extracted from the community exchange-flux
table for downstream statistical analyses.

Selection criteria
------------------
By default:

    median_range > 1e-6

    median_relative_range < 0.95

    fraction_variable >= 0.5

where:

    absolute flux scale =
        max(|minimum|, |maximum|)

    relative range =
        range / (absolute flux scale + 1e-9)

Usage
-----
python select_fva_reactions.py \
    --fva-dir /path/to/fva \
    --flux-file /path/to/community_exchange_fluxes.tsv \
    --output-dir /path/to/output

Outputs
-------
community_exchange_fluxes_selected.tsv
    Community exchange-flux table containing only FVA-selected reactions
    plus available metadata columns.

fva_selected_reactions_used.tsv
    Summary of selected reactions and their FVA statistics.

reaction_range_summary_rebuilt.tsv
    Complete reaction-level FVA summary.

Dependencies
------------
- Python >= 3.8
- pandas
"""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path

import pandas as pd


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Default FVA-selection thresholds
#
# These values reproduce the original downstream filtering analysis.
# ------------------------------------------------------------------------------

DEFAULT_ABS_RANGE_MIN = 1e-6
DEFAULT_REL_RANGE_MAX = 0.95
DEFAULT_MIN_FRAC_VARIABLE = 0.5

RELATIVE_RANGE_EPSILON = 1e-9

METADATA_COLUMNS = (
    "sample",
    "condition",
    "sample_dir",
)


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Select MICOM community exchange reactions for downstream "
            "analysis using sample-level FVA results."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--fva-dir",
        required=True,
        type=Path,
        help="Directory containing per-sample FVA CSV files.",
    )

    parser.add_argument(
        "--flux-file",
        required=True,
        type=Path,
        help=(
            "Community exchange-flux TSV from which selected reaction "
            "columns will be extracted."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which selected outputs will be written.",
    )

    parser.add_argument(
        "--abs-range-min",
        type=float,
        default=DEFAULT_ABS_RANGE_MIN,
        help=(
            "Minimum median absolute FVA range required for a reaction "
            "to pass."
        ),
    )

    parser.add_argument(
        "--rel-range-max",
        type=float,
        default=DEFAULT_REL_RANGE_MAX,
        help=(
            "Maximum median relative FVA range allowed for a reaction "
            "to pass."
        ),
    )

    parser.add_argument(
        "--min-frac-variable",
        type=float,
        default=DEFAULT_MIN_FRAC_VARIABLE,
        help=(
            "Minimum fraction of samples in which the reaction must "
            "be classified as variable."
        ),
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------

def validate_arguments(
    args: argparse.Namespace,
) -> None:
    """Validate paths and filtering thresholds."""
    args.fva_dir = (
        args.fva_dir
        .expanduser()
        .resolve()
    )

    args.flux_file = (
        args.flux_file
        .expanduser()
        .resolve()
    )

    args.output_dir = (
        args.output_dir
        .expanduser()
        .resolve()
    )

    if not args.fva_dir.is_dir():
        raise FileNotFoundError(
            "FVA directory does not exist or is not a directory: "
            f"{args.fva_dir}"
        )

    if not args.flux_file.is_file():
        raise FileNotFoundError(
            "Community flux file does not exist: "
            f"{args.flux_file}"
        )

    if args.flux_file.stat().st_size == 0:
        raise ValueError(
            f"Community flux file is empty: {args.flux_file}"
        )

    if args.abs_range_min < 0:
        raise ValueError(
            "--abs-range-min cannot be negative."
        )

    if args.rel_range_max < 0:
        raise ValueError(
            "--rel-range-max cannot be negative."
        )

    if not 0 <= args.min_frac_variable <= 1:
        raise ValueError(
            "--min-frac-variable must be between 0 and 1."
        )


# ------------------------------------------------------------------------------
# FVA loading
# ------------------------------------------------------------------------------

def load_fva_results(
    fva_dir: Path,
) -> pd.DataFrame:
    """
    Load and combine all per-sample FVA CSV files.

    The sample name is taken from the CSV filename and the condition is
    extracted as the text preceding the first underscore.
    """
    files = sorted(
        fva_dir.glob("*.csv")
    )

    if not files:
        raise FileNotFoundError(
            f"No FVA CSV files found in {fva_dir}"
        )

    dataframes = []

    for file_path in files:
        dataframe = pd.read_csv(
            file_path,
            index_col=0,
        )

        required_columns = {
            "minimum",
            "maximum",
            "range",
            "variable",
        }

        missing = (
            required_columns
            .difference(dataframe.columns)
        )

        if missing:
            raise ValueError(
                f"FVA file {file_path} is missing columns: "
                f"{sorted(missing)}"
            )

        sample = file_path.stem
        condition = sample.split("_")[0]

        dataframe["sample"] = sample
        dataframe["condition"] = condition

        dataframe["abs_flux_scale"] = (
            dataframe[
                ["minimum", "maximum"]
            ]
            .abs()
            .max(axis=1)
        )

        dataframe["relative_range"] = (
            dataframe["range"]
            / (
                dataframe["abs_flux_scale"]
                + RELATIVE_RANGE_EPSILON
            )
        )

        dataframes.append(
            dataframe
        )

    fva_all = pd.concat(
        dataframes
    )

    fva_all.index.name = (
        "reaction"
    )

    return fva_all


# ------------------------------------------------------------------------------
# Reaction-level FVA summary
# ------------------------------------------------------------------------------

def summarize_reactions(
    fva_all: pd.DataFrame,
    abs_range_min: float,
    rel_range_max: float,
    min_frac_variable: float,
) -> pd.DataFrame:
    """Calculate reaction-level FVA statistics and selection flags."""
    reaction_stats = (
        fva_all
        .groupby(level=0)
        .agg(
            median_range=(
                "range",
                "median",
            ),
            max_range=(
                "range",
                "max",
            ),
            mean_range=(
                "range",
                "mean",
            ),
            median_relative_range=(
                "relative_range",
                "median",
            ),
            max_relative_range=(
                "relative_range",
                "max",
            ),
            n_variable=(
                "variable",
                "sum",
            ),
            n_samples=(
                "range",
                "count",
            ),
            median_minimum=(
                "minimum",
                "median",
            ),
            median_maximum=(
                "maximum",
                "median",
            ),
        )
    )

    reaction_stats["frac_variable"] = (
        reaction_stats["n_variable"]
        / reaction_stats["n_samples"]
    )

    reaction_stats["median_center"] = (
        reaction_stats["median_minimum"]
        + reaction_stats["median_maximum"]
    ) / 2

    reaction_stats[
        "passes_abs_range_filter"
    ] = (
        reaction_stats["median_range"]
        > abs_range_min
    )

    reaction_stats[
        "passes_rel_range_filter"
    ] = (
        reaction_stats[
            "median_relative_range"
        ]
        < rel_range_max
    )

    reaction_stats[
        "passes_variability_filter"
    ] = (
        reaction_stats["frac_variable"]
        >= min_frac_variable
    )

    reaction_stats[
        "candidate_for_downstream"
    ] = (
        reaction_stats[
            "passes_abs_range_filter"
        ]
        & reaction_stats[
            "passes_rel_range_filter"
        ]
        & reaction_stats[
            "passes_variability_filter"
        ]
    )

    return reaction_stats


# ------------------------------------------------------------------------------
# Candidate selection
# ------------------------------------------------------------------------------

def select_candidates(
    reaction_stats: pd.DataFrame,
) -> list[str]:
    """
    Return reactions passing all filters.

    Candidate ordering reproduces the original script:
    first by median relative range and then by median absolute range.
    """
    return (
        reaction_stats[
            reaction_stats[
                "candidate_for_downstream"
            ]
        ]
        .sort_values(
            [
                "median_relative_range",
                "median_range",
            ]
        )
        .index
        .astype(str)
        .tolist()
    )


# ------------------------------------------------------------------------------
# Flux-table filtering
# ------------------------------------------------------------------------------

def select_flux_columns(
    flux_file: Path,
    candidates: list[str],
    reaction_stats: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract selected reaction columns from the community flux table.
    """
    flux_df = pd.read_csv(
        flux_file,
        sep="\t",
    )

    metadata_columns = [
        column
        for column in METADATA_COLUMNS
        if column in flux_df.columns
    ]

    selected_columns = []
    report_rows = []

    for reaction in candidates:
        selected_column = (
            reaction
            if reaction in flux_df.columns
            else None
        )

        if selected_column is not None:
            selected_columns.append(
                selected_column
            )

        report_rows.append(
            {
                "reaction": reaction,
                "median_range": reaction_stats.loc[
                    reaction,
                    "median_range",
                ],
                "median_relative_range": reaction_stats.loc[
                    reaction,
                    "median_relative_range",
                ],
                "frac_variable": reaction_stats.loc[
                    reaction,
                    "frac_variable",
                ],
                "median_center": reaction_stats.loc[
                    reaction,
                    "median_center",
                ],
                "selected_column": (
                    selected_column
                    if selected_column is not None
                    else "NOT_FOUND"
                ),
            }
        )

    selected = flux_df[
        metadata_columns
        + selected_columns
    ].copy()

    report = pd.DataFrame(
        report_rows
    )

    return (
        selected,
        report,
    )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main() -> None:
    """Run FVA-based reaction selection."""
    args = parse_args()

    validate_arguments(
        args
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    selected_flux_output = (
        args.output_dir
        / "community_exchange_fluxes_selected.tsv"
    )

    selected_report_output = (
        args.output_dir
        / "fva_selected_reactions_used.tsv"
    )

    reaction_summary_output = (
        args.fva_dir
        / "reaction_range_summary_rebuilt.tsv"
    )

    print("=== FVA-based reaction selection ===")
    print(f"Python version:       {platform.python_version()}")
    print(f"pandas version:       {pd.__version__}")
    print(f"FVA directory:        {args.fva_dir}")
    print(f"Flux file:            {args.flux_file}")
    print(f"Output directory:     {args.output_dir}")
    print(f"Absolute range min:   {args.abs_range_min}")
    print(f"Relative range max:   {args.rel_range_max}")
    print(
        "Minimum variable frac: "
        f"{args.min_frac_variable}"
    )
    print()

    fva_all = load_fva_results(
        args.fva_dir
    )

    reaction_stats = summarize_reactions(
        fva_all=fva_all,
        abs_range_min=args.abs_range_min,
        rel_range_max=args.rel_range_max,
        min_frac_variable=args.min_frac_variable,
    )

    candidates = select_candidates(
        reaction_stats
    )

    selected_fluxes, report = (
        select_flux_columns(
            flux_file=args.flux_file,
            candidates=candidates,
            reaction_stats=reaction_stats,
        )
    )

    selected_fluxes.to_csv(
        selected_flux_output,
        sep="\t",
        index=False,
    )

    report.to_csv(
        selected_report_output,
        sep="\t",
        index=False,
    )

    reaction_stats.to_csv(
        reaction_summary_output,
        sep="\t",
    )

    print("FVA candidates selected:")

    for reaction in candidates:
        print(
            f"  {reaction}"
        )

    print()
    print("Selected columns:")

    for column in selected_fluxes.columns:
        print(
            f"  {column}"
        )

    print()
    print("Saved:")
    print(
        selected_flux_output
    )
    print(
        selected_report_output
    )
    print(
        reaction_summary_output
    )


if __name__ == "__main__":
    main()