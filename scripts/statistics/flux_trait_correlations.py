#!/usr/bin/env python3
"""
Correlate FVA-selected MICOM exchange fluxes with plant traits.

This script evaluates associations between community-level exchange fluxes
and measured plant traits using Spearman rank correlation.

The analysis is restricted to the exchange reactions retained after the
FVA-based variability filtering step.

For each selected exchange reaction, correlations are calculated separately
for:

    - Maize yield
    - Total acid

Benjamini-Hochberg FDR correction is applied separately within each trait.

Variables are skipped when:

    - all flux values are missing,
    - the flux is constant across samples,
    - fewer than three paired observations are available,
    - the plant trait is constant.

Usage
-----
python flux_trait_correlations.py \
    --flux-file community_exchange_fluxes_selected.tsv \
    --trait-file plant_trait.tsv \
    --output-dir results

Outputs
-------
flux_trait_correlations.tsv
    Spearman correlation coefficients, raw p-values, and FDR-adjusted
    p-values.

flux_trait_merged.tsv
    Merged flux and plant-trait table used for the analysis.

flux_trait_correlations_skipped.tsv
    Flux-trait combinations excluded from testing and the reason for
    exclusion.

Input requirements
------------------
Flux table:
    Must contain a "sample" column and exchange reaction columns beginning
    with EX_.

Trait table:
    Must contain:
        SRR
        Total acid
        Maize yield

Dependencies
------------
- Python >= 3.8
- pandas
- scipy
- statsmodels
"""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path

import pandas as pd
import scipy
import statsmodels
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Analysis constants
# ------------------------------------------------------------------------------

TRAITS_TO_TEST = (
    "Maize yield",
    "Total acid",
)

FLUX_PREFIX = "EX_"

MIN_OBSERVATIONS = 3

EXCLUDE_COLUMNS = (
    "condition",
    "sample",
    "sample_dir",
    "SRR",
    "Sample",
    "Total acid",
    "Maize yield",
)


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Spearman correlations between FVA-selected MICOM "
            "exchange fluxes and plant traits."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--flux-file",
        required=True,
        type=Path,
        help=(
            "FVA-selected community exchange-flux TSV containing a "
            "'sample' column."
        ),
    )

    parser.add_argument(
        "--trait-file",
        required=True,
        type=Path,
        help=(
            "Plant-trait TSV containing SRR, Maize yield, and Total acid."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which analysis outputs will be written.",
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------

def validate_file(
    path: Path,
    description: str,
) -> Path:
    """Resolve and validate an input file."""
    path = path.expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(
            f"{description} does not exist or is not a file: {path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"{description} is empty: {path}"
        )

    return path


def load_data(
    flux_file: Path,
    trait_file: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read and validate flux and trait input tables."""
    flux = pd.read_csv(
        flux_file,
        sep="\t",
    )

    traits = pd.read_csv(
        trait_file,
        sep="\t",
    )

    flux.columns = [
        column.strip()
        for column in flux.columns
    ]

    traits.columns = [
        column.strip()
        for column in traits.columns
    ]

    required_flux_columns = {
        "sample",
    }

    missing_flux = (
        required_flux_columns
        .difference(flux.columns)
    )

    if missing_flux:
        raise ValueError(
            "Flux table is missing required columns: "
            f"{sorted(missing_flux)}"
        )

    required_trait_columns = {
        "SRR",
        "Total acid",
        "Maize yield",
    }

    missing_traits = (
        required_trait_columns
        .difference(traits.columns)
    )

    if missing_traits:
        raise ValueError(
            "Trait table is missing required columns: "
            f"{sorted(missing_traits)}"
        )

    traits["Total acid"] = pd.to_numeric(
        traits["Total acid"],
        errors="coerce",
    )

    traits["Maize yield"] = pd.to_numeric(
        traits["Maize yield"],
        errors="coerce",
    )

    return (
        flux,
        traits,
    )


# ------------------------------------------------------------------------------
# Flux selection
# ------------------------------------------------------------------------------

def select_flux_columns(
    dataframe: pd.DataFrame,
) -> list[str]:
    """
    Select numeric MICOM exchange-flux columns.

    Only columns beginning with EX_ are retained.
    """
    return [
        column
        for column in dataframe.columns
        if (
            column not in EXCLUDE_COLUMNS
            and column.startswith(FLUX_PREFIX)
            and pd.api.types.is_numeric_dtype(
                dataframe[column]
            )
        )
    ]


# ------------------------------------------------------------------------------
# Correlation analysis
# ------------------------------------------------------------------------------

def run_correlations(
    dataframe: pd.DataFrame,
    flux_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate flux-trait Spearman correlations.

    FDR correction is applied separately within each trait.
    """
    results = []
    skipped = []

    for trait in TRAITS_TO_TEST:
        y_values = pd.to_numeric(
            dataframe[trait],
            errors="coerce",
        )

        for flux_column in flux_columns:
            x_values = pd.to_numeric(
                dataframe[flux_column],
                errors="coerce",
            )

            if x_values.notna().sum() == 0:
                skipped.append(
                    {
                        "flux": flux_column,
                        "trait": trait,
                        "reason": "all NaN",
                    }
                )
                continue

            if (
                x_values
                .nunique(
                    dropna=True
                )
                <= 1
            ):
                skipped.append(
                    {
                        "flux": flux_column,
                        "trait": trait,
                        "reason": "constant across samples",
                    }
                )
                continue

            paired = pd.DataFrame(
                {
                    "x": x_values,
                    "y": y_values,
                }
            ).dropna()

            if len(paired) < MIN_OBSERVATIONS:
                skipped.append(
                    {
                        "flux": flux_column,
                        "trait": trait,
                        "reason": "too few observations",
                    }
                )
                continue

            if (
                paired["y"]
                .nunique(
                    dropna=True
                )
                <= 1
            ):
                skipped.append(
                    {
                        "flux": flux_column,
                        "trait": trait,
                        "reason": "trait constant",
                    }
                )
                continue

            rho, p_value = spearmanr(
                paired["x"],
                paired["y"],
            )

            results.append(
                {
                    "flux": flux_column,
                    "trait": trait,
                    "n": len(paired),
                    "rho": rho,
                    "p_value": p_value,
                }
            )

    results_df = pd.DataFrame(
        results,
        columns=[
            "flux",
            "trait",
            "n",
            "rho",
            "p_value",
        ],
    )

    skipped_df = pd.DataFrame(
        skipped,
        columns=[
            "flux",
            "trait",
            "reason",
        ],
    )

    if results_df.empty:
        raise SystemExit(
            "No valid correlations could be computed."
        )

    # --------------------------------------------------------------------------
    # Benjamini-Hochberg FDR correction
    #
    # Correction is intentionally performed separately within each trait.
    # --------------------------------------------------------------------------

    results_df["p_adj"] = pd.NA

    for trait in TRAITS_TO_TEST:
        mask = (
            results_df["trait"]
            == trait
        )

        if mask.sum() > 0:
            results_df.loc[
                mask,
                "p_adj",
            ] = multipletests(
                results_df.loc[
                    mask,
                    "p_value",
                ],
                method="fdr_bh",
            )[1]

    results_df["p_adj"] = pd.to_numeric(
        results_df["p_adj"],
        errors="coerce",
    )

    results_df = results_df.sort_values(
        [
            "trait",
            "p_adj",
            "p_value",
            "rho",
        ],
        ascending=[
            True,
            True,
            True,
            False,
        ],
    )

    return (
        results_df,
        skipped_df,
    )


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Run the flux-trait correlation analysis."""
    args = parse_args()

    flux_file = validate_file(
        args.flux_file,
        "Flux file",
    )

    trait_file = validate_file(
        args.trait_file,
        "Trait file",
    )

    output_dir = (
        args.output_dir
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    results_file = (
        output_dir
        / "flux_trait_correlations.tsv"
    )

    merged_file = (
        output_dir
        / "flux_trait_merged.tsv"
    )

    skipped_file = (
        output_dir
        / "flux_trait_correlations_skipped.tsv"
    )

    print("=== Flux-trait correlation analysis ===")
    print(
        f"Python version:      "
        f"{platform.python_version()}"
    )
    print(
        f"pandas version:      "
        f"{pd.__version__}"
    )
    print(
        f"SciPy version:       "
        f"{scipy.__version__}"
    )
    print(
        f"statsmodels version: "
        f"{statsmodels.__version__}"
    )
    print(
        f"Flux file:           "
        f"{flux_file}"
    )
    print(
        f"Trait file:          "
        f"{trait_file}"
    )
    print(
        f"Output directory:    "
        f"{output_dir}"
    )
    print()

    # --------------------------------------------------------------------------
    # Load data
    # --------------------------------------------------------------------------

    flux, traits = load_data(
        flux_file,
        trait_file,
    )

    # --------------------------------------------------------------------------
    # Merge by sample identifier
    #
    # MICOM table: sample
    # Trait table: SRR
    # --------------------------------------------------------------------------

    merged = flux.merge(
        traits,
        left_on="sample",
        right_on="SRR",
        how="inner",
    )

    if merged.empty:
        raise ValueError(
            "Flux and trait tables produced no matching samples "
            "when merging 'sample' with 'SRR'."
        )

    merged.to_csv(
        merged_file,
        sep="\t",
        index=False,
    )

    print(
        "Merged table shape:"
    )
    print(
        merged.shape
    )
    print()

    display_columns = [
        column
        for column in (
            "sample",
            "Sample",
            "condition",
            "Total acid",
            "Maize yield",
        )
        if column in merged.columns
    ]

    print(
        "Merged samples:"
    )

    print(
        merged[
            display_columns
        ].to_string(
            index=False
        )
    )

    # --------------------------------------------------------------------------
    # Select eligible exchange fluxes
    # --------------------------------------------------------------------------

    flux_columns = select_flux_columns(
        merged
    )

    if not flux_columns:
        raise ValueError(
            "No numeric EX_ exchange-flux columns were found."
        )

    print()
    print(
        f"Testing {len(flux_columns)} exchange flux "
        "variables per trait."
    )
    print()

    # --------------------------------------------------------------------------
    # Correlations
    # --------------------------------------------------------------------------

    (
        results,
        skipped,
    ) = run_correlations(
        merged,
        flux_columns,
    )

    # --------------------------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------------------------

    results.to_csv(
        results_file,
        sep="\t",
        index=False,
    )

    skipped.to_csv(
        skipped_file,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Console summary
    # --------------------------------------------------------------------------

    print(
        "Top correlations:"
    )
    print()

    print(
        results
        .head(20)
        .to_string(
            index=False
        )
    )

    print()
    print(
        f"Results saved to:\n{results_file}"
    )

    print(
        f"Merged table saved to:\n{merged_file}"
    )

    print(
        f"Skipped variables saved to:\n{skipped_file}"
    )

    print()
    print(
        f"Tested correlations: {len(results)}"
    )

    print(
        f"Skipped entries:     {len(skipped)}"
    )

    print()
    print(
        "FDR correction was applied separately "
        "within each trait."
    )


if __name__ == "__main__":
    main()