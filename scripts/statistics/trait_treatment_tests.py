#!/usr/bin/env python3
"""
Test plant-trait differences among fertilization treatments.

This script evaluates treatment-associated differences in the plant traits
"Total acid" and "Maize yield" across CK, NP, and NPM samples.

For each trait, the script performs:

1. A Kruskal-Wallis test across CK, NP, and NPM.
2. Pairwise two-sided Mann-Whitney U tests for:
       CK vs NP
       CK vs NPM
       NP vs NPM
3. Holm correction across the three pairwise comparisons separately
   for each trait.

Treatment labels are extracted from the "Sample" column using the pattern:

    NPM | NP | CK

Usage
-----
python trait_treatment_tests.py \
    --input plant_trait.tsv \
    --output-dir results

Outputs
-------
trait_treatment_kruskalHolm.tsv
    Kruskal-Wallis statistics and treatment-level descriptive statistics.

trait_pairwise_testsHolm.tsv
    Pairwise Mann-Whitney U tests with Holm-adjusted p-values.

Input requirements
------------------
The input TSV must contain:

    Sample
    Total acid
    Maize yield

An "SRR" column is displayed in the console when present.

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
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Analysis constants
# ------------------------------------------------------------------------------

NUMERIC_TRAITS = (
    "Total acid",
    "Maize yield",
)

TREATMENTS = (
    "CK",
    "NP",
    "NPM",
)

PAIRWISE_COMPARISONS = (
    ("CK", "NP"),
    ("CK", "NPM"),
    ("NP", "NPM"),
)

CONDITION_PATTERN = r"(NPM|NP|CK)"


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Test plant-trait differences among CK, NP, and NPM treatments "
            "using Kruskal-Wallis and pairwise Mann-Whitney U tests."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input plant-trait TSV file.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which result tables will be written.",
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Input loading
# ------------------------------------------------------------------------------

def load_traits(input_file: Path) -> pd.DataFrame:
    """Read and validate the plant-trait table."""
    if not input_file.is_file():
        raise FileNotFoundError(
            f"Input file does not exist: {input_file}"
        )

    if input_file.stat().st_size == 0:
        raise ValueError(
            f"Input file is empty: {input_file}"
        )

    traits = pd.read_csv(
        input_file,
        sep="\t",
    )

    if "Sample" not in traits.columns:
        raise ValueError(
            "Input file must contain a 'Sample' column."
        )

    # Preserve the original treatment-extraction rule.
    traits["condition"] = (
        traits["Sample"]
        .astype(str)
        .str.extract(
            CONDITION_PATTERN,
            expand=False,
        )
    )

    for trait in NUMERIC_TRAITS:
        if trait not in traits.columns:
            raise ValueError(
                f"Input file must contain '{trait}' column."
            )

        traits[trait] = pd.to_numeric(
            traits[trait],
            errors="coerce",
        )

    # Preserve original behavior: only rows lacking treatment assignment
    # are removed at this stage.
    traits = (
        traits
        .dropna(subset=["condition"])
        .copy()
    )

    return traits


# ------------------------------------------------------------------------------
# Statistical analysis
# ------------------------------------------------------------------------------

def run_tests(
    traits: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run Kruskal-Wallis and pairwise Mann-Whitney U tests."""
    kruskal_results = []
    pairwise_results = []

    for trait in NUMERIC_TRAITS:
        sub = (
            traits[
                ["condition", trait]
            ]
            .dropna()
            .copy()
        )

        ck = sub.loc[
            sub["condition"] == "CK",
            trait,
        ]

        np_group = sub.loc[
            sub["condition"] == "NP",
            trait,
        ]

        npm = sub.loc[
            sub["condition"] == "NPM",
            trait,
        ]

        if (
            len(ck) == 0
            or len(np_group) == 0
            or len(npm) == 0
        ):
            raise ValueError(
                "Missing one or more treatment groups "
                f"for trait '{trait}'."
            )

        statistic, p_value = kruskal(
            ck,
            np_group,
            npm,
        )

        kruskal_results.append(
            {
                "trait": trait,
                "n_total": len(sub),
                "n_CK": len(ck),
                "n_NP": len(np_group),
                "n_NPM": len(npm),
                "CK_mean": ck.mean(),
                "NP_mean": np_group.mean(),
                "NPM_mean": npm.mean(),
                "CK_median": ck.median(),
                "NP_median": np_group.median(),
                "NPM_median": npm.median(),
                "H_statistic": statistic,
                "p_value": p_value,
            }
        )

        print(trait)
        print(
            f"  CK  mean = {ck.mean():.3f} "
            f"| median = {ck.median():.3f}"
        )
        print(
            f"  NP  mean = {np_group.mean():.3f} "
            f"| median = {np_group.median():.3f}"
        )
        print(
            f"  NPM mean = {npm.mean():.3f} "
            f"| median = {npm.median():.3f}"
        )
        print(
            "  Kruskal-Wallis "
            f"H = {statistic:.4f}, "
            f"p = {p_value:.6f}"
        )

        trait_pairwise = []

        for group1, group2 in PAIRWISE_COMPARISONS:
            x = sub.loc[
                sub["condition"] == group1,
                trait,
            ]

            y = sub.loc[
                sub["condition"] == group2,
                trait,
            ]

            u_statistic, pairwise_p = (
                mannwhitneyu(
                    x,
                    y,
                    alternative="two-sided",
                )
            )

            trait_pairwise.append(
                {
                    "trait": trait,
                    "group1": group1,
                    "group2": group2,
                    "n_group1": len(x),
                    "n_group2": len(y),
                    "group1_mean": x.mean(),
                    "group2_mean": y.mean(),
                    "group1_median": x.median(),
                    "group2_median": y.median(),
                    "U_statistic": u_statistic,
                    "p_value": pairwise_p,
                }
            )

        # Holm correction is intentionally applied separately within each
        # trait across its three pairwise treatment comparisons.
        p_values = [
            result["p_value"]
            for result in trait_pairwise
        ]

        adjusted_p_values = multipletests(
            p_values,
            method="holm",
        )[1]

        for result, adjusted_p in zip(
            trait_pairwise,
            adjusted_p_values,
        ):
            result["p_adj_holm"] = adjusted_p

            print(
                f"    {result['group1']} vs "
                f"{result['group2']}: "
                f"U = {result['U_statistic']:.4f}, "
                f"p = {result['p_value']:.6f}, "
                f"p_adj = {adjusted_p:.6f}"
            )

        pairwise_results.extend(
            trait_pairwise
        )

        print()

    return (
        pd.DataFrame(kruskal_results),
        pd.DataFrame(pairwise_results),
    )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main() -> None:
    """Run plant-trait treatment analysis."""
    args = parse_args()

    input_file = (
        args.input
        .expanduser()
        .resolve()
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

    kruskal_output = (
        output_dir
        / "trait_treatment_kruskalHolm.tsv"
    )

    pairwise_output = (
        output_dir
        / "trait_pairwise_testsHolm.tsv"
    )

    print("=== Plant trait treatment analysis ===")
    print(f"Python version:      {platform.python_version()}")
    print(f"pandas version:      {pd.__version__}")
    print(f"SciPy version:       {scipy.__version__}")
    print(f"statsmodels version: {statsmodels.__version__}")
    print(f"Input:               {input_file}")
    print(f"Output directory:    {output_dir}")
    print()

    traits = load_traits(
        input_file
    )

    print("Loaded data:")

    display_columns = [
        column
        for column in (
            "SRR",
            "Sample",
            "condition",
            *NUMERIC_TRAITS,
        )
        if column in traits.columns
    ]

    print(
        traits[display_columns]
        .to_string(index=False)
    )
    print()

    (
        kruskal_results,
        pairwise_results,
    ) = run_tests(traits)

    kruskal_results.to_csv(
        kruskal_output,
        sep="\t",
        index=False,
    )

    pairwise_results.to_csv(
        pairwise_output,
        sep="\t",
        index=False,
    )

    print(
        "Kruskal-Wallis results saved to: "
        f"{kruskal_output}"
    )
    print(
        "Pairwise test results saved to:   "
        f"{pairwise_output}"
    )


if __name__ == "__main__":
    main()