#!/usr/bin/env python3
"""
Fit treatment-adjusted regression models for maize yield.

This script evaluates whether associations between community-level nutrient
uptake fluxes and maize yield remain after accounting for fertilization
treatment.

Three ordinary least-squares models are fitted:

    1. Treatment only
       maize_yield ~ C(condition)

    2. Treatment + phosphate uptake
       maize_yield ~ C(condition) + phosphate_uptake

    3. Treatment + ammonium uptake
       maize_yield ~ C(condition) + ammonium_uptake

The contribution of each flux is assessed by comparing the treatment-only
model with the corresponding treatment-plus-flux model using a nested
sequential F-test implemented with statsmodels ANOVA.

Separate phosphate and ammonium models are fitted to avoid including both
moderately correlated nutrient-uptake predictors in the same model.

Usage
-----
python regression_treatment_adjusted.py \
    --input analysis_table.tsv \
    --output-dir results

Input requirements
------------------
The input TSV must contain:

    condition
    maize_yield
    phosphate_uptake
    ammonium_uptake

Treatment levels are expected to be:

    CK
    NP
    NPM

Outputs
-------
regression_adjusted_terms.tsv
regression_adjusted_model_comparison.tsv
regression_model_comparison_phosphate.tsv
regression_model_comparison_ammonium.tsv
regression_adjusted_full_output.txt
regression_adjusted_model_comparison.tex

Dependencies
------------
- Python >= 3.8
- pandas
- statsmodels
"""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path

import pandas as pd
import statsmodels
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Analysis constants
# ------------------------------------------------------------------------------

TREATMENT_ORDER = (
    "CK",
    "NP",
    "NPM",
)

FORMULAS = {
    "treatment_only":
        "maize_yield ~ C(condition)",

    "treatment_plus_phosphate":
        "maize_yield ~ C(condition) + phosphate_uptake",

    "treatment_plus_ammonium":
        "maize_yield ~ C(condition) + ammonium_uptake",
}

ROW_LABELS = {
    "treatment_only":
        "Treatment only",

    "treatment_plus_phosphate":
        "+ Phosphate uptake",

    "treatment_plus_ammonium":
        "+ Ammonium uptake",
}


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Fit treatment-adjusted OLS models for maize yield and compare "
            "treatment-only models with treatment-plus-flux models."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help=(
            "Input analysis table containing condition, maize_yield, "
            "phosphate_uptake, and ammonium_uptake."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which regression outputs will be written.",
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


def load_analysis_table(
    input_file: Path,
) -> pd.DataFrame:
    """Read and validate the regression analysis table."""
    dataframe = pd.read_csv(
        input_file,
        sep="\t",
    )

    dataframe.columns = [
        column.strip()
        for column in dataframe.columns
    ]

    required_columns = {
        "condition",
        "maize_yield",
        "phosphate_uptake",
        "ammonium_uptake",
    }

    missing = (
        required_columns
        .difference(dataframe.columns)
    )

    if missing:
        raise ValueError(
            "Input table is missing required columns: "
            f"{sorted(missing)}"
        )

    dataframe["condition"] = pd.Categorical(
        dataframe["condition"],
        categories=list(TREATMENT_ORDER),
        ordered=False,
    )

    for column in (
        "maize_yield",
        "phosphate_uptake",
        "ammonium_uptake",
    ):
        dataframe[column] = pd.to_numeric(
            dataframe[column],
            errors="coerce",
        )

    return dataframe


# ------------------------------------------------------------------------------
# Model fitting
# ------------------------------------------------------------------------------

def fit_models(
    dataframe: pd.DataFrame,
):
    """Fit all predefined OLS regression models."""
    return {
        name: smf.ols(
            formula=formula,
            data=dataframe,
        ).fit()

        for name, formula
        in FORMULAS.items()
    }


# ------------------------------------------------------------------------------
# Model-level summary
# ------------------------------------------------------------------------------

def build_model_comparison(
    models,
) -> pd.DataFrame:
    """Build model-level comparison table."""
    base_r2 = (
        models["treatment_only"]
        .rsquared
    )

    base_adjusted_r2 = (
        models["treatment_only"]
        .rsquared_adj
    )

    rows = []

    for name, model in models.items():
        rows.append(
            {
                "model": name,
                "formula": FORMULAS[name],
                "n": int(model.nobs),
                "df_model": model.df_model,
                "df_resid": model.df_resid,
                "r_squared": model.rsquared,
                "adj_r_squared": model.rsquared_adj,
                "delta_r_squared_vs_treatment_only":
                    model.rsquared - base_r2,
                "delta_adj_r_squared_vs_treatment_only":
                    model.rsquared_adj - base_adjusted_r2,
                "aic": model.aic,
                "bic": model.bic,
                "f_statistic": model.fvalue,
                "model_p_value": model.f_pvalue,
            }
        )

    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------
# Term-level results
# ------------------------------------------------------------------------------

def build_term_table(
    models,
) -> pd.DataFrame:
    """Build coefficient and confidence-interval table."""
    rows = []

    for name, model in models.items():
        confidence = model.conf_int()

        confidence.columns = [
            "ci_lower",
            "ci_upper",
        ]

        for term in model.params.index:
            rows.append(
                {
                    "model": name,
                    "formula": FORMULAS[name],
                    "term": term,
                    "estimate": model.params[term],
                    "ci_lower":
                        confidence.loc[
                            term,
                            "ci_lower",
                        ],
                    "ci_upper":
                        confidence.loc[
                            term,
                            "ci_upper",
                        ],
                    "p_value":
                        model.pvalues[term],
                    "n":
                        int(model.nobs),
                    "df_resid":
                        model.df_resid,
                    "r_squared":
                        model.rsquared,
                    "adj_r_squared":
                        model.rsquared_adj,
                    "aic":
                        model.aic,
                    "bic":
                        model.bic,
                }
            )

    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------
# Formatting helpers
# ------------------------------------------------------------------------------

def format_p(
    value,
) -> str:
    """Format p-values for the LaTeX table."""
    if pd.isna(value):
        return "--"

    if value < 0.001:
        return "$<0.001$"

    return f"{value:.3f}"


def format_number(
    value,
    digits: int = 3,
) -> str:
    """Format numeric values for the LaTeX table."""
    if pd.isna(value):
        return "--"

    return f"{value:.{digits}f}"


# ------------------------------------------------------------------------------
# LaTeX output
# ------------------------------------------------------------------------------

def write_latex_table(
    model_comparison: pd.DataFrame,
    anova_phosphate: pd.DataFrame,
    anova_ammonium: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write main-text model comparison table in LaTeX."""
    flux_pvalues = {
        "treatment_only": None,

        "treatment_plus_phosphate":
            float(
                anova_phosphate.loc[
                    1,
                    "Pr(>F)",
                ]
            ),

        "treatment_plus_ammonium":
            float(
                anova_ammonium.loc[
                    1,
                    "Pr(>F)",
                ]
            ),
    }

    rows = []

    for _, row in model_comparison.iterrows():
        model_name = (
            row["model"]
        )

        rows.append(
            [
                ROW_LABELS[
                    model_name
                ],
                format_number(
                    row[
                        "adj_r_squared"
                    ]
                ),
                format_number(
                    row[
                        "delta_adj_r_squared_vs_treatment_only"
                    ]
                ),
                format_p(
                    flux_pvalues[
                        model_name
                    ]
                ),
            ]
        )

    latex_lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        (
            r"\caption{Comparison of treatment-adjusted regression models "
            r"for maize yield. The added flux p-value corresponds to the "
            r"nested ANOVA comparing the treatment-only model with the "
            r"corresponding treatment-plus-flux model.}"
        ),
        r"\label{tab:adjusted_regression}",
        r"\begin{tabular}{lrrr}",
        r"\toprule",
        (
            r"Model comparison & Adjusted $R^2$ & "
            r"$\Delta$ adjusted $R^2$ & Added flux $p$ \\"
        ),
        r"\midrule",
    ]

    for row in rows:
        latex_lines.append(
            f"{row[0]} & "
            f"{row[1]} & "
            f"{row[2]} & "
            f"{row[3]} \\\\"
        )

    latex_lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
        ]
    )

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "\n".join(
                latex_lines
            )
        )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main() -> None:
    """Run treatment-adjusted regression analysis."""
    args = parse_args()

    input_file = validate_file(
        args.input,
        "Regression analysis table",
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

    output_terms = (
        output_dir
        / "regression_adjusted_terms.tsv"
    )

    output_model_comparison = (
        output_dir
        / "regression_adjusted_model_comparison.tsv"
    )

    output_anova_phosphate = (
        output_dir
        / "regression_model_comparison_phosphate.tsv"
    )

    output_anova_ammonium = (
        output_dir
        / "regression_model_comparison_ammonium.tsv"
    )

    output_text = (
        output_dir
        / "regression_adjusted_full_output.txt"
    )

    output_latex = (
        output_dir
        / "regression_adjusted_model_comparison.tex"
    )

    print(
        "=== Treatment-adjusted regression analysis ==="
    )
    print(
        f"Python version:      "
        f"{platform.python_version()}"
    )
    print(
        f"pandas version:      "
        f"{pd.__version__}"
    )
    print(
        f"statsmodels version: "
        f"{statsmodels.__version__}"
    )
    print(
        f"Input:               "
        f"{input_file}"
    )
    print(
        f"Output directory:    "
        f"{output_dir}"
    )
    print()

    # --------------------------------------------------------------------------
    # Load data
    # --------------------------------------------------------------------------

    dataframe = load_analysis_table(
        input_file
    )

    # --------------------------------------------------------------------------
    # Fit OLS models
    # --------------------------------------------------------------------------

    models = fit_models(
        dataframe
    )

    # --------------------------------------------------------------------------
    # Model-level comparison
    # --------------------------------------------------------------------------

    model_comparison = build_model_comparison(
        models
    )

    model_comparison.to_csv(
        output_model_comparison,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Term-level coefficients
    # --------------------------------------------------------------------------

    terms = build_term_table(
        models
    )

    terms.to_csv(
        output_terms,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Nested model comparisons
    # --------------------------------------------------------------------------

    anova_phosphate = anova_lm(
        models[
            "treatment_only"
        ],
        models[
            "treatment_plus_phosphate"
        ],
    )

    anova_ammonium = anova_lm(
        models[
            "treatment_only"
        ],
        models[
            "treatment_plus_ammonium"
        ],
    )

    anova_phosphate.to_csv(
        output_anova_phosphate,
        sep="\t",
    )

    anova_ammonium.to_csv(
        output_anova_ammonium,
        sep="\t",
    )

    # --------------------------------------------------------------------------
    # Full text output
    # --------------------------------------------------------------------------

    full_output = []

    for name, model in models.items():
        full_output.append(
            "=" * 80
        )

        full_output.append(
            name
        )

        full_output.append(
            FORMULAS[name]
        )

        full_output.append(
            "=" * 80
        )

        full_output.append(
            str(
                model.summary()
            )
        )

        full_output.append(
            "\nType II ANOVA"
        )

        full_output.append(
            str(
                anova_lm(
                    model,
                    typ=2,
                )
            )
        )

        full_output.append(
            "\n"
        )

    full_output.append(
        "=" * 80
    )

    full_output.append(
        "Nested model comparison: "
        "treatment_only vs treatment_plus_phosphate"
    )

    full_output.append(
        "=" * 80
    )

    full_output.append(
        str(
            anova_phosphate
        )
    )

    full_output.append(
        "\n"
    )

    full_output.append(
        "=" * 80
    )

    full_output.append(
        "Nested model comparison: "
        "treatment_only vs treatment_plus_ammonium"
    )

    full_output.append(
        "=" * 80
    )

    full_output.append(
        str(
            anova_ammonium
        )
    )

    full_output.append(
        "\n"
    )

    with output_text.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "\n".join(
                full_output
            )
        )

    # --------------------------------------------------------------------------
    # LaTeX table
    # --------------------------------------------------------------------------

    write_latex_table(
        model_comparison=model_comparison,
        anova_phosphate=anova_phosphate,
        anova_ammonium=anova_ammonium,
        output_path=output_latex,
    )

    # --------------------------------------------------------------------------
    # Console summary
    # --------------------------------------------------------------------------

    print(
        f"Wrote term-level table: "
        f"{output_terms}"
    )

    print(
        f"Wrote model comparison table: "
        f"{output_model_comparison}"
    )

    print(
        f"Wrote phosphate nested comparison: "
        f"{output_anova_phosphate}"
    )

    print(
        f"Wrote ammonium nested comparison: "
        f"{output_anova_ammonium}"
    )

    print(
        f"Wrote full output: "
        f"{output_text}"
    )

    print(
        f"Wrote LaTeX table: "
        f"{output_latex}"
    )

    print()
    print(
        "Model comparison:"
    )

    print(
        model_comparison.to_string(
            index=False
        )
    )

    print()
    print(
        "Nested comparison: "
        "treatment only vs treatment + phosphate"
    )

    print(
        anova_phosphate
    )

    print()
    print(
        "Nested comparison: "
        "treatment only vs treatment + ammonium"
    )

    print(
        anova_ammonium
    )


if __name__ == "__main__":
    main()