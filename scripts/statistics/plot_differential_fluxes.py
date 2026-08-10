#!/usr/bin/env python3
"""
Plot treatment-associated differences in selected MICOM exchange fluxes.

The script reads:

1. Differential treatment-test results.
2. The FVA-filtered community exchange-flux table.

The top exchange reactions are selected according to:

    1. Lowest FDR-adjusted p-value
    2. Lowest raw p-value
    3. Highest Kruskal-Wallis H statistic

By default, the six highest-ranked exchange reactions are plotted.

Each panel shows:

- CK, NP, and NPM boxplots
- individual sample observations
- Kruskal-Wallis H statistic
- raw p-value
- FDR-adjusted p-value

Usage
-----
python plot_differential_fluxes.py \
    --diff-file differential_flux_treatment.tsv \
    --flux-file community_exchange_fluxes_selected.tsv \
    --output-dir figures

Outputs
-------
differential_exchange_fluxes.png
differential_exchange_fluxes.pdf
differential_exchange_fluxes_stats.tsv
differential_exchange_fluxes_selected.tsv

Dependencies
------------
- Python >= 3.8
- numpy
- pandas
- matplotlib
"""

from __future__ import annotations

import argparse
import math
import platform
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Reproducibility-critical plotting settings
# ------------------------------------------------------------------------------

TOP_N = 6

TREATMENT_ORDER = (
    "CK",
    "NP",
    "NPM",
)

JITTER_SEED = 42
JITTER_SD = 0.04

Y_UNIT = (
    r"Community exchange flux "
    r"(mmol gDW$^{-1}$ h$^{-1}$)"
)

COLORS = {
    "CK": "#4C72B0",
    "NP": "#DD8452",
    "NPM": "#55A868",
}


# ------------------------------------------------------------------------------
# Global figure formatting
# ------------------------------------------------------------------------------

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "Times",
            "Nimbus Roman No9 L",
            "Nimbus Roman",
            "DejaVu Serif",
        ],
        "mathtext.fontset": "stix",

        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,

        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 4,
        "ytick.major.size": 4,

        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        "figure.dpi": 300,
        "savefig.dpi": 300,
    }
)


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Plot the highest-ranked treatment-associated MICOM "
            "exchange fluxes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--diff-file",
        required=True,
        type=Path,
        help="Differential treatment-test results TSV.",
    )

    parser.add_argument(
        "--flux-file",
        required=True,
        type=Path,
        help="FVA-selected community exchange-flux TSV.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which figures and tables will be written.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=TOP_N,
        help=(
            "Number of highest-ranked exchange reactions to plot."
        ),
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

def choose_top_fluxes(
    diff_df: pd.DataFrame,
    top_n: int,
) -> list[str]:
    """
    Select top exchange fluxes according to differential-test ranking.

    Ranking is intentionally preserved from the original analysis:

    1. p_adj ascending
    2. p_value ascending
    3. H_statistic descending
    """
    required_columns = {
        "flux",
        "p_adj",
        "p_value",
        "H_statistic",
    }

    missing = (
        required_columns
        .difference(diff_df.columns)
    )

    if missing:
        raise ValueError(
            "Differential results file is missing columns: "
            f"{sorted(missing)}"
        )

    selected = diff_df.copy()

    selected = selected[
        selected["flux"]
        .astype(str)
        .str.startswith("EX_")
    ]

    selected = selected.sort_values(
        [
            "p_adj",
            "p_value",
            "H_statistic",
        ],
        ascending=[
            True,
            True,
            False,
        ],
        na_position="last",
    )

    return (
        selected["flux"]
        .drop_duplicates()
        .head(top_n)
        .tolist()
    )


def prettify_flux_name(
    flux: str,
) -> str:
    """Convert exchange reaction identifiers into readable labels."""
    labels = {
        "EX_h_m": "Proton exchange",
        "EX_h2o_m": "Water exchange",
        "EX_pi_m": "Phosphate uptake",
        "EX_nh4_m": "Ammonium uptake",
        "EX_ac_m": "Acetate exchange",
        "EX_o2_m": "Oxygen uptake",
        "EX_glc__D_m": "Glucose uptake",
        "EX_co2_m": "Carbon dioxide exchange",
        "EX_so4_m": "Sulfate uptake",
        "EX_gcald_m": "Glycolaldehyde exchange",
        "EX_cresol_m": "Cresol exchange",
    }

    return labels.get(
        flux,
        flux
        .replace("EX_", "")
        .replace("_m", "")
        .replace("__", "-")
        .replace("_", " ")
        .capitalize(),
    )


def p_to_text(
    p_value: float,
) -> str:
    """Format p-values for figure annotations."""
    if pd.isna(p_value):
        return r"$p = \mathrm{NA}$"

    if p_value < 0.001:
        return r"$p < 0.001$"

    return rf"$p = {p_value:.3f}$"


def style_axis(
    axis,
) -> None:
    """Apply consistent thesis figure styling."""
    for spine in axis.spines.values():
        spine.set_linewidth(
            1.0
        )

    axis.tick_params(
        axis="both",
        which="major",
        width=1.0,
        length=4,
    )

    axis.grid(
        False
    )


# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------

def validate_file(
    path: Path,
    description: str,
) -> Path:
    """Resolve and validate one input file."""
    path = (
        path
        .expanduser()
        .resolve()
    )

    if not path.is_file():
        raise FileNotFoundError(
            f"{description} does not exist: {path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"{description} is empty: {path}"
        )

    return path


# ------------------------------------------------------------------------------
# Main plotting workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Generate differential exchange-flux figure."""
    args = parse_args()

    if args.top_n < 1:
        raise ValueError(
            "--top-n must be at least 1."
        )

    diff_file = validate_file(
        args.diff_file,
        "Differential results file",
    )

    flux_file = validate_file(
        args.flux_file,
        "Flux table",
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

    print(
        "=== Differential exchange-flux figure ==="
    )
    print(
        f"Python version:     "
        f"{platform.python_version()}"
    )
    print(
        f"pandas version:     "
        f"{pd.__version__}"
    )
    print(
        f"NumPy version:      "
        f"{np.__version__}"
    )
    print(
        f"Matplotlib version: "
        f"{mpl.__version__}"
    )
    print(
        f"Differential file:  "
        f"{diff_file}"
    )
    print(
        f"Flux file:          "
        f"{flux_file}"
    )
    print(
        f"Output directory:   "
        f"{output_dir}"
    )
    print(
        f"Top N:              "
        f"{args.top_n}"
    )
    print(
        f"Jitter seed:        "
        f"{JITTER_SEED}"
    )
    print()

    # --------------------------------------------------------------------------
    # Read inputs
    # --------------------------------------------------------------------------

    diff_df = pd.read_csv(
        diff_file,
        sep="\t",
    )

    flux_df = pd.read_csv(
        flux_file,
        sep="\t",
    )

    if "condition" not in flux_df.columns:
        raise ValueError(
            "Flux table must contain a 'condition' column."
        )

    flux_df["condition"] = (
        flux_df["condition"]
        .astype(str)
        .str.strip()
    )

    unexpected_conditions = sorted(
        set(
            flux_df[
                "condition"
            ].dropna()
        ).difference(
            TREATMENT_ORDER
        )
    )

    if unexpected_conditions:
        print(
            "Warning: ignoring unexpected conditions: "
            f"{unexpected_conditions}"
        )

    top_fluxes = choose_top_fluxes(
        diff_df,
        args.top_n,
    )

    if not top_fluxes:
        raise ValueError(
            "No exchange fluxes were selected for plotting."
        )

    print(
        "Top exchange fluxes selected for plotting:"
    )

    for flux in top_fluxes:
        print(
            f" - {flux}"
        )

    # --------------------------------------------------------------------------
    # Figure dimensions
    # --------------------------------------------------------------------------

    n_fluxes = len(
        top_fluxes
    )

    ncols = 2

    nrows = math.ceil(
        n_fluxes
        / ncols
    )

    figure, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(
            9.0,
            3.5 * nrows,
        ),
        squeeze=False,
    )

    axes = (
        axes.flatten()
    )

    rng = np.random.default_rng(
        JITTER_SEED
    )

    stats_rows = []

    panel_labels = list(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    )

    # --------------------------------------------------------------------------
    # Plot each selected flux
    # --------------------------------------------------------------------------

    for index, flux in enumerate(
        top_fluxes
    ):
        axis = axes[index]

        if flux not in flux_df.columns:
            axis.set_visible(
                False
            )

            print(
                f"Warning: {flux} was selected but is "
                "not present in the flux table."
            )

            continue

        subset = flux_df[
            [
                "condition",
                flux,
            ]
        ].copy()

        subset[flux] = pd.to_numeric(
            subset[flux],
            errors="coerce",
        )

        subset = subset.dropna(
            subset=[
                "condition",
                flux,
            ]
        )

        grouped_data = [
            subset.loc[
                subset["condition"]
                == condition,
                flux,
            ].to_numpy()
            for condition
            in TREATMENT_ORDER
        ]

        empty_conditions = [
            condition
            for condition, values
            in zip(
                TREATMENT_ORDER,
                grouped_data,
            )
            if len(values) == 0
        ]

        if empty_conditions:
            axis.set_visible(
                False
            )

            print(
                f"Warning: {flux} has no valid values for: "
                f"{', '.join(empty_conditions)}"
            )

            continue

        # ----------------------------------------------------------------------
        # Boxplots
        # ----------------------------------------------------------------------

        boxplot = axis.boxplot(
            grouped_data,
            tick_labels=TREATMENT_ORDER,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            boxprops={
                "linewidth": 1.15,
                "edgecolor": "black",
            },
            whiskerprops={
                "linewidth": 1.15,
                "color": "black",
            },
            capprops={
                "linewidth": 1.15,
                "color": "black",
            },
            medianprops={
                "linewidth": 1.5,
                "color": "black",
            },
        )

        for box, condition in zip(
            boxplot["boxes"],
            TREATMENT_ORDER,
        ):
            box.set_facecolor(
                COLORS[condition]
            )

            box.set_alpha(
                0.28
            )

            box.set_edgecolor(
                "black"
            )

        # ----------------------------------------------------------------------
        # Individual observations
        # ----------------------------------------------------------------------

        for position, condition in enumerate(
            TREATMENT_ORDER,
            start=1,
        ):
            y_values = subset.loc[
                subset["condition"]
                == condition,
                flux,
            ].to_numpy()

            x_values = rng.normal(
                loc=position,
                scale=JITTER_SD,
                size=len(y_values),
            )

            axis.scatter(
                x_values,
                y_values,
                s=34,
                facecolor=COLORS[
                    condition
                ],
                edgecolor="black",
                linewidth=0.5,
                alpha=0.95,
                zorder=3,
            )

        # ----------------------------------------------------------------------
        # Statistical results
        # ----------------------------------------------------------------------

        result_rows = diff_df.loc[
            diff_df["flux"]
            == flux
        ]

        if result_rows.empty:
            h_statistic = np.nan
            p_value = np.nan
            p_adjusted = np.nan

        else:
            result_row = (
                result_rows.iloc[0]
            )

            h_statistic = pd.to_numeric(
                result_row[
                    "H_statistic"
                ],
                errors="coerce",
            )

            p_value = pd.to_numeric(
                result_row[
                    "p_value"
                ],
                errors="coerce",
            )

            p_adjusted = pd.to_numeric(
                result_row[
                    "p_adj"
                ],
                errors="coerce",
            )

        stats_rows.append(
            {
                "flux": flux,
                "label": prettify_flux_name(
                    flux
                ),
                "H_statistic": h_statistic,
                "p_value": p_value,
                "p_adj": p_adjusted,
            }
        )

        # ----------------------------------------------------------------------
        # Titles and annotations
        # ----------------------------------------------------------------------

        axis.set_title(
            prettify_flux_name(
                flux
            ),
            pad=7,
        )

        axis.text(
            -0.10,
            1.07,
            panel_labels[index],
            transform=axis.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

        if not pd.isna(
            p_adjusted
        ):
            annotation_text = (
                "Kruskal–Wallis\n"
                rf"$H = {h_statistic:.2f}$, "
                f"{p_to_text(p_value)}\n"
                rf"$p_{{\mathrm{{FDR}}}} = "
                f"{p_adjusted:.3f}$"
            )

        else:
            annotation_text = (
                "Kruskal–Wallis\n"
                rf"$H = {h_statistic:.2f}$, "
                f"{p_to_text(p_value)}\n"
                r"$p_{\mathrm{FDR}} = \mathrm{NA}$"
            )

        axis.text(
            0.04,
            0.96,
            annotation_text,
            transform=axis.transAxes,
            fontsize=8,
            va="top",
            ha="left",
            linespacing=1.15,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.85,
                "pad": 2.2,
            },
            zorder=5,
        )

        # ----------------------------------------------------------------------
        # Axis labels
        # ----------------------------------------------------------------------

        row_index = (
            index
            // ncols
        )

        is_bottom_row = (
            row_index
            == nrows - 1
        )

        if is_bottom_row:
            axis.set_xlabel(
                "Treatment"
            )
        else:
            axis.set_xlabel(
                ""
            )

        axis.set_ylabel(
            Y_UNIT
        )

        style_axis(
            axis
        )

    # Hide any unused axes.
    for axis in axes[
        n_fluxes:
    ]:
        axis.set_visible(
            False
        )

    # --------------------------------------------------------------------------
    # Layout and output
    # --------------------------------------------------------------------------

    figure.subplots_adjust(
        left=0.10,
        right=0.98,
        bottom=0.08,
        top=0.97,
        wspace=0.34,
        hspace=0.38,
    )

    png_path = (
        output_dir
        / "differential_exchange_fluxes.png"
    )

    pdf_path = (
        output_dir
        / "differential_exchange_fluxes.pdf"
    )

    stats_path = (
        output_dir
        / "differential_exchange_fluxes_stats.tsv"
    )

    selected_path = (
        output_dir
        / "differential_exchange_fluxes_selected.tsv"
    )

    figure.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    figure.savefig(
        pdf_path,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(
        figure
    )

    pd.DataFrame(
        stats_rows
    ).to_csv(
        stats_path,
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        {
            "flux": top_fluxes,
            "label": [
                prettify_flux_name(
                    flux
                )
                for flux
                in top_fluxes
            ],
        }
    ).to_csv(
        selected_path,
        sep="\t",
        index=False,
    )

    print(
        f"Saved: {png_path}"
    )
    print(
        f"Saved: {pdf_path}"
    )
    print(
        f"Saved stats table: "
        f"{stats_path}"
    )
    print(
        f"Saved selected fluxes: "
        f"{selected_path}"
    )


if __name__ == "__main__":
    main()