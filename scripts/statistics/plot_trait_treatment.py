#!/usr/bin/env python3
"""
Plot plant-trait distributions across fertilization treatments.

The script generates a two-panel figure showing:

A. Maize yield
B. Total acid content

for CK, NP, and NPM treatments.

Each panel contains:

- treatment-colored boxplots
- individual observations with reproducible horizontal jitter
- Kruskal-Wallis H statistic and p-value

The figure is written in both PDF and PNG formats.

Usage
-----
python plot_trait_treatment.py \
    --input plant_trait.tsv \
    --output-dir figures

Outputs
-------
plant_traits_combined.pdf
plant_traits_combined.png

Input requirements
------------------
The input TSV must contain:

    Sample
    Maize yield
    Total acid

Dependencies
------------
- Python >= 3.8
- numpy
- pandas
- scipy
- matplotlib
"""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import kruskal


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Analysis and plotting constants
# ------------------------------------------------------------------------------

TREATMENT_ORDER = (
    "CK",
    "NP",
    "NPM",
)

NUMERIC_TRAITS = (
    "Maize yield",
    "Total acid",
)

PANEL_LABELS = (
    "A",
    "B",
)

CONDITION_PATTERN = r"(NPM|NP|CK)"

JITTER_SEED = 42
JITTER_SD = 0.04

FIGURE_SIZE = (
    10.5,
    4.6,
)


# ------------------------------------------------------------------------------
# Global figure formatting
#
# These settings intentionally preserve the thesis figure formatting.
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
        "axes.labelsize": 10,
        "axes.titlesize": 10,
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
# Treatment colors
#
# These are fixed intentionally so figures remain visually consistent across
# analyses.
# ------------------------------------------------------------------------------

COLORS = {
    "CK": "#4C72B0",
    "NP": "#DD8452",
    "NPM": "#55A868",
}


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Plot maize yield and total acid distributions across "
            "CK, NP, and NPM treatments."
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
        help="Directory in which figure files will be written.",
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def p_to_text(p_value: float) -> str:
    """Format p-values for figure annotations."""
    if p_value < 0.001:
        return r"$p < 0.001$"

    return rf"$p = {p_value:.3f}$"


def load_traits(
    input_file: Path,
) -> pd.DataFrame:
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
                f"Input file must contain a "
                f"'{trait}' column."
            )

        traits[trait] = pd.to_numeric(
            traits[trait],
            errors="coerce",
        )

    traits = (
        traits
        .dropna(subset=["condition"])
        .copy()
    )

    return traits


# ------------------------------------------------------------------------------
# Main plotting workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Generate the combined plant-trait figure."""
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

    print("=== Plant trait figure ===")
    print(f"Python version:     {platform.python_version()}")
    print(f"pandas version:     {pd.__version__}")
    print(f"NumPy version:      {np.__version__}")
    print(f"SciPy version:      {scipy.__version__}")
    print(f"Matplotlib version: {mpl.__version__}")
    print(f"Input:              {input_file}")
    print(f"Output directory:   {output_dir}")
    print(f"Jitter seed:        {JITTER_SEED}")
    print()

    traits = load_traits(
        input_file
    )

    fig, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=FIGURE_SIZE,
    )

    # Fixed seed preserves the exact random jitter pattern between runs.
    rng = np.random.default_rng(
        JITTER_SEED
    )

    for ax, trait, panel_label in zip(
        axes,
        NUMERIC_TRAITS,
        PANEL_LABELS,
    ):
        sub = (
            traits[
                ["condition", trait]
            ]
            .dropna()
            .copy()
        )

        grouped_data = [
            sub.loc[
                sub["condition"] == condition,
                trait,
            ].to_numpy()
            for condition in TREATMENT_ORDER
        ]

        empty_groups = [
            condition
            for condition, values in zip(
                TREATMENT_ORDER,
                grouped_data,
            )
            if len(values) == 0
        ]

        if empty_groups:
            raise ValueError(
                f"No valid values found for "
                f"'{trait}' in: "
                f"{', '.join(empty_groups)}"
            )

        statistic, p_kw = kruskal(
            *grouped_data
        )

        # ----------------------------------------------------------------------
        # Boxplots
        # ----------------------------------------------------------------------

        boxplot = ax.boxplot(
            grouped_data,
            tick_labels=TREATMENT_ORDER,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            boxprops={
                "linewidth": 1.2,
                "edgecolor": "black",
            },
            whiskerprops={
                "linewidth": 1.2,
                "color": "black",
            },
            capprops={
                "linewidth": 1.2,
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
            y_values = sub.loc[
                sub["condition"] == condition,
                trait,
            ].to_numpy()

            x_values = rng.normal(
                loc=position,
                scale=JITTER_SD,
                size=len(y_values),
            )

            ax.scatter(
                x_values,
                y_values,
                s=36,
                facecolor=COLORS[condition],
                edgecolor="black",
                linewidth=0.55,
                alpha=0.95,
                zorder=3,
            )

        # ----------------------------------------------------------------------
        # Axis labels
        # ----------------------------------------------------------------------

        ax.set_xlabel(
            "Treatment"
        )

        if trait == "Maize yield":
            ax.set_ylabel(
                r"Maize yield (Mg ha$^{-1}$)"
            )
        else:
            ax.set_ylabel(
                "Total acid content"
            )

        # ----------------------------------------------------------------------
        # Panel label
        # ----------------------------------------------------------------------

        ax.text(
            -0.12,
            1.08,
            panel_label,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

        # ----------------------------------------------------------------------
        # Statistical annotation
        #
        # Panel A: upper left
        # Panel B: upper right
        # ----------------------------------------------------------------------

        if panel_label == "A":
            annotation_x = 0.04
            annotation_alignment = "left"
        else:
            annotation_x = 0.96
            annotation_alignment = "right"

        ax.text(
            annotation_x,
            0.96,
            (
                "Kruskal–Wallis\n"
                rf"$H = {statistic:.2f}$, "
                f"{p_to_text(p_kw)}"
            ),
            transform=ax.transAxes,
            fontsize=8.5,
            va="top",
            ha=annotation_alignment,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.85,
                "pad": 2.5,
            },
            zorder=5,
        )

        # ----------------------------------------------------------------------
        # Axis formatting
        # ----------------------------------------------------------------------

        for spine in ax.spines.values():
            spine.set_linewidth(
                1.0
            )

        ax.tick_params(
            axis="both",
            which="major",
            width=1.0,
            length=4,
        )

        ax.grid(
            False
        )

    # --------------------------------------------------------------------------
    # Layout and output
    # --------------------------------------------------------------------------

    fig.tight_layout()

    pdf_path = (
        output_dir
        / "plant_traits_combined.pdf"
    )

    png_path = (
        output_dir
        / "plant_traits_combined.png"
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        facecolor="white",
    )

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(
        fig
    )

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()