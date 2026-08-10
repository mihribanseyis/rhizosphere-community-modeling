#!/usr/bin/env python3
"""
Perform PCA on MICOM community-level exchange flux signatures.

This script reproduces the PCA workflow used for the thesis:

1. Reads a community exchange-flux signature table.
2. Retains community exchange reactions matching EX_*_m.
3. Removes all-NA variables.
4. Removes zero-variance variables.
5. Removes exact duplicate variables.
6. Excludes global activity metrics unless explicitly retained.
7. Z-score standardizes the remaining flux variables.
8. Performs PCA using NumPy singular value decomposition.
9. Tests treatment-associated differences for PC1-PC3 using
   Kruskal-Wallis tests with Benjamini-Hochberg FDR correction.
10. Optionally correlates PC1-PC4 with plant traits using Spearman
    rank correlation with Benjamini-Hochberg FDR correction.
11. Generates 2D PCA plots and, optionally, 3D PCA plots.

Important
---------
PCA is performed on the full set of eligible community exchange-flux
variables after the exclusions above.

It is NOT restricted to the FVA-selected reactions used for the
differential flux, flux-trait correlation, and regression analyses.

Usage
-----
PCA with treatment testing only:

python pca_fluxes.py \
    --flux qc_flux_signatures.tsv \
    --outdir pca_output

PCA with plant-trait correlations and 3D plots:

python pca_fluxes.py \
    --flux qc_flux_signatures.tsv \
    --traits plant_traits_clean.tsv \
    --outdir pca_output \
    --plot-3d

Input requirements
------------------
Flux table:
    Must contain:
        sample
        condition

    Community exchange variables are identified as:
        EX_*_m

Trait table:
    Must contain the sample identifier column and numeric trait columns.

Outputs
-------
pca_input_matrix.tsv
pca_input_matrix_scaled.tsv
pca_scores.tsv
pca_loadings.tsv
pca_explained_variance.tsv
pca_removed_columns.tsv
pc_treatment_kruskal.tsv
pc_trait_correlations.tsv          # only when --traits is supplied
pca_pc1_pc2.png
pca_pc1_pc2.pdf
pca_pc1_pc2_pc3.png               # only when --plot-3d
pca_pc1_pc2_pc3.pdf               # only when --plot-3d
summary.txt

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
from typing import List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.stats import kruskal, spearmanr


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Analysis constants
# ------------------------------------------------------------------------------

GLOBAL_ACTIVITY_METRICS = (
    "flux_nonzero",
    "flux_frac_nonzero",
)

TREATMENT_ORDER = (
    "CK",
    "NP",
    "NPM",
)

DEFAULT_PC_TRAIT_COUNT = 4
N_PCS_TREATMENT_TEST = 3


# ------------------------------------------------------------------------------
# Thesis plotting style
# ------------------------------------------------------------------------------

mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "Times",
            "Nimbus Roman",
        ],
        "font.size": 12,

        "axes.titlesize": 12,
        "axes.labelsize": 12,

        "xtick.labelsize": 12,
        "ytick.labelsize": 12,

        "legend.fontsize": 12,
        "legend.title_fontsize": 12,

        "figure.titlesize": 12,

        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
    }
)


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run PCA on MICOM community exchange-flux signatures, test "
            "PC1-PC3 across treatments, and optionally correlate PCs "
            "with plant traits."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--flux",
        required=True,
        type=Path,
        help="Community flux-signature table in TSV or CSV format.",
    )

    parser.add_argument(
        "--traits",
        required=False,
        type=Path,
        help=(
            "Optional plant-trait table. If omitted, PCA and treatment "
            "tests are still performed."
        ),
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--sample-col",
        default="sample",
        help="Sample identifier column.",
    )

    parser.add_argument(
        "--condition-col",
        default="condition",
        help="Treatment/condition column.",
    )

    parser.add_argument(
        "--sample-dir-col",
        default="sample_dir",
        help="Optional sample-directory metadata column.",
    )

    parser.add_argument(
        "--trait-cols",
        nargs="+",
        default=None,
        help=(
            "Plant trait columns to correlate with PCs. "
            "Default: all numeric trait columns."
        ),
    )

    parser.add_argument(
        "--pc-count",
        type=int,
        default=DEFAULT_PC_TRAIT_COUNT,
        help="Number of PCs used for PC-trait correlations.",
    )

    parser.add_argument(
        "--keep-global-metrics",
        action="store_true",
        help=(
            "Retain flux_nonzero and flux_frac_nonzero in PCA. "
            "By default these global activity metrics are excluded."
        ),
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure DPI.",
    )

    parser.add_argument(
        "--plot-3d",
        action="store_true",
        help="Generate PC1-PC2-PC3 3D PCA plots.",
    )

    parser.add_argument(
        "--elev",
        type=float,
        default=20.0,
        help="3D plot elevation angle.",
    )

    parser.add_argument(
        "--azim",
        type=float,
        default=45.0,
        help="3D plot azimuth angle.",
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# I/O helpers
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


def read_table_auto(
    path: Path,
) -> pd.DataFrame:
    """Read TSV, TXT, or CSV input based on extension."""
    suffix = path.suffix.lower()

    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(
            path,
            sep="\t",
        )

    if suffix == ".csv":
        return pd.read_csv(
            path
        )

    try:
        return pd.read_csv(
            path,
            sep="\t",
        )
    except Exception:
        return pd.read_csv(
            path
        )


# ------------------------------------------------------------------------------
# Statistical helpers
# ------------------------------------------------------------------------------

def bh_fdr(
    pvalues: List[float],
) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    pvalues_array = np.asarray(
        pvalues,
        dtype=float,
    )

    n_tests = len(
        pvalues_array
    )

    if n_tests == 0:
        return np.array(
            [],
            dtype=float,
        )

    order = np.argsort(
        pvalues_array
    )

    ranked = (
        pvalues_array[order]
    )

    adjusted = np.empty(
        n_tests,
        dtype=float,
    )

    previous = 1.0

    for index in range(
        n_tests - 1,
        -1,
        -1,
    ):
        value = (
            ranked[index]
            * n_tests
            / (index + 1)
        )

        previous = min(
            previous,
            value,
        )

        adjusted[index] = min(
            previous,
            1.0,
        )

    output = np.empty(
        n_tests,
        dtype=float,
    )

    output[order] = (
        adjusted
    )

    return output


def zscore_df(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Column-wise z-score standardization using ddof=0."""
    means = dataframe.mean(
        axis=0
    )

    stds = dataframe.std(
        axis=0,
        ddof=0,
    )

    if (stds == 0).any():
        zero_columns = (
            stds[
                stds == 0
            ]
            .index
            .tolist()
        )

        raise ValueError(
            "Zero-standard-deviation columns remain after filtering: "
            f"{zero_columns}"
        )

    return (
        dataframe - means
    ) / stds


# ------------------------------------------------------------------------------
# Flux-matrix preparation
# ------------------------------------------------------------------------------

def prepare_flux_matrix(
    flux_df: pd.DataFrame,
    sample_col: str,
    condition_col: str,
    sample_dir_col: str,
    keep_global_metrics: bool = False,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    List[str],
]:
    """
    Prepare the full community exchange-flux matrix for PCA.

    Community exchange reactions are identified as EX_*_m.

    Global activity metrics are excluded by default.
    """
    dataframe = flux_df.copy()

    required_columns = [
        sample_col,
        condition_col,
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in dataframe.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns in flux table: "
            f"{missing_columns}"
        )

    metadata_columns = [
        column
        for column in (
            sample_col,
            condition_col,
            sample_dir_col,
        )
        if column in dataframe.columns
    ]

    metadata = dataframe[
        metadata_columns
    ].copy()

    removed_columns: List[str] = []

    # --------------------------------------------------------------------------
    # Select full community exchange-flux space
    # --------------------------------------------------------------------------

    flux_columns = [
        column
        for column in dataframe.columns
        if (
            column.startswith("EX_")
            and column.endswith("_m")
        )
    ]

    if keep_global_metrics:
        for metric in GLOBAL_ACTIVITY_METRICS:
            if metric in dataframe.columns:
                flux_columns.append(
                    metric
                )

    if len(flux_columns) < 2:
        raise ValueError(
            "Fewer than 2 eligible community flux variables were found."
        )

    matrix = (
        dataframe[
            flux_columns
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
    )

    # --------------------------------------------------------------------------
    # Remove all-NA variables
    # --------------------------------------------------------------------------

    all_na_columns = (
        matrix.columns[
            matrix.isna().all()
        ]
        .tolist()
    )

    if all_na_columns:
        matrix = matrix.drop(
            columns=all_na_columns
        )

        removed_columns.extend(
            all_na_columns
        )

    # --------------------------------------------------------------------------
    # Remove zero-variance variables
    # --------------------------------------------------------------------------

    constant_columns = (
        matrix.columns[
            matrix.nunique(
                dropna=True
            ) <= 1
        ]
        .tolist()
    )

    if constant_columns:
        matrix = matrix.drop(
            columns=constant_columns
        )

        removed_columns.extend(
            constant_columns
        )

    # --------------------------------------------------------------------------
    # Remove exact duplicate variables
    #
    # Values are rounded to 12 decimal places before comparison, preserving
    # the original analysis implementation.
    # --------------------------------------------------------------------------

    duplicate_columns = []
    seen = {}

    for column in matrix.columns:
        key = tuple(
            np.round(
                matrix[column]
                .to_numpy(
                    dtype=float
                ),
                12,
            )
        )

        if key in seen:
            duplicate_columns.append(
                column
            )
        else:
            seen[key] = (
                column
            )

    if duplicate_columns:
        matrix = matrix.drop(
            columns=duplicate_columns
        )

        removed_columns.extend(
            duplicate_columns
        )

    # --------------------------------------------------------------------------
    # Missing-value validation
    # --------------------------------------------------------------------------

    if matrix.isna().any().any():
        na_columns = (
            matrix.columns[
                matrix.isna().any()
            ]
            .tolist()
        )

        raise ValueError(
            "Missing values remain in PCA matrix: "
            f"{na_columns}"
        )

    if matrix.shape[1] < 2:
        raise ValueError(
            "Fewer than 2 variables remained after PCA filtering."
        )

    # --------------------------------------------------------------------------
    # Align samples and metadata
    # --------------------------------------------------------------------------

    metadata = (
        metadata
        .drop_duplicates(
            subset=[sample_col]
        )
        .copy()
    )

    matrix = (
        pd.concat(
            [
                dataframe[
                    [sample_col]
                ],
                matrix,
            ],
            axis=1,
        )
        .drop_duplicates(
            subset=[sample_col]
        )
        .set_index(
            sample_col
        )
    )

    metadata = (
        metadata
        .set_index(
            sample_col
        )
        .loc[
            matrix.index
        ]
    )

    return (
        metadata,
        matrix,
        removed_columns,
    )


# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------

def run_pca(
    scaled_matrix: pd.DataFrame,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Run PCA using NumPy singular value decomposition."""
    matrix = scaled_matrix.to_numpy(
        dtype=float
    )

    if matrix.shape[0] < 2:
        raise ValueError(
            "PCA requires at least 2 samples."
        )

    u_matrix, singular_values, vt_matrix = np.linalg.svd(
        matrix,
        full_matrices=False,
    )

    scores = (
        u_matrix
        * singular_values
    )

    explained_variance = (
        singular_values ** 2
    ) / (
        matrix.shape[0] - 1
    )

    explained_ratio = (
        explained_variance
        / explained_variance.sum()
    )

    cumulative_ratio = np.cumsum(
        explained_ratio
    )

    score_columns = [
        f"PC{index}"
        for index in range(
            1,
            scores.shape[1] + 1,
        )
    ]

    loading_columns = [
        f"PC{index}"
        for index in range(
            1,
            vt_matrix.shape[0] + 1,
        )
    ]

    scores_df = pd.DataFrame(
        scores,
        index=scaled_matrix.index,
        columns=score_columns,
    )

    loadings_df = pd.DataFrame(
        vt_matrix.T,
        index=scaled_matrix.columns,
        columns=loading_columns,
    )

    variance_df = pd.DataFrame(
        {
            "PC": score_columns,
            "explained_variance": explained_variance,
            "explained_ratio": explained_ratio,
            "cumulative_ratio": cumulative_ratio,
        }
    )

    return (
        scores_df,
        loadings_df,
        variance_df,
    )


# ------------------------------------------------------------------------------
# Treatment tests
# ------------------------------------------------------------------------------

def kruskal_pcs(
    scores_df: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str,
    pcs: List[str],
) -> pd.DataFrame:
    """
    Test treatment-associated differences for selected PCs using
    Kruskal-Wallis followed by BH-FDR correction.
    """
    joined = scores_df.join(
        metadata[
            [condition_col]
        ]
    )

    results = []

    for pc in pcs:
        if pc not in joined.columns:
            continue

        groups = []
        labels = []

        for condition, subset in joined.groupby(
            condition_col
        ):
            values = (
                subset[pc]
                .dropna()
                .to_numpy()
            )

            if len(values) > 0:
                groups.append(
                    values
                )

                labels.append(
                    str(condition)
                )

        if len(groups) < 2:
            results.append(
                {
                    "PC": pc,
                    "H": np.nan,
                    "p_value": np.nan,
                    "groups": "",
                }
            )

            continue

        statistic, p_value = kruskal(
            *groups
        )

        results.append(
            {
                "PC": pc,
                "H": statistic,
                "p_value": p_value,
                "groups": ",".join(
                    labels
                ),
            }
        )

    output = pd.DataFrame(
        results
    )

    output["p_adj_bh"] = np.nan

    valid = (
        output["p_value"]
        .notna()
    )

    if valid.any():
        output.loc[
            valid,
            "p_adj_bh",
        ] = bh_fdr(
            output.loc[
                valid,
                "p_value",
            ].tolist()
        )

    return output


# ------------------------------------------------------------------------------
# Trait handling
# ------------------------------------------------------------------------------

def select_trait_columns(
    traits_df: pd.DataFrame,
    sample_col: str,
    requested: Optional[List[str]] = None,
) -> List[str]:
    """Determine which trait columns should be tested."""
    if sample_col not in traits_df.columns:
        raise ValueError(
            f"Trait table must contain '{sample_col}' column."
        )

    if requested is not None:
        missing_columns = [
            column
            for column in requested
            if column not in traits_df.columns
        ]

        if missing_columns:
            raise ValueError(
                "Requested trait columns were not found: "
                f"{missing_columns}"
            )

        return requested

    numeric_columns = [
        column
        for column in traits_df.columns
        if (
            column != sample_col
            and pd.api.types.is_numeric_dtype(
                traits_df[column]
            )
        )
    ]

    if not numeric_columns:
        raise ValueError(
            "No numeric trait columns available."
        )

    return numeric_columns


def correlate_pcs_with_traits(
    scores_df: pd.DataFrame,
    traits_df: pd.DataFrame,
    sample_col: str,
    trait_cols: List[str],
    pc_count: int,
) -> pd.DataFrame:
    """
    Calculate Spearman correlations between selected PCs and plant traits
    followed by BH-FDR correction.
    """
    pcs = [
        column
        for column in scores_df.columns
        if column.startswith("PC")
    ][:pc_count]

    merged = (
        scores_df
        .reset_index()
        .rename(
            columns={
                "index": sample_col
            }
        )
        .merge(
            traits_df[
                [sample_col]
                + trait_cols
            ],
            on=sample_col,
            how="inner",
        )
    )

    results = []

    for pc in pcs:
        for trait in trait_cols:
            subset = (
                merged[
                    [
                        pc,
                        trait,
                    ]
                ]
                .dropna()
            )

            n_samples = len(
                subset
            )

            if n_samples < 3:
                rho = np.nan
                p_value = np.nan
            else:
                rho, p_value = spearmanr(
                    subset[pc],
                    subset[trait],
                )

            results.append(
                {
                    "PC": pc,
                    "trait": trait,
                    "n": n_samples,
                    "rho": rho,
                    "p_value": p_value,
                }
            )

    output = pd.DataFrame(
        results
    )

    output["p_adj_bh"] = np.nan

    valid = (
        output["p_value"]
        .notna()
    )

    if valid.any():
        output.loc[
            valid,
            "p_adj_bh",
        ] = bh_fdr(
            output.loc[
                valid,
                "p_value",
            ].tolist()
        )

    return output.sort_values(
        [
            "p_adj_bh",
            "p_value",
            "PC",
            "trait",
        ],
        na_position="last",
    )


# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

def add_confidence_ellipse(
    axis: plt.Axes,
    x_values: np.ndarray,
    y_values: np.ndarray,
    n_std: float = 1.5,
    **kwargs,
) -> None:
    """Add covariance-based ellipse around a treatment group."""
    from matplotlib.patches import Ellipse

    if len(x_values) < 3:
        return

    covariance = np.cov(
        x_values,
        y_values,
    )

    if np.linalg.matrix_rank(
        covariance
    ) < 2:
        return

    eigenvalues, eigenvectors = np.linalg.eigh(
        covariance
    )

    order = (
        eigenvalues
        .argsort()[::-1]
    )

    eigenvalues = (
        eigenvalues[order]
    )

    eigenvectors = (
        eigenvectors[:, order]
    )

    angle = np.degrees(
        np.arctan2(
            *eigenvectors[:, 0][::-1]
        )
    )

    width, height = (
        2
        * n_std
        * np.sqrt(
            eigenvalues
        )
    )

    ellipse = Ellipse(
        xy=(
            np.mean(
                x_values
            ),
            np.mean(
                y_values
            ),
        ),
        width=width,
        height=height,
        angle=angle,
        fill=False,
        **kwargs,
    )

    axis.add_patch(
        ellipse
    )


def plot_pca_2d(
    scores_df: pd.DataFrame,
    metadata: pd.DataFrame,
    variance_df: pd.DataFrame,
    condition_col: str,
    out_png: Path,
    out_pdf: Path,
    dpi: int,
) -> None:
    """Generate PC1 vs PC2 PCA plot."""
    plot_df = scores_df.join(
        metadata[
            [condition_col]
        ]
    )

    colormap = plt.get_cmap(
        "tab10"
    )

    color_map = {
        "CK": colormap(0),
        "NP": colormap(1),
        "NPM": colormap(2),
    }

    unique_conditions = [
        condition
        for condition in TREATMENT_ORDER
        if condition in plot_df[
            condition_col
        ].unique()
    ]

    figure, axis = plt.subplots(
        figsize=(
            6.3,
            5.0,
        )
    )

    for condition in unique_conditions:
        subset = plot_df[
            plot_df[
                condition_col
            ] == condition
        ]

        axis.scatter(
            subset["PC1"],
            subset["PC2"],
            s=70,
            label=condition,
            color=color_map[
                condition
            ],
            edgecolor="black",
            linewidth=0.4,
            alpha=0.95,
            zorder=3,
        )

        add_confidence_ellipse(
            axis,
            subset[
                "PC1"
            ].to_numpy(),
            subset[
                "PC2"
            ].to_numpy(),
            edgecolor=color_map[
                condition
            ],
            linewidth=1.5,
            alpha=0.20,
        )

    pc1_percent = (
        variance_df.loc[
            variance_df["PC"]
            == "PC1",
            "explained_ratio",
        ]
        .iloc[0]
        * 100
    )

    pc2_percent = (
        variance_df.loc[
            variance_df["PC"]
            == "PC2",
            "explained_ratio",
        ]
        .iloc[0]
        * 100
    )

    axis.set_xlabel(
        f"PC1 ({pc1_percent:.1f}% explained variance)"
    )

    axis.set_ylabel(
        f"PC2 ({pc2_percent:.1f}% explained variance)"
    )

    axis.axhline(
        0,
        color="lightgray",
        linewidth=0.8,
        zorder=0,
    )

    axis.axvline(
        0,
        color="lightgray",
        linewidth=0.8,
        zorder=0,
    )

    axis.legend(
        frameon=False,
        loc="best",
        handletextpad=0.4,
        borderpad=0.2,
    )

    axis.grid(
        False
    )

    figure.tight_layout()

    figure.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
    )

    figure.savefig(
        out_pdf,
        bbox_inches="tight",
        transparent=False,
    )

    plt.close(
        figure
    )


def plot_pca_3d(
    scores_df: pd.DataFrame,
    metadata: pd.DataFrame,
    variance_df: pd.DataFrame,
    condition_col: str,
    out_png: Path,
    out_pdf: Path,
    dpi: int,
    elev: float,
    azim: float,
) -> None:
    """Generate PC1-PC2-PC3 3D PCA plot and rotated views."""
    required_pcs = {
        "PC1",
        "PC2",
        "PC3",
    }

    if not required_pcs.issubset(
        scores_df.columns
    ):
        missing_pcs = (
            required_pcs
            - set(
                scores_df.columns
            )
        )

        print(
            "WARNING: cannot plot 3D PCA; missing PCs: "
            f"{missing_pcs}"
        )

        return

    plot_df = scores_df.join(
        metadata[
            [condition_col]
        ]
    )

    unique_conditions = [
        condition
        for condition in TREATMENT_ORDER
        if condition in plot_df[
            condition_col
        ].unique()
    ]

    colormap = plt.get_cmap(
        "tab10"
    )

    color_map = {
        "CK": colormap(0),
        "NP": colormap(1),
        "NPM": colormap(2),
    }

    marker_map = {
        "CK": "o",
        "NP": "s",
        "NPM": "^",
    }

    def explained_percent(
        pc_name: str,
    ) -> float:
        row = variance_df.loc[
            variance_df["PC"]
            == pc_name,
            "explained_ratio",
        ]

        return (
            row.iloc[0] * 100
            if len(row) > 0
            else 0.0
        )

    figure = plt.figure(
        figsize=(
            9,
            7,
        )
    )

    axis = figure.add_subplot(
        111,
        projection="3d",
    )

    for condition in unique_conditions:
        subset = plot_df[
            plot_df[
                condition_col
            ] == condition
        ]

        axis.scatter(
            subset["PC1"],
            subset["PC2"],
            subset["PC3"],
            s=80,
            label=condition,
            color=color_map[
                condition
            ],
            marker=marker_map[
                condition
            ],
            depthshade=True,
            alpha=0.85,
        )

        for sample_id, row in subset.iterrows():
            axis.text(
                row["PC1"],
                row["PC2"],
                row["PC3"],
                s=str(
                    sample_id
                )[-3:],
                fontsize=6,
                color=color_map[
                    condition
                ],
                alpha=0.7,
            )

    axis.set_xlabel(
        f"PC1 ({explained_percent('PC1'):.1f}%)",
        labelpad=8,
    )

    axis.set_ylabel(
        f"PC2 ({explained_percent('PC2'):.1f}%)",
        labelpad=8,
    )

    axis.set_zlabel(
        f"PC3 ({explained_percent('PC3'):.1f}%)",
        labelpad=8,
    )

    axis.set_title(
        "3D PCA of microbial metabolic flux signatures"
    )

    axis.legend(
        frameon=False,
        loc="best",
        handletextpad=0.4,
        borderpad=0.2,
    )

    axis.view_init(
        elev=elev,
        azim=azim,
    )

    figure.tight_layout()

    figure.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
    )

    figure.savefig(
        out_pdf,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    rotated_views = (
        (30, 0),
        (30, 90),
        (30, 180),
        (90, 0),
    )

    for elevation, azimuth in rotated_views:
        rotated_figure = plt.figure(
            figsize=(
                9,
                7,
            )
        )

        rotated_axis = rotated_figure.add_subplot(
            111,
            projection="3d",
        )

        for condition in unique_conditions:
            subset = plot_df[
                plot_df[
                    condition_col
                ] == condition
            ]

            rotated_axis.scatter(
                subset["PC1"],
                subset["PC2"],
                subset["PC3"],
                s=80,
                label=condition,
                color=color_map[
                    condition
                ],
                marker=marker_map[
                    condition
                ],
                depthshade=True,
                alpha=0.85,
            )

        rotated_axis.set_xlabel(
            f"PC1 ({explained_percent('PC1'):.1f}%)"
        )

        rotated_axis.set_ylabel(
            f"PC2 ({explained_percent('PC2'):.1f}%)"
        )

        rotated_axis.set_zlabel(
            f"PC3 ({explained_percent('PC3'):.1f}%)"
        )

        rotated_axis.set_title(
            f"3D PCA — elev={elevation} azim={azimuth}"
        )

        rotated_axis.legend(
            frameon=True,
            loc="upper left",
        )

        rotated_axis.view_init(
            elev=elevation,
            azim=azimuth,
        )

        rotated_output = out_png.with_name(
            f"{out_png.stem}_elev{elevation}_azim{azimuth}.png"
        )

        rotated_figure.tight_layout()

        rotated_figure.savefig(
            rotated_output,
            dpi=dpi,
            bbox_inches="tight",
        )

        plt.close(
            rotated_figure
        )


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main() -> None:
    """Run PCA, treatment testing, and optional PC-trait correlations."""
    args = parse_args()

    if args.pc_count < 1:
        raise ValueError(
            "--pc-count must be at least 1."
        )

    flux_path = validate_file(
        args.flux,
        "Flux input file",
    )

    trait_path = None

    if args.traits is not None:
        trait_path = validate_file(
            args.traits,
            "Trait input file",
        )

    output_dir = (
        args.outdir
        .expanduser()
        .resolve()
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=== PCA of MICOM community exchange fluxes ===")
    print(f"Python version:     {platform.python_version()}")
    print(f"NumPy version:      {np.__version__}")
    print(f"pandas version:     {pd.__version__}")
    print(f"SciPy version:      {scipy.__version__}")
    print(f"Matplotlib version: {mpl.__version__}")
    print(f"Flux input:         {flux_path}")
    print(f"Trait input:        {trait_path or 'not provided'}")
    print(f"Output directory:   {output_dir}")
    print(
        "Global metrics:     "
        + (
            "retained"
            if args.keep_global_metrics
            else "excluded"
        )
    )
    print()

    # --------------------------------------------------------------------------
    # Prepare full exchange-flux matrix
    # --------------------------------------------------------------------------

    flux_df = read_table_auto(
        flux_path
    )

    (
        metadata,
        flux_matrix,
        removed_columns,
    ) = prepare_flux_matrix(
        flux_df=flux_df,
        sample_col=args.sample_col,
        condition_col=args.condition_col,
        sample_dir_col=args.sample_dir_col,
        keep_global_metrics=args.keep_global_metrics,
    )

    print(
        f"Samples used:       {flux_matrix.shape[0]}"
    )

    print(
        f"Flux variables:     {flux_matrix.shape[1]}"
    )

    print(
        f"Variables removed:  {len(removed_columns)}"
    )

    # --------------------------------------------------------------------------
    # Z-score standardization and PCA
    # --------------------------------------------------------------------------

    scaled_matrix = zscore_df(
        flux_matrix
    )

    (
        scores_df,
        loadings_df,
        variance_df,
    ) = run_pca(
        scaled_matrix
    )

    scores_with_condition = (
        scores_df.join(
            metadata[
                [args.condition_col]
            ]
        )
    )

    # --------------------------------------------------------------------------
    # Save PCA inputs and outputs
    # --------------------------------------------------------------------------

    flux_matrix.to_csv(
        output_dir
        / "pca_input_matrix.tsv",
        sep="\t",
    )

    scaled_matrix.to_csv(
        output_dir
        / "pca_input_matrix_scaled.tsv",
        sep="\t",
    )

    scores_with_condition.to_csv(
        output_dir
        / "pca_scores.tsv",
        sep="\t",
    )

    loadings_df.to_csv(
        output_dir
        / "pca_loadings.tsv",
        sep="\t",
    )

    variance_df.to_csv(
        output_dir
        / "pca_explained_variance.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        {
            "removed_column": removed_columns
        }
    ).to_csv(
        output_dir
        / "pca_removed_columns.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Treatment tests: first three PCs
    # --------------------------------------------------------------------------

    pcs_to_test = [
        column
        for column in scores_df.columns
        if column.startswith("PC")
    ][:N_PCS_TREATMENT_TEST]

    treatment_results = kruskal_pcs(
        scores_df=scores_df,
        metadata=metadata,
        condition_col=args.condition_col,
        pcs=pcs_to_test,
    )

    treatment_results.to_csv(
        output_dir
        / "pc_treatment_kruskal.tsv",
        sep="\t",
        index=False,
    )

    print()
    print(
        "Kruskal-Wallis tests for PC1-PC3:"
    )

    print(
        treatment_results.to_string(
            index=False
        )
    )

    # --------------------------------------------------------------------------
    # 2D PCA plot
    # --------------------------------------------------------------------------

    plot_pca_2d(
        scores_df=scores_df,
        metadata=metadata,
        variance_df=variance_df,
        condition_col=args.condition_col,
        out_png=output_dir
        / "pca_pc1_pc2.png",
        out_pdf=output_dir
        / "pca_pc1_pc2.pdf",
        dpi=args.dpi,
    )

    # --------------------------------------------------------------------------
    # Optional 3D PCA
    # --------------------------------------------------------------------------

    if args.plot_3d:
        plot_pca_3d(
            scores_df=scores_df,
            metadata=metadata,
            variance_df=variance_df,
            condition_col=args.condition_col,
            out_png=output_dir
            / "pca_pc1_pc2_pc3.png",
            out_pdf=output_dir
            / "pca_pc1_pc2_pc3.pdf",
            dpi=args.dpi,
            elev=args.elev,
            azim=args.azim,
        )

    # --------------------------------------------------------------------------
    # Optional PC-trait correlations: first four PCs by default
    # --------------------------------------------------------------------------

    if trait_path is not None:
        traits_df = read_table_auto(
            trait_path
        )

        trait_columns = select_trait_columns(
            traits_df=traits_df,
            sample_col=args.sample_col,
            requested=args.trait_cols,
        )

        correlation_results = correlate_pcs_with_traits(
            scores_df=scores_df,
            traits_df=traits_df,
            sample_col=args.sample_col,
            trait_cols=trait_columns,
            pc_count=args.pc_count,
        )

        correlation_results.to_csv(
            output_dir
            / "pc_trait_correlations.tsv",
            sep="\t",
            index=False,
        )

        print()
        print(
            "PC-trait correlations:"
        )

        print(
            correlation_results.to_string(
                index=False
            )
        )

    # --------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------

    def explained_percent(
        pc_name: str,
    ) -> float:
        row = variance_df.loc[
            variance_df["PC"]
            == pc_name,
            "explained_ratio",
        ]

        return (
            row.iloc[0] * 100
            if len(row) > 0
            else 0.0
        )

    cumulative_pc3 = (
        variance_df.loc[
            variance_df["PC"]
            == "PC3",
            "cumulative_ratio",
        ]
    )

    cumulative_pc3_percent = (
        cumulative_pc3.iloc[0] * 100
        if len(cumulative_pc3) > 0
        else 0.0
    )

    summary_lines = [
        "PCA analysis summary",
        "====================",
        f"Flux input:                {flux_path}",
        f"Trait input:               {trait_path or 'not provided'}",
        f"Output directory:          {output_dir}",
        f"Samples used:              {flux_matrix.shape[0]}",
        f"Variables used in PCA:     {flux_matrix.shape[1]}",
        f"Variables removed:         {len(removed_columns)}",
        f"PC1 explained variance:    {explained_percent('PC1'):.3f}%",
        f"PC2 explained variance:    {explained_percent('PC2'):.3f}%",
        f"PC3 explained variance:    {explained_percent('PC3'):.3f}%",
        f"PC1+PC2+PC3 cumulative:    {cumulative_pc3_percent:.3f}%",
    ]

    for _, row in treatment_results.iterrows():
        if pd.notna(
            row["p_value"]
        ):
            summary_lines.append(
                f"{row['PC']} Kruskal-Wallis: "
                f"H={row['H']:.4f}, "
                f"p={row['p_value']:.4g}, "
                f"FDR={row['p_adj_bh']:.4g}"
            )

    with (
        output_dir
        / "summary.txt"
    ).open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "\n".join(
                summary_lines
            )
            + "\n"
        )

    print()
    print(
        f"Done. Outputs written to: {output_dir}"
    )


if __name__ == "__main__":
    main()