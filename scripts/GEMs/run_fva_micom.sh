#!/usr/bin/env python3
"""
Run MICOM-conditioned flux variability analysis.

This script performs flux variability analysis (FVA) on serialized MICOM
community models while preserving the cooperative-tradeoff solution according
to one of three assumptions.

Assumptions
-----------
objective_only
    Run MICOM cooperative tradeoff to obtain a reference solution, but do not
    explicitly fix community or taxon growth before FVA.

    This is the broadest feasible-space interpretation and is closest to
    conventional COBRA-style FVA on the MICOM community model.

fix_community_growth
    Fix the MICOM community growth variable to the cooperative-tradeoff
    reference value, within the specified tolerance.

    Individual taxon growth rates remain unconstrained.

fix_taxon_growth
    Fix both community growth and each taxon's growth rate to the MICOM
    cooperative-tradeoff reference solution, within the specified tolerance.

    This is the strictest interpretation and is closest to performing FVA
    around the original MICOM tradeoff solution.

Reaction-selection modes
------------------------
active
    Use the reactions selected in a cross_check_reaction_summary.tsv file.

    The file must contain reaction identifiers as its row index and a Boolean
    selection column, such as selected_by_threshold.

all
    Analyze all community medium exchange reactions whose identifiers begin
    with EX_ and end with _m.

internal
    Analyze all reactions except community medium exchanges.

Expected MICOM directory structure
----------------------------------
BASE_PATH/
├── CK/
│   ├── SRR16095329/
│   │   ├── community.pickle
│   │   └── medium_applied.tsv
│   └── ...
├── NP/
│   └── ...
└── NPM/
    └── ...

Usage
-----
Active-reaction FVA:

python run_fva_micom.py \
    --base_path /path/to/micom/output \
    --output_dir /path/to/fva/results \
    --selection_summary /path/to/cross_check_reaction_summary.tsv \
    --selection_col selected_by_threshold \
    --assumption fix_taxon_growth \
    --reactions active \
    --tradeoff 0.5 \
    --fraction 1.0 \
    --solver gurobi

Run all community exchange reactions:

python run_fva_micom.py \
    --base_path /path/to/micom/output \
    --output_dir /path/to/fva/results \
    --assumption fix_taxon_growth \
    --reactions all \
    --tradeoff 0.5 \
    --fraction 1.0 \
    --solver gurobi

Run one sample:

python run_fva_micom.py \
    --base_path /path/to/micom/output \
    --output_dir /path/to/fva/results \
    --selection_summary /path/to/cross_check_reaction_summary.tsv \
    --only_sample CK_SRR16095329

Outputs
-------
One CSV file is written per sample:

    OUTPUT_DIR/CK_329.csv
    OUTPUT_DIR/NP_324.csv
    OUTPUT_DIR/NPM_333.csv

Each sample file contains:

    minimum
    maximum
    range
    fixed_zero
    variable

A combined summary is written to:

    OUTPUT_DIR/fva_run_summary.tsv

Dependencies
------------
- Python >= 3.8
- numpy
- pandas
- cobra
- micom
- A supported optimization solver
"""

from __future__ import annotations

import argparse
import pickle
import platform
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cobra
import micom
import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Study design
#
# These sample assignments are part of the experimental design rather than
# machine-specific configuration and are intentionally preserved.
# ------------------------------------------------------------------------------

SAMPLES = {
    "CK": [
        "SRR16095329",
        "SRR16095330",
        "SRR16095331",
        "SRR16095332",
    ],
    "NP": [
        "SRR16095324",
        "SRR16095325",
        "SRR16095326",
        "SRR16095327",
    ],
    "NPM": [
        "SRR16095333",
        "SRR16095334",
        "SRR16095335",
        "SRR16095336",
    ],
}

VALID_ASSUMPTIONS = (
    "objective_only",
    "fix_community_growth",
    "fix_taxon_growth",
)

VALID_REACTION_MODES = (
    "active",
    "all",
    "internal",
)

VALID_SOLVERS = (
    "gurobi",
    "cplex",
    "glpk",
)

DEFAULT_SELECTION_COLUMN = "selected_by_threshold"
DEFAULT_TRADEOFF = 0.5
DEFAULT_FRACTION_OF_OPTIMUM = 1.0
DEFAULT_ASSUMPTION = "fix_taxon_growth"
DEFAULT_REACTION_MODE = "active"
DEFAULT_SOLVER = "gurobi"
DEFAULT_PROCESSES = 1
DEFAULT_TOLERANCE = 1e-6


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run flux variability analysis on MICOM community models while "
            "conditioning the feasible space on a cooperative-tradeoff "
            "reference solution."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--base_path",
        required=True,
        type=Path,
        help=(
            "Base directory containing condition/sample MICOM output "
            "directories and community.pickle files."
        ),
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="Directory in which FVA result files will be written.",
    )

    parser.add_argument(
        "--selection_summary",
        default=None,
        type=Path,
        help=(
            "Path to cross_check_reaction_summary.tsv used for reproducible "
            "active-reaction selection. Required when --reactions active."
        ),
    )

    parser.add_argument(
        "--selection_col",
        default=DEFAULT_SELECTION_COLUMN,
        help=(
            "Boolean column in --selection_summary indicating which "
            "reactions should be retained."
        ),
    )

    parser.add_argument(
        "--tradeoff",
        type=float,
        default=DEFAULT_TRADEOFF,
        help=(
            "MICOM cooperative-tradeoff fraction used to compute the "
            "reference solution."
        ),
    )

    parser.add_argument(
        "--fraction",
        type=float,
        default=DEFAULT_FRACTION_OF_OPTIMUM,
        help=(
            "COBRApy FVA fraction_of_optimum. Use 1.0 when fixing the "
            "tradeoff solution. Lower values broaden the feasible space."
        ),
    )

    parser.add_argument(
        "--assumption",
        choices=VALID_ASSUMPTIONS,
        default=DEFAULT_ASSUMPTION,
        help=(
            "How strictly the MICOM cooperative-tradeoff solution is "
            "preserved during FVA."
        ),
    )

    parser.add_argument(
        "--pfba_factor",
        type=float,
        default=None,
        help=(
            "Optional COBRApy pFBA factor. Values greater than 1 constrain "
            "total absolute flux near the parsimonious solution. Examples: "
            "1.05 or 1.1."
        ),
    )

    parser.add_argument(
        "--reactions",
        choices=VALID_REACTION_MODES,
        default=DEFAULT_REACTION_MODE,
        help=(
            "Reaction-selection mode: active uses the selection-summary "
            "file; all uses all community EX_*_m exchanges; internal uses "
            "all non-medium reactions."
        ),
    )

    parser.add_argument(
        "--solver",
        choices=VALID_SOLVERS,
        default=DEFAULT_SOLVER,
        help="Optimization solver used for MICOM and FVA.",
    )

    parser.add_argument(
        "--only_sample",
        default=None,
        help=(
            "Run only one condition/sample combination, for example "
            "CK_SRR16095329."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing per-sample FVA result files.",
    )

    parser.add_argument(
        "--loopless",
        action="store_true",
        help=(
            "Run loopless FVA. This can be extremely slow for large "
            "community models."
        ),
    )

    parser.add_argument(
        "--processes",
        type=int,
        default=DEFAULT_PROCESSES,
        help="Number of processes passed to COBRApy FVA.",
    )

    parser.add_argument(
        "--tol",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=(
            "Tolerance used for fixed growth bounds, zero-flux calls, "
            "and variable-range classification."
        ),
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Argument validation
# ------------------------------------------------------------------------------

def validate_arguments(args: argparse.Namespace) -> None:
    """Validate command-line arguments and resolve paths."""
    args.base_path = args.base_path.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    if args.selection_summary is not None:
        args.selection_summary = (
            args.selection_summary
            .expanduser()
            .resolve()
        )

    if not args.base_path.is_dir():
        raise FileNotFoundError(
            "MICOM base directory does not exist or is not a directory: "
            f"{args.base_path}"
        )

    if args.reactions == "active":
        if args.selection_summary is None:
            raise ValueError(
                "--selection_summary is required when --reactions active."
            )

        if not args.selection_summary.is_file():
            raise FileNotFoundError(
                "Reaction-selection summary does not exist or is not a file: "
                f"{args.selection_summary}"
            )

        if args.selection_summary.stat().st_size == 0:
            raise ValueError(
                "Reaction-selection summary is empty: "
                f"{args.selection_summary}"
            )

    if not 0 < args.tradeoff <= 1:
        raise ValueError(
            "--tradeoff must be greater than 0 and less than or equal to 1."
        )

    if not 0 < args.fraction <= 1:
        raise ValueError(
            "--fraction must be greater than 0 and less than or equal to 1."
        )

    if args.pfba_factor is not None and args.pfba_factor <= 1:
        raise ValueError(
            "--pfba_factor must be greater than 1 when provided."
        )

    if args.processes < 1:
        raise ValueError("--processes must be at least 1.")

    if args.tol <= 0:
        raise ValueError("--tol must be greater than 0.")

    if args.only_sample is not None:
        valid_sample_keys = {
            f"{condition}_{sample_id}"
            for condition, sample_ids in SAMPLES.items()
            for sample_id in sample_ids
        }

        if args.only_sample not in valid_sample_keys:
            raise ValueError(
                f"Unknown --only_sample value: {args.only_sample}. "
                "Expected one of: "
                + ", ".join(sorted(valid_sample_keys))
            )


# ------------------------------------------------------------------------------
# MICOM input loading
# ------------------------------------------------------------------------------

def load_community(pickle_path: Path):
    """Load a serialized MICOM Community object."""
    with pickle_path.open("rb") as handle:
        community = pickle.load(handle)

    return community


def apply_medium(
    community,
    medium_path: Path,
) -> None:
    """
    Read and apply a sample-specific MICOM medium.

    The medium file must contain:

        exchange
        uptake
    """
    medium_table = pd.read_csv(
        medium_path,
        sep="\t",
    )

    medium_table.columns = [
        str(column).strip()
        for column in medium_table.columns
    ]

    required_columns = {
        "exchange",
        "uptake",
    }

    if not required_columns.issubset(
        medium_table.columns
    ):
        raise ValueError(
            f"{medium_path} must contain columns: "
            f"{sorted(required_columns)}"
        )

    medium_table["exchange"] = (
        medium_table["exchange"]
        .astype(str)
        .str.strip()
    )

    medium_table["uptake"] = pd.to_numeric(
        medium_table["uptake"],
        errors="raise",
    )

    community.medium = dict(
        zip(
            medium_table["exchange"],
            medium_table["uptake"],
        )
    )


# ------------------------------------------------------------------------------
# Reaction selection
# ------------------------------------------------------------------------------

def load_selected_active_reactions(
    summary_path: Path,
    selection_col: str,
) -> List[str]:
    """
    Load selected exchange reactions from a reaction-summary table.

    The file is expected to contain reaction identifiers in the first column,
    which is read as the dataframe index. The selected reactions are indicated
    by a Boolean column such as selected_by_threshold.

    String values "true" and "false" are accepted and converted to Boolean.
    """
    if not summary_path.is_file():
        raise FileNotFoundError(
            f"Selection summary file not found: {summary_path}"
        )

    selection_table = pd.read_csv(
        summary_path,
        sep="\t",
        index_col=0,
    )

    if selection_col not in selection_table.columns:
        raise ValueError(
            f"Column '{selection_col}' was not found in {summary_path}. "
            f"Available columns: {list(selection_table.columns)}"
        )

    selected = selection_table[selection_col]

    if not pd.api.types.is_bool_dtype(selected):
        selected = (
            selected
            .astype(str)
            .str.strip()
            .str.lower()
            .map(
                {
                    "true": True,
                    "false": False,
                }
            )
        )

        if selected.isna().any():
            problematic_values = (
                selection_table.loc[
                    selected.isna(),
                    selection_col,
                ]
                .unique()
                .tolist()
            )

            raise ValueError(
                f"Could not parse Boolean values from column "
                f"'{selection_col}' in {summary_path}. "
                f"Problematic values: {problematic_values}"
            )

    reaction_ids = (
        selection_table
        .index[selected]
        .astype(str)
        .tolist()
    )

    if not reaction_ids:
        raise ValueError(
            f"No reactions were selected from {summary_path} using "
            f"column '{selection_col}'."
        )

    return reaction_ids


def get_target_reactions(
    community,
    mode: str = DEFAULT_REACTION_MODE,
    summary_path: Optional[Path] = None,
    selection_col: str = DEFAULT_SELECTION_COLUMN,
):
    """Return the community reactions selected for FVA."""
    if mode == "active":
        if summary_path is None:
            raise ValueError(
                "summary_path must be provided when mode='active'."
            )

        selected_ids = set(
            load_selected_active_reactions(
                summary_path,
                selection_col,
            )
        )

        target_reactions = [
            reaction
            for reaction in community.reactions
            if reaction.id in selected_ids
        ]

        found_ids = {
            reaction.id
            for reaction in target_reactions
        }

        missing_ids = (
            selected_ids
            - found_ids
        )

        if missing_ids:
            missing_preview = sorted(
                missing_ids
            )[:10]

            suffix = (
                " ..."
                if len(missing_ids) > 10
                else ""
            )

            print(
                "Warning: "
                f"{len(missing_ids)} selected reactions were not found "
                "in the community model: "
                f"{missing_preview}{suffix}"
            )

    elif mode == "all":
        target_reactions = [
            reaction
            for reaction in community.reactions
            if (
                reaction.id.startswith("EX_")
                and reaction.id.endswith("_m")
            )
        ]

    elif mode == "internal":
        target_reactions = [
            reaction
            for reaction in community.reactions
            if not (
                reaction.id.startswith("EX_")
                and reaction.id.endswith("_m")
            )
        ]

    else:
        raise ValueError(
            f"Unknown reaction mode: {mode}"
        )

    if not target_reactions:
        raise ValueError(
            f"No reactions found for mode={mode}."
        )

    return target_reactions


# ------------------------------------------------------------------------------
# Cooperative-tradeoff reference
# ------------------------------------------------------------------------------

def solve_tradeoff_reference(
    community,
    tradeoff: float,
    tolerance: float,
):
    """
    Solve MICOM cooperative tradeoff and return the reference solution.

    Returns
    -------
    solution
        MICOM cooperative-tradeoff solution.

    community_growth
        Reference community growth rate.

    taxon_growth
        Mapping from taxon identifier to reference taxon growth rate.
    """
    solution = community.cooperative_tradeoff(
        fraction=tradeoff,
        fluxes=False,
        pfba=False,
        atol=tolerance,
        rtol=tolerance,
    )

    community_growth = float(
        solution.growth_rate
    )

    if community_growth < tolerance:
        raise RuntimeError(
            "Community growth is near zero: "
            f"{community_growth}"
        )

    members = solution.members.copy()

    if "growth_rate" not in members.columns:
        raise RuntimeError(
            "MICOM solution.members has no 'growth_rate' column."
        )

    taxon_growth = (
        members["growth_rate"]
        .to_dict()
    )

    return (
        solution,
        community_growth,
        taxon_growth,
    )


# ------------------------------------------------------------------------------
# FVA conditioning
# ------------------------------------------------------------------------------

def set_community_growth_bounds(
    community,
    community_growth: float,
    tolerance: float,
) -> None:
    """
    Fix community growth to approximately the reference value.

    MICOM stores community growth in:

        community.variables.community_objective
    """
    community_objective = (
        community.variables.community_objective
    )

    lower_bound = max(
        0.0,
        community_growth - tolerance,
    )

    upper_bound = (
        community_growth + tolerance
    )

    community_objective.lb = lower_bound
    community_objective.ub = upper_bound


def set_taxon_growth_bounds(
    community,
    taxon_growth: Dict[str, float],
    tolerance: float,
) -> None:
    """
    Fix each taxon's growth rate to its reference tradeoff value.

    MICOM stores taxon growth constraints using names of the form:

        objective_<taxon>
    """
    missing_constraints = []

    for taxon in community.taxa:
        constraint_name = (
            f"objective_{taxon}"
        )

        if constraint_name not in community.constraints:
            missing_constraints.append(
                constraint_name
            )
            continue

        growth_rate = float(
            taxon_growth.get(
                taxon,
                0.0,
            )
        )

        lower_bound = max(
            0.0,
            growth_rate - tolerance,
        )

        upper_bound = (
            growth_rate + tolerance
        )

        constraint = (
            community.constraints[
                constraint_name
            ]
        )

        constraint.lb = lower_bound
        constraint.ub = upper_bound

    if missing_constraints:
        preview = missing_constraints[:10]

        suffix = (
            " ..."
            if len(missing_constraints) > 10
            else ""
        )

        raise RuntimeError(
            "Missing taxon growth constraints in MICOM model: "
            + ", ".join(preview)
            + suffix
        )


def restore_fba_objective(
    community,
) -> None:
    """
    Restore the MICOM community growth objective before running FVA.
    """
    community.objective = (
        community.scale
        * community.variables.community_objective
    )

    community.objective_direction = "max"


# ------------------------------------------------------------------------------
# FVA
# ------------------------------------------------------------------------------

def run_conditioned_fva(
    pickle_path: Path,
    medium_path: Path,
    solver: str,
    tradeoff: float,
    fraction: float,
    assumption: str,
    reaction_mode: str,
    summary_path: Optional[Path],
    selection_col: str,
    pfba_factor: Optional[float],
    loopless: bool,
    processes: int,
    tolerance: float,
) -> Tuple[pd.DataFrame, float, int]:
    """
    Run cooperative-tradeoff-conditioned FVA for one MICOM community.
    """
    community = load_community(
        pickle_path
    )

    community.solver = solver

    apply_medium(
        community,
        medium_path,
    )

    target_reactions = get_target_reactions(
        community,
        mode=reaction_mode,
        summary_path=summary_path,
        selection_col=selection_col,
    )

    (
        _,
        reference_growth,
        taxon_growth,
    ) = solve_tradeoff_reference(
        community,
        tradeoff,
        tolerance,
    )

    with community:
        restore_fba_objective(
            community
        )

        if assumption in {
            "fix_community_growth",
            "fix_taxon_growth",
        }:
            set_community_growth_bounds(
                community,
                reference_growth,
                tolerance,
            )

        if assumption == "fix_taxon_growth":
            set_taxon_growth_bounds(
                community,
                taxon_growth,
                tolerance,
            )

        fva_results = flux_variability_analysis(
            community,
            reaction_list=target_reactions,
            fraction_of_optimum=fraction,
            pfba_factor=pfba_factor,
            loopless=loopless,
            processes=processes,
        )

    fva_results = fva_results.copy()

    fva_results.index.name = "reaction"

    fva_results["range"] = (
        fva_results["maximum"]
        - fva_results["minimum"]
    )

    fva_results["fixed_zero"] = (
        (
            fva_results["maximum"].abs()
            <= tolerance
        )
        & (
            fva_results["minimum"].abs()
            <= tolerance
        )
    )

    fva_results["variable"] = (
        fva_results["range"].abs()
        > tolerance
    )

    return (
        fva_results,
        reference_growth,
        len(target_reactions),
    )


def summarize_fva(
    fva_results: pd.DataFrame,
    tolerance: float,
) -> Dict[str, int]:
    """Summarize the FVA results for one sample."""
    return {
        "n_reactions": int(
            len(fva_results)
        ),
        "n_variable": int(
            (
                fva_results["range"].abs()
                > tolerance
            ).sum()
        ),
        "n_fixed_zero": int(
            fva_results["fixed_zero"].sum()
        ),
        "n_nan": int(
            fva_results[
                ["minimum", "maximum"]
            ]
            .isna()
            .any(axis=1)
            .sum()
        ),
    }


# ------------------------------------------------------------------------------
# Summary-row helper
# ------------------------------------------------------------------------------

def make_summary_row(
    sample: str,
    condition: str,
    status: str,
    runtime_sec: float,
    reference_growth: float = np.nan,
    n_reactions: float = np.nan,
    n_variable: float = np.nan,
    n_fixed_zero: float = np.nan,
    n_nan: float = np.nan,
) -> dict:
    """Create one standardized FVA summary row."""
    return {
        "sample": sample,
        "condition": condition,
        "status": status,
        "mu_ref": reference_growth,
        "n_reactions": n_reactions,
        "n_variable": n_variable,
        "n_fixed_zero": n_fixed_zero,
        "n_nan": n_nan,
        "runtime_sec": runtime_sec,
    }


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Run MICOM-conditioned FVA for all selected samples."""
    total_start = time.perf_counter()

    args = parse_args()
    validate_arguments(args)

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print("=" * 72)
    print("MICOM-conditioned FVA")
    print("=" * 72)
    print(f"Python version     : {platform.python_version()}")
    print(f"MICOM version      : {micom.__version__}")
    print(f"COBRApy version    : {cobra.__version__}")
    print(f"pandas version     : {pd.__version__}")
    print(f"NumPy version      : {np.__version__}")
    print(f"base_path          : {args.base_path}")
    print(f"output_dir         : {args.output_dir}")
    print(f"assumption         : {args.assumption}")
    print(f"tradeoff           : {args.tradeoff}")
    print(f"fraction           : {args.fraction}")
    print(f"pfba_factor        : {args.pfba_factor}")
    print(f"reactions          : {args.reactions}")
    print(
        "selection_summary  : "
        f"{args.selection_summary or 'not used'}"
    )
    print(f"selection_col      : {args.selection_col}")
    print(f"solver             : {args.solver}")
    print(f"loopless           : {args.loopless}")
    print(f"processes          : {args.processes}")
    print(f"tol                : {args.tol}")
    print(
        "only_sample        : "
        f"{args.only_sample or 'all 12'}"
    )
    print("=" * 72)
    print()

    if args.reactions == "active":
        selected_ids = (
            load_selected_active_reactions(
                args.selection_summary,
                args.selection_col,
            )
        )

        print(
            f"Loaded {len(selected_ids)} selected active reactions "
            "from the summary file."
        )
        print()

    summary_rows = []

    for condition, sample_ids in SAMPLES.items():
        for sample_id in sample_ids:
            sample_key = (
                f"{condition}_{sample_id}"
            )

            if (
                args.only_sample
                and sample_key != args.only_sample
            ):
                continue

            pickle_path = (
                args.base_path
                / condition
                / sample_id
                / "community.pickle"
            )

            medium_path = (
                args.base_path
                / condition
                / sample_id
                / "medium_applied.tsv"
            )

            short_id = sample_id[-3:]

            output_csv = (
                args.output_dir
                / f"{condition}_{short_id}.csv"
            )

            print()
            print(f"[{sample_key}]")
            print(f"  pickle: {pickle_path}")
            print(f"  medium: {medium_path}")
            print(f"  output: {output_csv}")

            sample_start = (
                time.perf_counter()
            )

            if (
                output_csv.exists()
                and not args.overwrite
            ):
                runtime_sec = (
                    time.perf_counter()
                    - sample_start
                )

                print(
                    "  already exists, skipping"
                )

                print(
                    f"  runtime_sec   : "
                    f"{runtime_sec:.2f}"
                )

                summary_rows.append(
                    make_summary_row(
                        sample=sample_key,
                        condition=condition,
                        status="ALREADY_DONE",
                        runtime_sec=runtime_sec,
                    )
                )

                continue

            if not pickle_path.is_file():
                runtime_sec = (
                    time.perf_counter()
                    - sample_start
                )

                print(
                    f"  missing pickle: "
                    f"{pickle_path}"
                )

                print(
                    f"  runtime_sec   : "
                    f"{runtime_sec:.2f}"
                )

                summary_rows.append(
                    make_summary_row(
                        sample=sample_key,
                        condition=condition,
                        status="SKIP_NO_PICKLE",
                        runtime_sec=runtime_sec,
                    )
                )

                continue

            if not medium_path.is_file():
                runtime_sec = (
                    time.perf_counter()
                    - sample_start
                )

                print(
                    f"  missing medium file: "
                    f"{medium_path}"
                )

                print(
                    f"  runtime_sec   : "
                    f"{runtime_sec:.2f}"
                )

                summary_rows.append(
                    make_summary_row(
                        sample=sample_key,
                        condition=condition,
                        status="SKIP_NO_MEDIUM",
                        runtime_sec=runtime_sec,
                    )
                )

                continue

            try:
                (
                    fva_results,
                    reference_growth,
                    _,
                ) = run_conditioned_fva(
                    pickle_path=pickle_path,
                    medium_path=medium_path,
                    solver=args.solver,
                    tradeoff=args.tradeoff,
                    fraction=args.fraction,
                    assumption=args.assumption,
                    reaction_mode=args.reactions,
                    summary_path=args.selection_summary,
                    selection_col=args.selection_col,
                    pfba_factor=args.pfba_factor,
                    loopless=args.loopless,
                    processes=args.processes,
                    tolerance=args.tol,
                )

                fva_results.to_csv(
                    output_csv
                )

                statistics = summarize_fva(
                    fva_results,
                    args.tol,
                )

                runtime_sec = (
                    time.perf_counter()
                    - sample_start
                )

                print("  OK")
                print(
                    f"  mu_ref        : "
                    f"{reference_growth:.8f}"
                )
                print(
                    f"  n_reactions   : "
                    f"{statistics['n_reactions']}"
                )
                print(
                    f"  n_variable    : "
                    f"{statistics['n_variable']}"
                )
                print(
                    f"  n_fixed_zero  : "
                    f"{statistics['n_fixed_zero']}"
                )
                print(
                    f"  n_nan         : "
                    f"{statistics['n_nan']}"
                )
                print(
                    f"  runtime_sec   : "
                    f"{runtime_sec:.2f}"
                )

                summary_rows.append(
                    make_summary_row(
                        sample=sample_key,
                        condition=condition,
                        status="OK",
                        runtime_sec=runtime_sec,
                        reference_growth=reference_growth,
                        n_reactions=statistics[
                            "n_reactions"
                        ],
                        n_variable=statistics[
                            "n_variable"
                        ],
                        n_fixed_zero=statistics[
                            "n_fixed_zero"
                        ],
                        n_nan=statistics[
                            "n_nan"
                        ],
                    )
                )

            except Exception as error:
                runtime_sec = (
                    time.perf_counter()
                    - sample_start
                )

                print(
                    f"  FAIL: {error}"
                )

                print(
                    f"  runtime_sec   : "
                    f"{runtime_sec:.2f}"
                )

                traceback.print_exc()

                summary_rows.append(
                    make_summary_row(
                        sample=sample_key,
                        condition=condition,
                        status=f"FAIL: {error}",
                        runtime_sec=runtime_sec,
                    )
                )

    summary_table = pd.DataFrame(
        summary_rows,
        columns=[
            "sample",
            "condition",
            "status",
            "mu_ref",
            "n_reactions",
            "n_variable",
            "n_fixed_zero",
            "n_nan",
            "runtime_sec",
        ],
    )

    summary_path = (
        args.output_dir
        / "fva_run_summary.tsv"
    )

    summary_table.to_csv(
        summary_path,
        sep="\t",
        index=False,
    )

    total_runtime = (
        time.perf_counter()
        - total_start
    )

    successful_samples = int(
        (
            summary_table["status"]
            == "OK"
        ).sum()
    )

    print()
    print("=" * 72)
    print(
        "Done. Successful samples: "
        f"{successful_samples} / "
        f"{len(summary_table)}"
    )
    print(f"Summary: {summary_path}")
    print(
        f"Total runtime: "
        f"{total_runtime:.2f} sec"
    )

    if not summary_table.empty:
        print(
            summary_table.to_string(
                index=False
            )
        )

    print("=" * 72)


if __name__ == "__main__":
    main()