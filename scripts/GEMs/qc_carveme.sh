#!/usr/bin/env python3
"""
Perform basic structural and optimization QC on CarveMe SBML models.

For each XML model, the script reports:

- Whether the SBML file can be read
- Optimization status
- Objective value under the model's existing constraints
- Numbers of reactions, metabolites, genes, and boundary reactions
- Objective reaction identifiers
- Presence of a biomass- or growth-like reaction
- QC flags for:
    - Small models: fewer than 400 reactions
    - Missing biomass/growth-like reactions
    - No detectable growth: objective value <= 1e-6
    - Non-optimal optimization
    - SBML read failures

The QC thresholds and calculations reproduce the original analysis.

Usage
-----
python qc_carveme.py \
    --models-dir /path/to/models \
    --output-dir /path/to/qc

Outputs
-------
qc_models.tsv
    Per-model QC results.

qc_summary.tsv
    Dataset-level summary statistics.

qc_flagged.tsv
    Models meeting at least one QC flag.

Dependencies
------------
- pandas
- cobrapy
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import platform
import sys
from pathlib import Path
from typing import Any

import cobra
import pandas as pd
from cobra.io import read_sbml_model


SMALL_MODEL_REACTION_THRESHOLD = 400
GROWTH_THRESHOLD = 1e-6


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Perform structural and optimization QC on CarveMe SBML models."
        )
    )

    parser.add_argument(
        "--models-dir",
        required=True,
        type=Path,
        help="Directory containing CarveMe XML models.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which QC tables will be written.",
    )

    return parser.parse_args()


def safe_float(value: Any) -> float:
    """Convert a value to float, returning NaN if conversion fails."""
    try:
        return float(value)
    except Exception:
        return float("nan")


def get_objective_rxn_ids(model: cobra.Model) -> list[str]:
    """Return unique reaction identifiers contributing to the objective."""
    ids: list[str] = []

    try:
        for reaction in list(model.objective.keys()):
            if hasattr(reaction, "id"):
                ids.append(reaction.id)
    except Exception:
        pass

    if not ids:
        try:
            objective_string = str(model.objective).lower()

            if "growth" in objective_string:
                ids.append("Growth")
            elif "biomass" in objective_string:
                ids.append("biomass")
        except Exception:
            pass

    unique_ids: list[str] = []
    seen: set[str] = set()

    for reaction_id in ids:
        if reaction_id not in seen:
            seen.add(reaction_id)
            unique_ids.append(reaction_id)

    return unique_ids


def has_biomass_like_rxn(
    model: cobra.Model,
    objective_rxn_ids: list[str],
) -> bool:
    """Determine whether the model contains a biomass- or growth-like reaction."""
    for reaction_id in objective_rxn_ids:
        reaction_id_lower = reaction_id.lower()

        if (
            "growth" in reaction_id_lower
            or "biomass" in reaction_id_lower
        ):
            return True

    for reaction in model.reactions:
        reaction_id = reaction.id.lower()
        reaction_name = (reaction.name or "").lower()

        if (
            "biomass" in reaction_id
            or "growth" in reaction_id
            or "biomass" in reaction_name
            or "growth" in reaction_name
        ):
            return True

    return False


def count_exchanges(model: cobra.Model) -> int:
    """
    Count boundary reactions using the same fallback order as the original QC.
    """
    try:
        return len(model.boundary)
    except Exception:
        pass

    try:
        return len(model.exchanges)
    except Exception:
        pass

    return sum(
        1
        for reaction in model.reactions
        if reaction.id.upper().startswith("EX_")
    )


def make_read_failure_row(
    model_id: str,
    relative_path: str,
    error: Exception,
) -> dict[str, Any]:
    """Construct a QC row for a model that could not be read."""
    return {
        "model_id": model_id,
        "path": relative_path,
        "read_ok": False,
        "opt_status": "READ_FAIL",
        "growth": float("nan"),
        "n_rxn": float("nan"),
        "n_met": float("nan"),
        "n_genes": float("nan"),
        "n_ex": float("nan"),
        "objective_rxns": "",
        "biomass_like": False,
        "flag_small": True,
        "flag_no_biomass": True,
        "flag_no_growth": True,
        "flag_not_optimal": True,
        "error": str(error),
    }


def inspect_model(
    model_path: Path,
    models_dir: Path,
) -> dict[str, Any]:
    """Read, optimize, and summarize one SBML model."""
    model_id = model_path.stem
    relative_path = os.path.relpath(model_path, models_dir)

    try:
        model = read_sbml_model(str(model_path))
    except Exception as error:
        return make_read_failure_row(
            model_id=model_id,
            relative_path=relative_path,
            error=error,
        )

    n_rxn = len(model.reactions)
    n_met = len(model.metabolites)
    n_genes = len(model.genes)
    n_ex = count_exchanges(model)

    objective_rxns = get_objective_rxn_ids(model)
    biomass_like = has_biomass_like_rxn(model, objective_rxns)

    growth = float("nan")
    opt_status = "OPT_FAIL"
    opt_error = ""

    try:
        solution = model.optimize()
        opt_status = solution.status

        if opt_status == "optimal":
            growth = safe_float(solution.objective_value)

    except Exception as error:
        opt_status = "OPT_FAIL"
        opt_error = str(error)

    flag_small = n_rxn < SMALL_MODEL_REACTION_THRESHOLD
    flag_not_optimal = opt_status != "optimal"

    flag_no_growth = flag_not_optimal or (
        not math.isnan(growth)
        and growth <= GROWTH_THRESHOLD
    )

    flag_no_biomass = not biomass_like

    return {
        "model_id": model_id,
        "path": relative_path,
        "read_ok": True,
        "opt_status": opt_status,
        "growth": growth,
        "n_rxn": n_rxn,
        "n_met": n_met,
        "n_genes": n_genes,
        "n_ex": n_ex,
        "objective_rxns": ",".join(objective_rxns),
        "biomass_like": biomass_like,
        "flag_small": flag_small,
        "flag_no_biomass": flag_no_biomass,
        "flag_no_growth": flag_no_growth,
        "flag_not_optimal": flag_not_optimal,
        "error": opt_error,
    }


def main() -> None:
    """Run model QC and write output tables."""
    args = parse_args()

    models_dir = args.models_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not models_dir.is_dir():
        raise SystemExit(
            f"Models directory does not exist or is not a directory: "
            f"{models_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    output_models = output_dir / "qc_models.tsv"
    output_summary = output_dir / "qc_summary.tsv"
    output_flagged = output_dir / "qc_flagged.tsv"

    # This intentionally reproduces the original non-recursive *.xml search.
    model_paths = sorted(
        Path(path)
        for path in glob.glob(str(models_dir / "*.xml"))
    )

    if not model_paths:
        raise SystemExit(
            f"No XML models found under: {models_dir}/*.xml"
        )

    print("=== CarveMe model QC ===")
    print(f"Python version:   {platform.python_version()}")
    print(f"COBRApy version:  {cobra.__version__}")
    print(f"pandas version:   {pd.__version__}")
    print(f"Models directory: {models_dir}")
    print(f"Output directory: {output_dir}")
    print(
        "Small-model threshold: "
        f"< {SMALL_MODEL_REACTION_THRESHOLD} reactions"
    )
    print(f"Growth threshold: <= {GROWTH_THRESHOLD}")
    print()

    rows = [
        inspect_model(
            model_path=model_path,
            models_dir=models_dir,
        )
        for model_path in model_paths
    ]

    results = pd.DataFrame(rows)

    output_columns = [
        "model_id",
        "read_ok",
        "opt_status",
        "growth",
        "n_rxn",
        "n_met",
        "n_genes",
        "n_ex",
        "objective_rxns",
        "biomass_like",
        "flag_small",
        "flag_no_biomass",
        "flag_no_growth",
        "flag_not_optimal",
        "path",
        "error",
    ]

    results[output_columns].sort_values("model_id").to_csv(
        output_models,
        sep="\t",
        index=False,
    )

    readable_models = results[results["read_ok"]].copy()

    summary = pd.DataFrame(
        [
            {
                "n_models": len(results),
                "n_read_ok": int(results["read_ok"].sum()),
                "n_optimal": int(
                    (results["opt_status"] == "optimal").sum()
                ),
                "median_growth": (
                    readable_models["growth"].median()
                    if len(readable_models)
                    else float("nan")
                ),
                "min_growth": (
                    readable_models["growth"].min()
                    if len(readable_models)
                    else float("nan")
                ),
                "max_growth": (
                    readable_models["growth"].max()
                    if len(readable_models)
                    else float("nan")
                ),
                "median_rxn": (
                    readable_models["n_rxn"].median()
                    if len(readable_models)
                    else float("nan")
                ),
                "median_met": (
                    readable_models["n_met"].median()
                    if len(readable_models)
                    else float("nan")
                ),
                "median_genes": (
                    readable_models["n_genes"].median()
                    if len(readable_models)
                    else float("nan")
                ),
                "median_ex": (
                    readable_models["n_ex"].median()
                    if len(readable_models)
                    else float("nan")
                ),
                "n_small": int(results["flag_small"].sum()),
                "n_no_biomass": int(
                    results["flag_no_biomass"].sum()
                ),
                "n_no_growth": int(
                    results["flag_no_growth"].sum()
                ),
                "n_not_optimal": int(
                    results["flag_not_optimal"].sum()
                ),
                "n_read_fail": int(
                    (~results["read_ok"]).sum()
                ),
            }
        ]
    )

    summary.to_csv(
        output_summary,
        sep="\t",
        index=False,
    )

    flagged = results[
        (~results["read_ok"])
        | results["flag_small"]
        | results["flag_no_biomass"]
        | results["flag_no_growth"]
        | results["flag_not_optimal"]
    ].copy()

    flagged[output_columns].sort_values("model_id").to_csv(
        output_flagged,
        sep="\t",
        index=False,
    )

    print("=== CarveMe QC summary ===")
    print(f"Models scanned: {len(results)}")
    print(f"Read OK:        {int(results['read_ok'].sum())}")
    print(
        "Optimal:        "
        f"{int((results['opt_status'] == 'optimal').sum())}"
    )
    print(f"Flagged:        {len(flagged)}")
    print()
    print(f"Per-model table: {output_models}")
    print(f"Summary:         {output_summary}")
    print(f"Flagged models:  {output_flagged}")
    print()

    print(summary.round(4).to_string(index=False))

    if len(flagged) > 0:
        print("\nFlagged models:")

        display_columns = [
            "model_id",
            "opt_status",
            "growth",
            "n_rxn",
            "n_met",
            "n_genes",
            "flag_small",
            "flag_no_biomass",
            "flag_no_growth",
            "flag_not_optimal",
        ]

        print(
            flagged[display_columns]
            .sort_values("model_id")
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()