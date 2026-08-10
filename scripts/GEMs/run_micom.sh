#!/usr/bin/env python3
"""
Run sample-specific MICOM community simulations.

This script constructs one microbial community per sample from:

- A universal set of genome-scale metabolic models in SBML/XML format
- Sample-specific MAG relative abundances
- A sample-to-treatment mapping
- Treatment-specific MICOM medium files for CK, NP, and NPM

For each sample, the script:

1. Matches abundance-table MAG identifiers to available SBML models.
2. Optionally excludes models with fewer than a specified number of
   exchange reactions.
3. Removes zero-abundance taxa.
4. Applies a minimum relative-abundance threshold.
5. Renormalizes the remaining abundances to sum to one.
6. Constructs a MICOM Community.
7. Applies the treatment-specific medium.
8. Runs cooperative tradeoff optimization.
9. Writes member growth rates, exchange fluxes, the community object,
   debugging information, and a combined run summary.

Usage
-----
python run_micom.py \
    --models_dir /path/to/models \
    --abundance_long /path/to/abundance_micom.tsv \
    --sample_to_condition /path/to/sample_to_condition.tsv \
    --medium_ck /path/to/CK_micom_medium.tsv \
    --medium_np /path/to/NP_micom_medium.tsv \
    --medium_npm /path/to/NPM_micom_medium.tsv \
    --outdir /path/to/output

Optional example
----------------
python run_micom.py \
    --models_dir /path/to/models \
    --abundance_long /path/to/abundance_micom.tsv \
    --sample_to_condition /path/to/sample_to_condition.tsv \
    --medium_ck /path/to/CK_micom_medium.tsv \
    --medium_np /path/to/NP_micom_medium.tsv \
    --medium_npm /path/to/NPM_micom_medium.tsv \
    --outdir /path/to/output \
    --tradeoff 0.5 \
    --solver gurobi \
    --min_abundance 1e-4 \
    --min_exchanges 0 \
    --debug_exchanges

Input formats
-------------
Abundance table:
    sample    id    abundance

Sample-to-condition table:
    sample    condition

Condition values must be:
    CK
    NP
    NPM

MICOM medium tables:
    exchange    uptake

Outputs
-------
For each sample:

    OUTDIR/<condition>/<sample>/
        community_manifest.tsv
        medium_applied.tsv
        member_growth_rates.csv
        exchange_fluxes.csv
        community.pickle
        debug.log

When --debug_exchanges is used:

        debug_exchanges.txt

Dataset-level outputs:

    OUTDIR/micom_run_summary.tsv

A combined log is written by default to:

    MODELS_DIR/micom_combined.log

Dependencies
------------
- Python >= 3.8
- pandas
- micom
- cobra
"""

from __future__ import annotations

import argparse
import logging
import os
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import cobra
import micom
import pandas as pd
from micom import Community


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

VALID_CONDITIONS = ("CK", "NP", "NPM")

DEFAULT_TRADEOFF = 0.5
DEFAULT_SOLVER = "gurobi"
DEFAULT_MIN_ABUNDANCE = 1e-4
DEFAULT_MIN_EXCHANGES = 0
DEFAULT_COMBINED_LOG_NAME = "micom_combined.log"

GROWTH_THRESHOLD = 1e-6


# ------------------------------------------------------------------------------
# Logging configuration
# ------------------------------------------------------------------------------

logging.getLogger("cobra").setLevel(logging.ERROR)
logging.getLogger("micom").setLevel(logging.ERROR)


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run MICOM separately for each sample using universal models "
            "and treatment-specific MICOM medium files."
        )
    )

    parser.add_argument(
        "--models_dir",
        required=True,
        type=Path,
        help="Directory containing universal SBML models as *.xml files.",
    )

    parser.add_argument(
        "--abundance_long",
        required=True,
        type=Path,
        help=(
            "TSV or CSV abundance table with columns: "
            "sample, id, abundance."
        ),
    )

    parser.add_argument(
        "--sample_to_condition",
        required=True,
        type=Path,
        help=(
            "TSV file with columns: sample, condition. "
            "Condition must be CK, NP, or NPM."
        ),
    )

    parser.add_argument(
        "--only_sample",
        default=None,
        help=(
            "Run only one sample ID, for example SRR16095324. "
            "By default, all samples in the abundance table are run."
        ),
    )

    parser.add_argument(
        "--medium_ck",
        required=True,
        type=Path,
        help="CK MICOM-ready medium TSV with columns: exchange, uptake.",
    )

    parser.add_argument(
        "--medium_np",
        required=True,
        type=Path,
        help="NP MICOM-ready medium TSV with columns: exchange, uptake.",
    )

    parser.add_argument(
        "--medium_npm",
        required=True,
        type=Path,
        help="NPM MICOM-ready medium TSV with columns: exchange, uptake.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Base output directory.",
    )

    parser.add_argument(
        "--tradeoff",
        type=float,
        default=DEFAULT_TRADEOFF,
        help=(
            "Cooperative tradeoff fraction. "
            f"Default: {DEFAULT_TRADEOFF}."
        ),
    )

    parser.add_argument(
        "--solver",
        default=DEFAULT_SOLVER,
        help=(
            "Optimization solver used by MICOM. "
            f"Default: {DEFAULT_SOLVER}."
        ),
    )

    parser.add_argument(
        "--min_abundance",
        type=float,
        default=DEFAULT_MIN_ABUNDANCE,
        help=(
            "Minimum relative abundance retained before renormalization. "
            f"Default: {DEFAULT_MIN_ABUNDANCE}."
        ),
    )

    parser.add_argument(
        "--min_exchanges",
        type=int,
        default=DEFAULT_MIN_EXCHANGES,
        help=(
            "Exclude models with fewer than this number of EX_ reactions. "
            "Default: 0, meaning no exchange-count filtering. "
            "A value such as 150 may be used to exclude incomplete models."
        ),
    )

    parser.add_argument(
        "--dry_run",
        action="store_true",
        help=(
            "Build and write community manifests without constructing "
            "or optimizing MICOM communities."
        ),
    )

    parser.add_argument(
        "--debug_exchanges",
        action="store_true",
        help=(
            "Write previews of community exchange IDs and medium-match counts."
        ),
    )

    parser.add_argument(
        "--combined_log_name",
        default=DEFAULT_COMBINED_LOG_NAME,
        help=(
            "Combined log filename written under --models_dir. "
            f"Default: {DEFAULT_COMBINED_LOG_NAME}."
        ),
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def model_id_from_xml(xml_path: Path) -> str:
    """Return the model identifier derived from an XML filename."""
    return xml_path.stem


def validate_file(path: Path, description: str) -> None:
    """Validate that a required input file exists and is not empty."""
    if not path.is_file():
        raise FileNotFoundError(
            f"{description} does not exist or is not a file: {path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(f"{description} is empty: {path}")


def read_micom_medium_tsv(path: Path) -> Dict[str, float]:
    """
    Read a MICOM medium table and return exchange-to-uptake mapping.

    The table must contain the columns:

        exchange
        uptake

    Rows with non-positive uptake values are removed, reproducing the
    original behavior.
    """
    validate_file(path, "Medium TSV")

    medium_table = pd.read_csv(path, sep="\t")
    medium_table.columns = [
        column.strip()
        for column in medium_table.columns
    ]

    required_columns = {"exchange", "uptake"}

    if not required_columns.issubset(medium_table.columns):
        raise ValueError(
            "Medium TSV must contain columns: exchange, uptake "
            f"({path})"
        )

    medium_table["exchange"] = (
        medium_table["exchange"]
        .astype(str)
        .str.strip()
    )

    medium_table["uptake"] = medium_table["uptake"].astype(float)

    medium_table = medium_table[
        medium_table["uptake"] > 0
    ].copy()

    return dict(
        zip(
            medium_table["exchange"],
            medium_table["uptake"],
        )
    )


def count_exchanges_in_xml(xml_path: Path) -> int:
    """
    Count lines containing 'EX_' in an XML model file.

    This intentionally preserves the original exchange-counting method.
    It is a text-based count and does not parse the SBML model with COBRApy.
    """
    count = 0

    with xml_path.open("r", errors="replace") as handle:
        for line in handle:
            if "EX_" in line:
                count += 1

    return count


def build_model_index(
    models_dir: Path,
    min_exchanges: int = DEFAULT_MIN_EXCHANGES,
) -> Dict[str, str]:
    """
    Build a mapping from model ID to XML path.

    Models with fewer than min_exchanges text occurrences of 'EX_' are
    optionally excluded.
    """
    model_files = sorted(models_dir.glob("*.xml"))

    if not model_files:
        raise FileNotFoundError(
            f"No XML models found under: {models_dir}/*.xml"
        )

    model_index: Dict[str, str] = {}
    excluded: List[tuple[str, int]] = []

    for model_path in model_files:
        model_id = model_id_from_xml(model_path)

        if min_exchanges > 0:
            exchange_count = count_exchanges_in_xml(model_path)

            if exchange_count < min_exchanges:
                excluded.append((model_id, exchange_count))
                continue

        model_index[model_id] = str(model_path)

    if excluded:
        print(
            f"[MICOM] Excluded {len(excluded)} models with "
            f"<{min_exchanges} exchange reactions:",
            flush=True,
        )

        for model_id, exchange_count in sorted(
            excluded,
            key=lambda item: item[1],
        ):
            print(
                f"         {model_id}: "
                f"{exchange_count} exchanges",
                flush=True,
            )

    return model_index


def build_manifest_for_sample(
    sample_df: pd.DataFrame,
    id_to_file: Dict[str, str],
    min_abundance: float,
) -> pd.DataFrame:
    """
    Build a MICOM manifest for one sample.

    The original workflow is preserved:

    1. Keep taxa with available model files.
    2. Remove non-positive abundances.
    3. Apply the minimum abundance threshold.
    4. Renormalize retained abundances to sum to one.
    """
    present = sample_df[
        sample_df["id"].isin(id_to_file)
    ].copy()

    present["file"] = present["id"].map(id_to_file)

    present["abundance"] = present["abundance"].astype(float)

    present = present[
        present["abundance"] > 0
    ].copy()

    if min_abundance > 0:
        present = present[
            present["abundance"] >= min_abundance
        ].copy()

    if present.empty:
        return present

    present["abundance"] = (
        present["abundance"]
        / present["abundance"].sum()
    )

    return (
        present[["id", "file", "abundance"]]
        .reset_index(drop=True)
    )


def log_line(
    sample_log_path: Path,
    combined_lines: List[str],
    line: str,
    also_print: bool = False,
) -> None:
    """
    Append a timestamped message to the sample log and combined log buffer.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message = f"[{timestamp}] {line}"

    with sample_log_path.open("a") as handle:
        handle.write(message + "\n")

    combined_lines.append(message)

    if also_print:
        print(line, flush=True)


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate command-line inputs before starting the analysis."""
    args.models_dir = args.models_dir.expanduser().resolve()
    args.abundance_long = args.abundance_long.expanduser().resolve()
    args.sample_to_condition = (
        args.sample_to_condition.expanduser().resolve()
    )
    args.medium_ck = args.medium_ck.expanduser().resolve()
    args.medium_np = args.medium_np.expanduser().resolve()
    args.medium_npm = args.medium_npm.expanduser().resolve()
    args.outdir = args.outdir.expanduser().resolve()

    if not args.models_dir.is_dir():
        raise FileNotFoundError(
            "Models directory does not exist or is not a directory: "
            f"{args.models_dir}"
        )

    validate_file(
        args.abundance_long,
        "Abundance table",
    )

    validate_file(
        args.sample_to_condition,
        "Sample-to-condition table",
    )

    validate_file(
        args.medium_ck,
        "CK medium table",
    )

    validate_file(
        args.medium_np,
        "NP medium table",
    )

    validate_file(
        args.medium_npm,
        "NPM medium table",
    )

    if args.tradeoff <= 0 or args.tradeoff > 1:
        raise ValueError(
            "--tradeoff must be greater than 0 and less than or equal to 1."
        )

    if args.min_abundance < 0:
        raise ValueError("--min_abundance cannot be negative.")

    if args.min_exchanges < 0:
        raise ValueError("--min_exchanges cannot be negative.")

    if not args.combined_log_name.strip():
        raise ValueError("--combined_log_name cannot be empty.")


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Run MICOM community simulations."""
    args = parse_args()
    validate_arguments(args)

    args.outdir.mkdir(parents=True, exist_ok=True)

    combined_lines: List[str] = []

    combined_log_path = (
        args.models_dir
        / args.combined_log_name
    )

    print("=== MICOM community simulation ===")
    print(f"Python version:   {platform.python_version()}")
    print(f"MICOM version:    {micom.__version__}")
    print(f"COBRApy version:  {cobra.__version__}")
    print(f"pandas version:   {pd.__version__}")
    print(f"Models directory: {args.models_dir}")
    print(f"Output directory: {args.outdir}")
    print(f"Solver:           {args.solver}")
    print(f"Tradeoff:         {args.tradeoff}")
    print(f"Min abundance:    {args.min_abundance}")
    print(f"Min exchanges:    {args.min_exchanges}")
    print(f"Dry run:          {args.dry_run}")
    print()

    # --------------------------------------------------------------------------
    # Build model index once
    # --------------------------------------------------------------------------

    id_to_file = build_model_index(
        args.models_dir,
        min_exchanges=args.min_exchanges,
    )

    if not id_to_file:
        raise SystemExit(
            "No models remained after applying model filtering."
        )

    # --------------------------------------------------------------------------
    # Load abundance table
    # --------------------------------------------------------------------------

    abundance = pd.read_csv(
        args.abundance_long,
        sep=None,
        engine="python",
    )

    abundance.columns = [
        column.strip()
        for column in abundance.columns
    ]

    required_abundance_columns = {
        "sample",
        "id",
        "abundance",
    }

    if not required_abundance_columns.issubset(
        abundance.columns
    ):
        raise ValueError(
            "Abundance file must contain columns: "
            "sample, id, abundance"
        )

    abundance["sample"] = abundance["sample"].astype(str)
    abundance["id"] = abundance["id"].astype(str)
    abundance["abundance"] = (
        abundance["abundance"].astype(float)
    )

    # --------------------------------------------------------------------------
    # Load sample-to-condition mapping
    # --------------------------------------------------------------------------

    sample_conditions = pd.read_csv(
        args.sample_to_condition,
        sep="\t",
    )

    sample_conditions.columns = [
        column.strip()
        for column in sample_conditions.columns
    ]

    required_mapping_columns = {
        "sample",
        "condition",
    }

    if not required_mapping_columns.issubset(
        sample_conditions.columns
    ):
        raise ValueError(
            "Sample-to-condition table must contain columns: "
            "sample, condition"
        )

    sample_conditions["sample"] = (
        sample_conditions["sample"].astype(str)
    )

    sample_conditions["condition"] = (
        sample_conditions["condition"]
        .astype(str)
        .str.strip()
    )

    condition_map = dict(
        zip(
            sample_conditions["sample"],
            sample_conditions["condition"],
        )
    )

    # --------------------------------------------------------------------------
    # Load media once
    # --------------------------------------------------------------------------

    base_media = {
        "CK": read_micom_medium_tsv(args.medium_ck),
        "NP": read_micom_medium_tsv(args.medium_np),
        "NPM": read_micom_medium_tsv(args.medium_npm),
    }

    # --------------------------------------------------------------------------
    # Select samples
    # --------------------------------------------------------------------------

    samples = sorted(
        abundance["sample"].unique()
    )

    if args.only_sample is not None:
        samples = [str(args.only_sample)]

    print(
        "[MICOM] "
        f"samples={len(samples)} "
        f"models={len(id_to_file)} "
        f"min_abundance={args.min_abundance} "
        f"min_exchanges={args.min_exchanges}",
        flush=True,
    )

    results: List[dict] = []

    # --------------------------------------------------------------------------
    # Run each sample
    # --------------------------------------------------------------------------

    for sample in samples:
        condition = condition_map.get(sample)

        if condition not in VALID_CONDITIONS:
            message = (
                "SKIP: sample missing or invalid in "
                "sample_to_condition"
            )

            print(
                f"[{sample}] {message}",
                flush=True,
            )

            results.append(
                {
                    "sample": sample,
                    "status": "SKIP",
                    "n_taxa": 0,
                    "message": message,
                }
            )

            continue

        sample_abundance = abundance[
            abundance["sample"] == sample
        ][["id", "abundance"]].copy()

        manifest = build_manifest_for_sample(
            sample_abundance,
            id_to_file,
            args.min_abundance,
        )

        sample_out = (
            args.outdir
            / condition
            / sample
        )

        sample_out.mkdir(
            parents=True,
            exist_ok=True,
        )

        sample_log_path = (
            sample_out
            / "debug.log"
        )

        with sample_log_path.open("w") as handle:
            handle.write(
                f"# MICOM debug log for "
                f"{sample} ({condition})\n"
            )

        log_line(
            sample_log_path,
            combined_lines,
            (
                f"[{sample}] START "
                f"condition={condition} "
                f"solver={args.solver} "
                f"tradeoff={args.tradeoff} "
                f"min_exchanges={args.min_exchanges}"
            ),
            also_print=True,
        )

        if manifest.empty or len(manifest) < 2:
            message = (
                "SKIP: <2 taxa after filtering/matching"
            )

            log_line(
                sample_log_path,
                combined_lines,
                f"[{sample}] {message}",
                also_print=True,
            )

            results.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "status": "SKIP",
                    "n_taxa": len(manifest),
                    "message": message,
                }
            )

            continue

        manifest_path = (
            sample_out
            / "community_manifest.tsv"
        )

        manifest.to_csv(
            manifest_path,
            sep="\t",
            index=False,
        )

        if args.dry_run:
            message = (
                f"DRY_RUN: condition={condition} "
                f"taxa={len(manifest)}"
            )

            log_line(
                sample_log_path,
                combined_lines,
                f"[{sample}] {message}",
                also_print=True,
            )

            results.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "status": "DRY_RUN",
                    "n_taxa": len(manifest),
                }
            )

            continue

        try:
            log_line(
                sample_log_path,
                combined_lines,
                (
                    f"[{sample}] Building Community(...) "
                    f"taxa={len(manifest)}"
                ),
                also_print=True,
            )

            community = Community(
                manifest,
                solver=args.solver,
            )

            community.id = (
                f"{sample}_{condition}"
            )

            medium = base_media[condition]

            community_exchange_ids = [
                reaction.id
                for reaction in community.exchanges
            ]

            community_exchange_set = set(
                community_exchange_ids
            )

            if args.debug_exchanges:
                medium_keys = list(
                    medium.keys()
                )

                medium_matches = [
                    exchange_id
                    for exchange_id in medium_keys
                    if exchange_id in community_exchange_set
                ]

                exchange_debug_path = (
                    sample_out
                    / "debug_exchanges.txt"
                )

                with exchange_debug_path.open("w") as handle:
                    handle.write(
                        "First 30 community exchanges:\n"
                    )

                    handle.write(
                        str(community_exchange_ids[:30])
                        + "\n\n"
                    )

                    handle.write(
                        "First 30 medium exchanges:\n"
                    )

                    handle.write(
                        str(medium_keys[:30])
                        + "\n\n"
                    )

                    handle.write(
                        "Medium matches in community: "
                        f"{len(medium_matches)} / "
                        f"{len(medium_keys)}\n"
                    )

                    handle.write(
                        "Example matches (up to 50):\n"
                    )

                    handle.write(
                        str(medium_matches[:50])
                        + "\n"
                    )

                log_line(
                    sample_log_path,
                    combined_lines,
                    (
                        f"[{sample}] DEBUG medium matches "
                        f"{len(medium_matches)}/"
                        f"{len(medium_keys)} "
                        "(see debug_exchanges.txt)"
                    ),
                )

            matching_medium_ids = (
                set(medium.keys())
                .intersection(community_exchange_set)
            )

            if len(matching_medium_ids) == 0:
                exchange_preview = (
                    community_exchange_ids[:30]
                )

                raise RuntimeError(
                    "No ID from the MICOM-ready medium could be "
                    "found in the community exchange reactions.\n"
                    "Example community exchanges: "
                    f"{exchange_preview[:10]}"
                )

            medium_applied = pd.DataFrame(
                {
                    "exchange": list(medium.keys()),
                    "uptake": list(medium.values()),
                }
            )

            medium_applied.to_csv(
                sample_out / "medium_applied.tsv",
                sep="\t",
                index=False,
            )

            community.medium = medium

            log_line(
                sample_log_path,
                combined_lines,
                (
                    f"[{sample}] Solving "
                    "cooperative_tradeoff(...) "
                    "(this may take a while)"
                ),
                also_print=True,
            )

            solution = community.cooperative_tradeoff(
                fraction=args.tradeoff,
                fluxes=True,
                pfba=False,
            )

            log_line(
                sample_log_path,
                combined_lines,
                (
                    f"[{sample}] "
                    "cooperative_tradeoff finished"
                ),
                also_print=True,
            )

            solution.members.to_csv(
                sample_out
                / "member_growth_rates.csv"
            )

            solution.fluxes.to_csv(
                sample_out
                / "exchange_fluxes.csv"
            )

            community.to_pickle(
                sample_out
                / "community.pickle"
            )

            growth_rates = solution.members.get(
                "growth_rate",
                pd.Series(dtype=float),
            )

            growers = int(
                (
                    growth_rates
                    > GROWTH_THRESHOLD
                ).sum()
            )

            log_line(
                sample_log_path,
                combined_lines,
                f"[{sample}] OK growers={growers}",
                also_print=True,
            )

            results.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "status": "OK",
                    "n_taxa": len(manifest),
                    "growers": growers,
                }
            )

        except Exception as error:
            log_line(
                sample_log_path,
                combined_lines,
                (
                    f"[{sample}] FAIL "
                    f"{type(error).__name__}: "
                    f"{error}"
                ),
                also_print=True,
            )

            results.append(
                {
                    "sample": sample,
                    "condition": condition,
                    "status": "FAIL",
                    "n_taxa": len(manifest),
                    "message": str(error),
                }
            )

    # --------------------------------------------------------------------------
    # Write combined log
    # --------------------------------------------------------------------------

    with combined_log_path.open("w") as handle:
        handle.write("# MICOM combined log\n")

        for line in combined_lines:
            handle.write(line + "\n")

    print(
        f"[MICOM] Wrote combined log: "
        f"{combined_log_path}",
        flush=True,
    )

    # --------------------------------------------------------------------------
    # Write run summary
    # --------------------------------------------------------------------------

    summary = pd.DataFrame(results)

    summary_path = (
        args.outdir
        / "micom_run_summary.tsv"
    )

    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
    )

    print(
        f"[MICOM] Wrote run summary: "
        f"{summary_path}",
        flush=True,
    )


if __name__ == "__main__":
    main()
