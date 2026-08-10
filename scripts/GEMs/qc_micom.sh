#!/usr/bin/env python3
"""
Perform quality-control checks on MICOM community simulation outputs.

The script evaluates:

1. Model uniqueness
   - Calculates a SHA-256 checksum for every XML model.
   - Identifies byte-identical duplicate model files.

2. Community manifest integrity
   - Checks for required manifest files.
   - Counts taxa.
   - Verifies that relative abundances sum approximately to one.
   - Records the five most abundant taxa.
   - Generates a stable content fingerprint for each manifest.

3. Applied-medium consistency
   - Checks for required medium files.
   - Counts medium exchanges.
   - Generates a stable content fingerprint for each medium.
   - Flags suspicious exchange identifiers beginning with "EXs".

4. Member growth
   - Counts community members.
   - Counts growers with growth rates greater than 1e-6.
   - Counts members with growth rates less than or equal to 1e-9.
   - Reports minimum, median, and maximum growth rates.

5. Exchange-flux activity
   - Counts numeric flux entries.
   - Counts non-zero entries.
   - Calculates the fraction of non-zero entries.
   - Calculates total uptake and secretion across the flux table.

6. Output completeness
   - Checks for expected MICOM files and per-sample debug logs.

Usage
-----
The default directory structure is:

    PROJECT_DIR/
    ├── models/
    │   └── *.xml
    └── output/
        ├── CK/
        ├── NP/
        └── NPM/

Run with the default structure:

    python qc_micom.py \
        --project_dir /path/to/micom_project

Use explicit model and output directories:

    python qc_micom.py \
        --project_dir /path/to/micom_project \
        --models_dir /path/to/models \
        --output_dir /path/to/output

Outputs
-------
The following files are written under --output_dir:

    qc_model_hashes.tsv
        SHA-256 checksums for all model XML files.

    qc_summary.tsv
        Per-sample MICOM QC summary.

    qc_flags.tsv
        All detected QC issues.

Expected per-sample files
-------------------------
    community_manifest.tsv
    medium_applied.tsv
    member_growth_rates.csv
    exchange_fluxes.csv
    debug.log

Dependencies
------------
- Python >= 3.8
- pandas
"""

from __future__ import annotations

import argparse
import hashlib
import platform
import sys
from pathlib import Path
from typing import Any

import pandas as pd


# ------------------------------------------------------------------------------
# Python version requirement
# ------------------------------------------------------------------------------

if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required.")


# ------------------------------------------------------------------------------
# Reproducibility-critical constants
# ------------------------------------------------------------------------------

CONDITIONS = ("CK", "NP", "NPM")

ABUNDANCE_SUM_MIN = 0.999
ABUNDANCE_SUM_MAX = 1.001

ZERO_GROWTH_THRESHOLD = 1e-9
GROWER_THRESHOLD = 1e-6

HASH_CHUNK_SIZE = 1024 * 1024


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Perform integrity, growth, flux, medium, and model-uniqueness "
            "QC on MICOM simulation outputs."
        )
    )

    parser.add_argument(
        "--project_dir",
        required=True,
        type=Path,
        help=(
            "MICOM project directory. By default, models are read from "
            "<project_dir>/models and outputs from <project_dir>/output."
        ),
    )

    parser.add_argument(
        "--output_dir",
        default=None,
        type=Path,
        help=(
            "MICOM output directory. "
            "Default: <project_dir>/output."
        ),
    )

    parser.add_argument(
        "--models_dir",
        default=None,
        type=Path,
        help=(
            "Directory containing model XML files. "
            "Default: <project_dir>/models."
        ),
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------
# Hashing functions
# ------------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    """
    Calculate the SHA-256 checksum of a file.

    The file is read in 1 MiB chunks, preserving the original implementation.
    """
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(HASH_CHUNK_SIZE),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def fingerprint_df(dataframe: pd.DataFrame) -> str:
    """
    Calculate a stable SHA-256 fingerprint of dataframe content.

    The fingerprint procedure intentionally preserves the original method:

    1. Convert column names to strings.
    2. Sort columns alphabetically.
    3. Sort rows by all columns.
    4. Reset the row index.
    5. Serialize as CSV without the index.
    6. Hash the UTF-8 encoded CSV representation.
    """
    normalized = dataframe.copy()

    normalized.columns = [
        str(column)
        for column in normalized.columns
    ]

    normalized = normalized.sort_index(axis=1)

    normalized = normalized.sort_values(
        list(normalized.columns)
    ).reset_index(drop=True)

    serialized = normalized.to_csv(index=False)

    return hashlib.sha256(
        serialized.encode("utf-8")
    ).hexdigest()


# ------------------------------------------------------------------------------
# Flux QC
# ------------------------------------------------------------------------------

def flux_summary(fluxes_path: Path) -> dict[str, int | float]:
    """
    Summarize a MICOM exchange-flux table.

    MICOM solution flux tables commonly contain reactions as rows and community
    members as columns. The script removes a reaction/index column when named
    either "reaction" or "Unnamed: 0", then summarizes all remaining numeric
    entries.

    Sign convention:
        Negative exchange flux = uptake
        Positive exchange flux = secretion
    """
    fluxes = pd.read_csv(fluxes_path)

    for candidate in ("reaction", "Unnamed: 0"):
        if candidate in fluxes.columns:
            fluxes = fluxes.drop(columns=[candidate])
            break

    numeric_fluxes = fluxes.select_dtypes(
        include=["number"]
    )

    if numeric_fluxes.shape[1] == 0:
        return {
            "nonzero": 0,
            "total_entries": 0,
            "frac_nonzero": 0.0,
            "total_uptake": 0.0,
            "total_secretion": 0.0,
        }

    values = numeric_fluxes.to_numpy()

    total_entries = values.size
    nonzero = int((values != 0).sum())

    # Uptake is represented by negative exchange flux.
    total_uptake = float(
        values[values < 0].sum()
    )

    total_secretion = float(
        values[values > 0].sum()
    )

    return {
        "nonzero": nonzero,
        "total_entries": int(total_entries),
        "frac_nonzero": (
            nonzero / total_entries
            if total_entries
            else 0.0
        ),
        "total_uptake": total_uptake,
        "total_secretion": total_secretion,
    }


# ------------------------------------------------------------------------------
# Growth QC
# ------------------------------------------------------------------------------

def growth_summary(
    growth_path: Path,
) -> dict[str, int | float | None]:
    """
    Summarize MICOM member growth rates.

    The preferred column is "growth_rate". If it is absent, the final table
    column is used as a fallback, preserving the original behavior.
    """
    growth_table = pd.read_csv(growth_path)

    if "growth_rate" in growth_table.columns:
        growth_rates = pd.to_numeric(
            growth_table["growth_rate"],
            errors="coerce",
        ).dropna()
    else:
        growth_rates = pd.to_numeric(
            growth_table.iloc[:, -1],
            errors="coerce",
        ).dropna()

    if growth_rates.empty:
        return {
            "growers": 0,
            "n_members": 0,
            "gr_min": None,
            "gr_med": None,
            "gr_max": None,
            "n_zero": None,
        }

    n_members = int(len(growth_rates))

    n_zero = int(
        (
            growth_rates
            <= ZERO_GROWTH_THRESHOLD
        ).sum()
    )

    growers = int(
        (
            growth_rates
            > GROWER_THRESHOLD
        ).sum()
    )

    return {
        "growers": growers,
        "n_members": n_members,
        "gr_min": float(growth_rates.min()),
        "gr_med": float(growth_rates.median()),
        "gr_max": float(growth_rates.max()),
        "n_zero": n_zero,
    }


# ------------------------------------------------------------------------------
# Flag helper
# ------------------------------------------------------------------------------

def add_flag(
    flags: list[dict[str, str]],
    sample: str,
    condition: str,
    issue: str,
) -> None:
    """Append one standardized QC flag."""
    flags.append(
        {
            "sample": sample,
            "condition": condition,
            "issue": issue,
        }
    )


# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------

def resolve_directories(
    args: argparse.Namespace,
) -> tuple[Path, Path, Path]:
    """Resolve project, output, and model directories."""
    project_dir = args.project_dir.expanduser().resolve()

    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else project_dir / "output"
    )

    models_dir = (
        args.models_dir.expanduser().resolve()
        if args.models_dir is not None
        else project_dir / "models"
    )

    if not project_dir.is_dir():
        raise FileNotFoundError(
            "Project directory does not exist or is not a directory: "
            f"{project_dir}"
        )

    if not output_dir.is_dir():
        raise FileNotFoundError(
            "MICOM output directory does not exist or is not a directory: "
            f"{output_dir}"
        )

    if not models_dir.is_dir():
        raise FileNotFoundError(
            "Models directory does not exist or is not a directory: "
            f"{models_dir}"
        )

    return project_dir, output_dir, models_dir


# ------------------------------------------------------------------------------
# Model checksum QC
# ------------------------------------------------------------------------------

def inspect_model_hashes(
    models_dir: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, set[str], Path]:
    """
    Calculate model checksums and identify duplicate model files.
    """
    model_paths = sorted(
        models_dir.glob("*.xml")
    )

    if not model_paths:
        raise FileNotFoundError(
            f"No XML models found under: {models_dir}/*.xml"
        )

    model_hash_rows = []

    for model_path in model_paths:
        model_hash_rows.append(
            {
                "model_file": model_path.name,
                "sha256": sha256_file(model_path),
            }
        )

    model_hashes = pd.DataFrame(
        model_hash_rows,
        columns=["model_file", "sha256"],
    )

    model_hash_output = (
        output_dir
        / "qc_model_hashes.tsv"
    )

    model_hashes.to_csv(
        model_hash_output,
        sep="\t",
        index=False,
    )

    duplicate_counts = (
        model_hashes["sha256"]
        .value_counts()
    )

    duplicate_hashes = set(
        duplicate_counts[
            duplicate_counts > 1
        ].index
    )

    return (
        model_hashes,
        duplicate_hashes,
        model_hash_output,
    )


# ------------------------------------------------------------------------------
# Per-sample QC
# ------------------------------------------------------------------------------

def inspect_sample(
    sample_dir: Path,
    sample: str,
    condition: str,
    flags: list[dict[str, str]],
) -> dict[str, Any]:
    """Inspect all expected MICOM outputs for one sample."""
    manifest_path = (
        sample_dir
        / "community_manifest.tsv"
    )

    medium_path = (
        sample_dir
        / "medium_applied.tsv"
    )

    growth_path = (
        sample_dir
        / "member_growth_rates.csv"
    )

    fluxes_path = (
        sample_dir
        / "exchange_fluxes.csv"
    )

    debug_path = (
        sample_dir
        / "debug.log"
    )

    row: dict[str, Any] = {
        "sample": sample,
        "condition": condition,
        "sample_dir": str(sample_dir),
    }

    # --------------------------------------------------------------------------
    # Community manifest checks
    # --------------------------------------------------------------------------

    if manifest_path.exists():
        manifest = pd.read_csv(
            manifest_path,
            sep="\t",
        )

        row["n_taxa"] = int(
            len(manifest)
        )

        if "abundance" in manifest.columns:
            abundance_sum = float(
                pd.to_numeric(
                    manifest["abundance"],
                    errors="coerce",
                ).sum()
            )

            row["abundance_sum"] = (
                abundance_sum
            )

            if not (
                ABUNDANCE_SUM_MIN
                <= abundance_sum
                <= ABUNDANCE_SUM_MAX
            ):
                add_flag(
                    flags,
                    sample,
                    condition,
                    (
                        "Manifest abundance sum not ~1 "
                        f"(sum={abundance_sum:.6f})"
                    ),
                )
        else:
            add_flag(
                flags,
                sample,
                condition,
                "Manifest missing 'abundance' column",
            )

        if {"id", "abundance"}.issubset(
            manifest.columns
        ):
            top_taxa = (
                manifest
                .sort_values(
                    "abundance",
                    ascending=False,
                )
                .head(5)[["id", "abundance"]]
            )

            row["top5_taxa"] = ";".join(
                [
                    f"{record.id}:{record.abundance:.4g}"
                    for record in top_taxa.itertuples(
                        index=False
                    )
                ]
            )

        row["manifest_hash"] = (
            fingerprint_df(manifest)
        )

    else:
        add_flag(
            flags,
            sample,
            condition,
            "Missing community_manifest.tsv",
        )

    # --------------------------------------------------------------------------
    # Applied-medium checks
    # --------------------------------------------------------------------------

    if medium_path.exists():
        medium = pd.read_csv(
            medium_path,
            sep="\t",
        )

        row["n_medium_exchanges"] = int(
            len(medium)
        )

        row["medium_hash"] = (
            fingerprint_df(medium)
        )

        if "exchange" in medium.columns:
            suspicious = medium[
                medium["exchange"]
                .astype(str)
                .str.contains(
                    r"^EXs",
                    regex=True,
                )
            ]

            if len(suspicious) > 0:
                suspicious_ids = (
                    suspicious["exchange"]
                    .tolist()[:5]
                )

                add_flag(
                    flags,
                    sample,
                    condition,
                    (
                        "Suspicious medium exchange IDs: "
                        f"{suspicious_ids} "
                        "(check typos)"
                    ),
                )

    else:
        add_flag(
            flags,
            sample,
            condition,
            "Missing medium_applied.tsv",
        )

    # --------------------------------------------------------------------------
    # Member-growth checks
    # --------------------------------------------------------------------------

    if growth_path.exists():
        growth_stats = growth_summary(
            growth_path
        )

        row.update(growth_stats)

        if growth_stats["growers"] == 0:
            add_flag(
                flags,
                sample,
                condition,
                (
                    "No growers (>1e-6) in "
                    "member_growth_rates"
                ),
            )

        if (
            growth_stats["n_zero"] is not None
            and growth_stats["n_zero"]
            == growth_stats["n_members"]
        ):
            add_flag(
                flags,
                sample,
                condition,
                "All growth rates are ~0",
            )

    else:
        add_flag(
            flags,
            sample,
            condition,
            "Missing member_growth_rates.csv",
        )

    # --------------------------------------------------------------------------
    # Exchange-flux checks
    # --------------------------------------------------------------------------

    if fluxes_path.exists():
        flux_stats = flux_summary(
            fluxes_path
        )

        row.update(
            {
                f"flux_{key}": value
                for key, value
                in flux_stats.items()
            }
        )

        if flux_stats["total_entries"] == 0:
            add_flag(
                flags,
                sample,
                condition,
                (
                    "exchange_fluxes.csv "
                    "has no numeric entries"
                ),
            )

        if flux_stats["nonzero"] == 0:
            add_flag(
                flags,
                sample,
                condition,
                (
                    "All exchange flux "
                    "entries are zero"
                ),
            )

    else:
        add_flag(
            flags,
            sample,
            condition,
            "Missing exchange_fluxes.csv",
        )

    # --------------------------------------------------------------------------
    # Debug-log presence
    # --------------------------------------------------------------------------

    row["has_debug_log"] = (
        debug_path.exists()
    )

    return row


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

def main() -> None:
    """Run MICOM output quality control."""
    args = parse_args()

    (
        project_dir,
        output_dir,
        models_dir,
    ) = resolve_directories(args)

    print("=== MICOM quality control ===")
    print(f"Python version:   {platform.python_version()}")
    print(f"pandas version:   {pd.__version__}")
    print(f"Project directory: {project_dir}")
    print(f"Models directory:  {models_dir}")
    print(f"Output directory:  {output_dir}")
    print(f"Conditions:        {', '.join(CONDITIONS)}")
    print(
        "Abundance range:  "
        f"{ABUNDANCE_SUM_MIN}–{ABUNDANCE_SUM_MAX}"
    )
    print(
        "Zero threshold:   "
        f"<= {ZERO_GROWTH_THRESHOLD}"
    )
    print(
        "Grower threshold: "
        f"> {GROWER_THRESHOLD}"
    )
    print()

    # --------------------------------------------------------------------------
    # Model uniqueness via SHA-256
    # --------------------------------------------------------------------------

    (
        model_hashes,
        duplicate_hashes,
        model_hash_output,
    ) = inspect_model_hashes(
        models_dir=models_dir,
        output_dir=output_dir,
    )

    # --------------------------------------------------------------------------
    # Per-sample QC
    # --------------------------------------------------------------------------

    rows: list[dict[str, Any]] = []
    flags: list[dict[str, str]] = []

    for condition in CONDITIONS:
        condition_dir = (
            output_dir
            / condition
        )

        if not condition_dir.is_dir():
            add_flag(
                flags,
                sample="",
                condition=condition,
                issue=(
                    "Missing condition folder: "
                    f"{condition_dir}"
                ),
            )

            continue

        sample_dirs = sorted(
            path
            for path in condition_dir.iterdir()
            if path.is_dir()
        )

        for sample_dir in sample_dirs:
            sample = sample_dir.name

            row = inspect_sample(
                sample_dir=sample_dir,
                sample=sample,
                condition=condition,
                flags=flags,
            )

            rows.append(row)

    # --------------------------------------------------------------------------
    # Write per-sample QC summary
    # --------------------------------------------------------------------------

    if rows:
        qc_summary = (
            pd.DataFrame(rows)
            .sort_values(
                ["condition", "sample"]
            )
        )
    else:
        qc_summary = pd.DataFrame(
            columns=[
                "sample",
                "condition",
                "sample_dir",
            ]
        )

    qc_summary_output = (
        output_dir
        / "qc_summary.tsv"
    )

    qc_summary.to_csv(
        qc_summary_output,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Write flags table
    # --------------------------------------------------------------------------

    qc_flags = pd.DataFrame(
        flags,
        columns=[
            "sample",
            "condition",
            "issue",
        ],
    )

    qc_flags_output = (
        output_dir
        / "qc_flags.tsv"
    )

    qc_flags.to_csv(
        qc_flags_output,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------------------------
    # Human-readable summary
    # --------------------------------------------------------------------------

    print(f"[QC] Wrote: {qc_summary_output}")
    print(f"[QC] Wrote: {qc_flags_output}")
    print(f"[QC] Wrote: {model_hash_output}")

    if duplicate_hashes:
        duplicates = (
            model_hashes[
                model_hashes["sha256"]
                .isin(duplicate_hashes)
            ]
            .sort_values("sha256")
        )

        print(
            "\n[QC] WARNING: Duplicate model XMLs "
            "detected (byte-identical):"
        )

        print(
            duplicates.to_string(
                index=False
            )
        )

    else:
        print(
            "\n[QC] OK: No duplicate model XMLs "
            "detected by SHA256."
        )

    if not qc_flags.empty:
        print("\n[QC] Flags (first 20):")

        print(
            qc_flags
            .head(20)
            .to_string(index=False)
        )

    else:
        print("\n[QC] OK: No flags raised.")


if __name__ == "__main__":
    main()

