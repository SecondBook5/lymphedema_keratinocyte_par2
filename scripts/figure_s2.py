#!/usr/bin/env python3
"""Regenerate the keratinocyte supplementary state-analysis figure set."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate the keratinocyte supplementary figure and state-analysis panels."
    )
    parser.add_argument(
        "--cohort",
        default="manuscript_cohort",
        help="Tracked cohort stem used to resolve the default annotated object.",
    )
    parser.add_argument(
        "--annotated",
        help="Annotated h5ad path. Defaults to work/<cohort>/results/adata_annotated.h5ad.",
    )
    parser.add_argument("--gene", default="F2RL1", help="Gene symbol. Default: F2RL1")
    parser.add_argument(
        "--cell-type",
        default="Keratinocytes",
        help="Annotated cell type. Default: Keratinocytes",
    )
    parser.add_argument(
        "--patient-ids",
        nargs="+",
        type=int,
        default=[1, 4, 6],
        help="Patient IDs to include for the expression summary. Default: 1 4 6",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.4,
        help="Leiden resolution for the state-level rerun. Default: 0.4",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Optional filename prefix for exported figure assets.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotated = (
        Path(args.annotated).expanduser().resolve()
        if args.annotated
        else (REPO_ROOT / "work" / args.cohort / "results" / "adata_annotated.h5ad").resolve()
    )

    supplement_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "keratinocyte_f2rl1_analysis.py"),
        "--adata",
        str(annotated),
        "--gene",
        args.gene,
        "--cell-type",
        args.cell_type,
        "--figure-set",
        "supplement",
        "--figures-only",
    ]
    if args.prefix:
        supplement_cmd.extend(["--prefix", args.prefix])
    if args.patient_ids:
        supplement_cmd.extend(["--patient-ids", *[str(pid) for pid in args.patient_ids]])

    state_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "keratinocyte_state_analysis.py"),
        "--adata",
        str(annotated),
        "--gene",
        args.gene,
        "--cell-type",
        args.cell_type,
        "--resolution",
        str(args.resolution),
        "--figure-set",
        "all",
        "--figures-only",
    ]
    if args.prefix:
        state_cmd.extend(["--prefix", args.prefix])

    print("Regenerating Figure S2 from:", annotated)
    subprocess.run(supplement_cmd, check=True)
    subprocess.run(state_cmd, check=True)


if __name__ == "__main__":
    main()
