#!/usr/bin/env python3
"""Regenerate the keratinocyte state decomposition figure."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Regenerate the keratinocyte state decomposition figure.")
    parser.add_argument("--cohort", default="manuscript_cohort", help="Tracked cohort stem used to resolve the default annotated object.")
    parser.add_argument("--annotated", help="Annotated h5ad path. Defaults to work/<cohort>/results/adata_annotated.h5ad.")
    parser.add_argument("--gene", default="F2RL1", help="Gene symbol. Default: F2RL1")
    parser.add_argument("--cell-type", default="Keratinocytes", help="Annotated cell type. Default: Keratinocytes")
    parser.add_argument("--resolution", type=float, default=0.4, help="Leiden resolution for the state-level rerun. Default: 0.4")
    parser.add_argument("--prefix", default="", help="Optional filename prefix for exported figure assets.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    annotated = (
        Path(args.annotated).expanduser().resolve()
        if args.annotated
        else (REPO_ROOT / "work" / args.cohort / "results" / "adata_annotated.h5ad").resolve()
    )
    cmd = [
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
        "decomposition",
        "--figures-only",
    ]
    if args.prefix:
        cmd.extend(["--prefix", args.prefix])
    print("Regenerating state decomposition from:", annotated)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
