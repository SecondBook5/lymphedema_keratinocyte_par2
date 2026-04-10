#!/usr/bin/env python3
"""Regenerate a specific standalone panel from Figure G."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent


def parse_args(panel_required: bool = True, show_panel: bool = True) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate a specific standalone panel from Figure G."
    )
    parser.add_argument(
        "--panel",
        choices=["A", "B", "C", "D", "table"],
        required=panel_required,
        help="Standalone Figure G panel to render." if show_panel else argparse.SUPPRESS,
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
        help="Patient IDs to include. Default: 1 4 6",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Optional filename prefix for exported figure assets.",
    )
    return parser.parse_args()


def run_panel(panel: str, args: argparse.Namespace) -> None:
    annotated = (
        Path(args.annotated).expanduser().resolve()
        if args.annotated
        else (REPO_ROOT / "work" / args.cohort / "results" / "adata_annotated.h5ad").resolve()
    )

    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "keratinocyte_f2rl1_analysis.py"),
        "--adata",
        str(annotated),
        "--gene",
        args.gene,
        "--cell-type",
        args.cell_type,
        "--figure-set",
        "main_panel",
        "--panel-set",
        panel,
        "--figures-only",
    ]
    if args.prefix:
        cmd.extend(["--prefix", args.prefix])
    if args.patient_ids:
        cmd.extend(["--patient-ids", *[str(pid) for pid in args.patient_ids]])

    print(f"Regenerating Figure G panel {panel} from: {annotated}")
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    run_panel(args.panel, args)


if __name__ == "__main__":
    main()
