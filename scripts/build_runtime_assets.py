#!/usr/bin/env python3
"""Build paired_sc runtime assets for the keratinocyte paper repo."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFESTS_DIR = REPO_ROOT / "data" / "manifests"
DEFAULT_COHORT = "manuscript_cohort"


def pts_command() -> list[str]:
    exe_path = Path(sys.executable).resolve()
    candidate = exe_path.parent / "Scripts" / "paired-sc.exe"
    if candidate.exists():
        return [str(candidate)]
    return [sys.executable, "-m", "paired_sc.cli.app"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build runtime assets for a keratinocyte cohort rerun.")
    parser.add_argument("--cohort", default=DEFAULT_COHORT, help="Manifest stem under data/manifests.")
    parser.add_argument("--matrix-root", default=os.environ.get("SCRNA_MATRIX_ROOT"), help="Directory containing H5 files.")
    parser.add_argument("--workdir", default=None, help="Work root for runtime assets and package outputs.")
    parser.add_argument("--validate", action="store_true", help="Run `paired-sc validate` after writing the assets.")
    return parser.parse_args()


def default_matrix_root() -> Path | None:
    candidates = [
        os.environ.get("SCRNA_MATRIX_ROOT"),
        str(REPO_ROOT.parent / "Lymphadema_Mast_Cell_scrna-seq_Analysis" / "data" / "h5_files"),
    ]
    for raw in candidates:
        if not raw:
            continue
        path = Path(raw).expanduser().resolve()
        if path.exists():
            return path
    return None


def build_runtime_manifest(cohort: str, matrix_root: Path, runtime_dir: Path) -> Path:
    source = MANIFESTS_DIR / f"{cohort}.csv"
    if not source.exists():
        raise FileNotFoundError(f"Unknown cohort manifest: {source}")

    cohort_df = pd.read_csv(source).copy()
    runtime_df = pd.DataFrame(
        {
            "sample_id": cohort_df["sample_id"],
            "donor_id": cohort_df["donor_id"],
            "condition": cohort_df["condition"],
            "batch": cohort_df["batch"],
            "matrix_h5": cohort_df["matrix_filename"].map(lambda x: str((matrix_root / str(x)).resolve())),
            "sample_name": cohort_df["sample_name"],
            "replicate_group": cohort_df["replicate_group"],
            "metrics_csv": pd.NA,
            "metadata_json": pd.NA,
        }
    )

    missing = [path for path in runtime_df["matrix_h5"] if not Path(path).exists()]
    if missing:
        preview = "\n".join(missing[:5])
        raise FileNotFoundError(f"Missing matrix files under {matrix_root}:\n{preview}")

    manifest_path = runtime_dir / "manifest.csv"
    runtime_df.to_csv(manifest_path, index=False)
    return manifest_path


def build_project_yaml(runtime_dir: Path) -> Path:
    payload = {
        "project_name": "lymphedema_keratinocyte_par2",
        "project_subtitle": "Paired keratinocyte F2RL1 / PAR2 rerun",
        "condition_key": "condition",
        "case_condition": "LE",
        "control_condition": "Normal",
        "donor_key": "donor_id",
        "batch_key": "batch",
        "sample_key": "sample_id",
        "replicate_key": "replicate_group",
        "condition_colors": {
            "Normal": "#1B4F8A",
            "LE": "#C41E3A",
        },
        "qc": {
            "min_genes": 200,
            "min_counts": 500,
            "max_counts": None,
            "max_pct_mt": 20.0,
            "min_cells_per_gene": 10,
        },
        "preprocess": {
            "target_sum": 10000.0,
            "hvg_flavor": "seurat_v3",
            "n_top_genes": 3000,
            "scale_max_value": 10.0,
            "pca_n_comps": 50,
            "harmony_max_iter": 20,
            "neighbors_n_neighbors": 15,
            "neighbors_n_pcs": 30,
            "umap_min_dist": 0.3,
            "umap_spread": 1.0,
        },
        "annotation": {
            "backend": "celltypist",
            "model": "Adult_Human_Skin.pkl",
            "primary_leiden_resolution": 0.5,
            "leiden_resolutions": [0.3, 0.5, 0.8, 1.2],
            "cell_type_key": "cell_type",
            "label_map_path": str((REPO_ROOT / "scripts" / "skin_label_adapter.py").resolve()),
        },
        "domains": {
            "enabled": ["liana", "magic", "trajectory", "latent", "regulatory"],
        },
        "outputs": {
            "results_dir": "results",
            "figures_dir": "figures",
            "reports_dir": "reports",
            "logs_dir": "logs",
        },
    }
    project_path = runtime_dir / "project.yaml"
    project_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return project_path


def main() -> None:
    args = parse_args()
    matrix_root = Path(args.matrix_root).expanduser().resolve() if args.matrix_root else default_matrix_root()
    if matrix_root is None:
        raise ValueError("Provide --matrix-root or set SCRNA_MATRIX_ROOT to a directory containing the H5 files.")

    workdir = Path(args.workdir).expanduser().resolve() if args.workdir else (REPO_ROOT / "work" / args.cohort).resolve()
    runtime_dir = workdir / "runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = build_runtime_manifest(args.cohort, matrix_root, runtime_dir)
    project_path = build_project_yaml(runtime_dir)

    print(f"Built runtime project: {project_path}")
    print(f"Built runtime manifest: {manifest_path}")

    if args.validate:
        subprocess.run(
            pts_command()
            + [
                "validate",
                "--project",
                str(project_path),
                "--manifest",
                str(manifest_path),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()

