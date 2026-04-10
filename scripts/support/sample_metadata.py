#!/usr/bin/env python3
"""Local manifest-backed metadata helpers for the paper repo."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MANIFESTS_DIR = PROJECT_ROOT / "data" / "manifests"


def normalize_sample_name(sample_name: object) -> Optional[str]:
    """Collapse replicate-style sample labels onto a shared base name."""
    if pd.isna(sample_name):
        return None
    sample = str(sample_name).strip()
    if not sample:
        return None
    return re.sub(r"([_-])v\d+$", "", sample, flags=re.IGNORECASE)


def _patient_id_from_donor_id(donor_id: object) -> pd._libs.missing.NAType | int:
    if pd.isna(donor_id):
        return pd.NA
    match = re.search(r"(\d+)", str(donor_id))
    if not match:
        return pd.NA
    return int(match.group(1))


def load_patient_manifest() -> pd.DataFrame:
    """Load and combine the curated cohort manifests shipped with this repo."""
    manifest_paths = sorted(MANIFESTS_DIR.glob("*.csv"))
    if not manifest_paths:
        raise FileNotFoundError(f"No cohort manifests found under {MANIFESTS_DIR}")

    frames = [pd.read_csv(path) for path in manifest_paths]
    manifest = pd.concat(frames, ignore_index=True).drop_duplicates("sample_id").copy()
    manifest["sample_base"] = manifest["replicate_group"].map(normalize_sample_name)
    manifest["patient_id"] = manifest["donor_id"].map(_patient_id_from_donor_id).astype("Int64")
    manifest["patient_label"] = manifest["donor_id"].astype(str)

    cond_counts = (
        manifest.dropna(subset=["patient_id"])
        .groupby("patient_id")["condition"]
        .nunique()
        .rename("n_conditions")
        .reset_index()
    )
    manifest = manifest.merge(cond_counts, on="patient_id", how="left")
    manifest["paired_patient"] = manifest["n_conditions"].fillna(0).ge(2)
    return manifest


def attach_patient_metadata(obs: pd.DataFrame, sample_col: str = "sample_id") -> pd.DataFrame:
    """Append donor/patient metadata onto an obs dataframe."""
    if sample_col not in obs.columns:
        raise KeyError(f"Missing required sample column: {sample_col}")

    manifest = load_patient_manifest()
    keep_cols = [
        "sample_id",
        "sample_name",
        "sample_base",
        "batch",
        "condition",
        "donor_id",
        "patient_id",
        "patient_label",
        "paired_patient",
        "replicate_group",
    ]
    sample_meta = manifest[keep_cols].drop_duplicates("sample_id")

    out = obs.copy()
    overlap_cols = [c for c in sample_meta.columns if c in out.columns and c != sample_col]
    if overlap_cols:
        out = out.drop(columns=overlap_cols)

    out = out.merge(sample_meta, left_on=sample_col, right_on="sample_id", how="left")
    if sample_col != "sample_id":
        out = out.drop(columns=["sample_id"])
    return out
