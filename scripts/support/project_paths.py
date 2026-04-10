#!/usr/bin/env python3
"""Shared project path helpers for cohort-specific reruns."""

from __future__ import annotations

import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"


def _resolve_env_path(env_name: str, default: Path) -> Path:
    raw = os.environ.get(env_name)
    if not raw:
        return default
    return Path(raw).expanduser().resolve()


RESULTS_DIR = _resolve_env_path("SCRNA_RESULTS_DIR", PROJECT_ROOT / "work" / "results")
FIGURES_DIR = _resolve_env_path("SCRNA_FIGURES_DIR", PROJECT_ROOT / "work" / "figures")
DOCS_EXPORTS_DIR = _resolve_env_path(
    "SCRNA_DOCS_EXPORTS_DIR",
    PROJECT_ROOT / "work" / "reports",
)
SAMPLE_MANIFEST_PATH = _resolve_env_path(
    "SCRNA_SAMPLE_MANIFEST",
    DATA_DIR / "manifests" / "manuscript_cohort.csv",
)
ANALYSIS_CONFIG_PATH = _resolve_env_path(
    "SCRNA_ANALYSIS_CONFIG",
    DATA_DIR / "analysis_config.json",
)


def ensure_dir(path: Path) -> Path:
    """Create a directory if needed and return the path."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def figure_dir(*parts: str) -> Path:
    """Return a figure subdirectory rooted at the active figures directory."""
    return ensure_dir(FIGURES_DIR.joinpath(*parts))


def docs_exports_dir() -> Path:
    """Return the active docs export directory."""
    return ensure_dir(DOCS_EXPORTS_DIR)

