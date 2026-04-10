#!/usr/bin/env python3
"""Study-level configuration for reusable two-condition scRNA-seq analyses."""

from __future__ import annotations

import json
from pathlib import Path

from support.project_paths import ANALYSIS_CONFIG_PATH


DEFAULT_CONFIG = {
    "project_title": "Keratinocyte F2RL1/PAR2 in Secondary Lymphedema",
    "project_subtitle": "Exact paired keratinocyte scRNA-seq analysis",
    "condition_key": "condition",
    "case_condition": "LE",
    "control_condition": "Normal",
    "condition_colors": {
        "LE": "#C41E3A",
        "Normal": "#1B4F8A",
    },
}


def _load_config_file(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _merge_config(defaults: dict, override: dict) -> dict:
    merged = dict(defaults)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            nested = dict(merged[key])
            nested.update(value)
            merged[key] = nested
        else:
            merged[key] = value
    return merged


def _slugify(label: str) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(label).strip())
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_") or "condition"


CONFIG_PATH = ANALYSIS_CONFIG_PATH
CONFIG = _merge_config(DEFAULT_CONFIG, _load_config_file(CONFIG_PATH))

PROJECT_TITLE = CONFIG["project_title"]
PROJECT_SUBTITLE = CONFIG["project_subtitle"]
CONDITION_KEY = CONFIG["condition_key"]
CASE_CONDITION = CONFIG["case_condition"]
CONTROL_CONDITION = CONFIG["control_condition"]
CONDITION_ORDER = [CONTROL_CONDITION, CASE_CONDITION]

CONDITION_COLORS = {
    CONTROL_CONDITION: CONFIG["condition_colors"].get(CONTROL_CONDITION, "#1B4F8A"),
    CASE_CONDITION: CONFIG["condition_colors"].get(CASE_CONDITION, "#C41E3A"),
}


def get_condition_order(values=None):
    """Return the configured condition order, optionally filtered to observed values."""
    if values is None:
        return CONDITION_ORDER.copy()
    observed = {str(v) for v in values if v is not None}
    ordered = [c for c in CONDITION_ORDER if c in observed]
    extras = sorted(observed - set(ordered))
    return ordered + extras


def contrast_label(joiner: str = "vs") -> str:
    """Human-readable condition contrast label."""
    return f"{CASE_CONDITION} {joiner} {CONTROL_CONDITION}"


def delta_label() -> str:
    """Human-readable delta label."""
    return f"{CASE_CONDITION} - {CONTROL_CONDITION}"


def contrast_stem() -> str:
    """Filesystem-safe condition contrast stem."""
    return f"{_slugify(CASE_CONDITION)}_vs_{_slugify(CONTROL_CONDITION)}"


