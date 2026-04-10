#!/usr/bin/env python3
"""F2RL1/PAR2 keratinocyte analysis."""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import sparse, stats


THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parent
FIGURES_OUT = PROJECT_ROOT / "figures"
COMPOSITE_OUT = FIGURES_OUT / "composite"
PANELS_OUT = FIGURES_OUT / "panels"
TABLES_OUT = PROJECT_ROOT / "tables"
NOTES_OUT = PROJECT_ROOT / "notes"
DOCS_OUT = PROJECT_ROOT / "docs" / "generated"
DEFAULT_ADATA = PROJECT_ROOT / "work" / "manuscript_cohort" / "results" / "adata_annotated.h5ad"
DISPLAY_CONDITION_ORDER = ["Normal", "LE"]
DEFAULT_PREFIX = ""

FIGURES_OUT.mkdir(parents=True, exist_ok=True)
COMPOSITE_OUT.mkdir(parents=True, exist_ok=True)
PANELS_OUT.mkdir(parents=True, exist_ok=True)
TABLES_OUT.mkdir(parents=True, exist_ok=True)
NOTES_OUT.mkdir(parents=True, exist_ok=True)
DOCS_OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(THIS_DIR.resolve()))
from support.sample_metadata import attach_patient_metadata
from support import viz_config as viz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Advanced F2RL1/PAR2 keratinocyte analysis."
    )
    parser.add_argument(
        "--adata",
        default=str(DEFAULT_ADATA),
        help="Annotated AnnData file to analyze.",
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
        default=DEFAULT_PREFIX,
        help="Optional prefix for exported filenames inside this package.",
    )
    parser.add_argument(
        "--figures-only",
        action="store_true",
        help="Render figures only and skip tracked tables, notes, and report exports.",
    )
    parser.add_argument(
        "--figure-set",
        choices=["all", "main_panel", "supplement"],
        default="all",
        help="Subset of figure outputs to render. Default: all",
    )
    parser.add_argument(
        "--panel-set",
        choices=["all", "A", "B", "C", "D", "table", "E", "F", "G", "H"],
        default="all",
        help="Render a specific standalone panel only. Default: all",
    )
    return parser.parse_args()


def artifact_stem(name: str, prefix: str = "") -> str:
    prefix = prefix.strip().strip("_")
    return f"{prefix}_{name}" if prefix else name


def panel_selected(panel_set: str, panel_key: str) -> bool:
    return panel_set == "all" or panel_set == panel_key


def analysis_slug(gene: str, cell_type: str) -> str:
    return (
        f"{gene.lower()}_{cell_type.lower()}"
        .replace("/", "_")
        .replace(" ", "_")
        .replace("+", "plus")
        .replace("-", "_")
    )


def save_figure_exports(
    fig,
    alias_stem: str | None,
    export_stem: str,
    export_dir: Path,
    prefix: str,
    formats: tuple[str, ...] = ("pdf", "png"),
    tight: bool = False,
) -> None:
    targets: list[Path] = []
    if alias_stem:
        targets.append(FIGURES_OUT / artifact_stem(alias_stem, prefix))
    targets.append(export_dir / artifact_stem(export_stem, prefix))

    seen: set[str] = set()
    for target in targets:
        key = str(target)
        if key in seen:
            continue
        viz.save_figure(fig, target, formats=formats, tight=tight)
        seen.add(key)


def extract_vector(matrix_like) -> np.ndarray:
    if sparse.issparse(matrix_like):
        return matrix_like.toarray().ravel()
    return np.asarray(matrix_like).ravel()


def extract_counts(adata: sc.AnnData, gene: str) -> tuple[np.ndarray, str]:
    if "counts" in adata.layers and gene in adata.var_names:
        values = adata.layers["counts"][:, adata.var_names.get_loc(gene)]
        source = "layers[counts]"
    elif adata.raw is not None and gene in adata.raw.var_names:
        values = adata.raw[:, gene].X
        source = "raw.X"
    elif gene in adata.var_names:
        values = adata[:, gene].X
        source = "X"
    else:
        raise KeyError(f"Gene not found: {gene}")
    return extract_vector(values), source


def extract_normalized_expr(adata: sc.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene not found in adata.var_names: {gene}")
    return extract_vector(adata[:, gene].X)


def safe_paired_t(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return np.nan
    return float(stats.ttest_rel(x, y).pvalue)


def safe_wilcoxon(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.allclose(x - y, 0):
        return np.nan
    try:
        return float(
            stats.wilcoxon(
                x,
                y,
                alternative="two-sided",
                zero_method="wilcox",
                mode="exact",
            ).pvalue
        )
    except Exception:
        return np.nan


def safe_sign_test(deltas: np.ndarray) -> float:
    vals = np.asarray(deltas, dtype=float)
    vals = vals[~np.isclose(vals, 0)]
    if len(vals) == 0:
        return np.nan
    return float(stats.binomtest(int(np.sum(vals > 0)), len(vals), 0.5, alternative="two-sided").pvalue)


def short_cell_type_label(label: str) -> str:
    mapping = {
        "CD4+ T cells": "CD4 T cells",
        "CD8+ T cells": "CD8 T cells",
        "Regulatory T cells": "Treg cells",
        "Dendritic cells": "Dendritic",
        "Fibroblasts": "Fibroblasts",
        "Keratinocytes": "Keratinocytes",
        "Macrophages": "Macrophages",
        "Endothelial": "Endothelial",
        "Pericytes": "Pericytes",
        "Plasma cells": "Plasma cells",
        "Mast cells": "Mast cells",
        "Melanocytes": "Melanocytes",
    }
    return mapping.get(label, label)


def summarise_metrics(
    adata: sc.AnnData,
    gene: str,
    cell_type: str,
    patient_ids: list[int],
) -> dict[str, object]:
    obs_full = attach_patient_metadata(
        adata.obs[["sample_id", "condition", "cell_type"]].copy()
    )
    if obs_full["patient_id"].isna().any():
        missing = obs_full.loc[obs_full["patient_id"].isna(), "sample_id"].drop_duplicates().tolist()
        raise ValueError(f"Unmapped patient metadata for samples: {missing}")

    full_mask = obs_full["patient_id"].isin(patient_ids).to_numpy(dtype=bool)
    ct_mask = (
        obs_full["cell_type"].eq(cell_type)
        & obs_full["patient_id"].isin(patient_ids)
    ).to_numpy(dtype=bool)

    adata_ct = adata[ct_mask, :].copy()
    obs_ct = obs_full.loc[ct_mask].copy().reset_index(drop=True)
    if adata_ct.n_obs == 0:
        raise ValueError(f"No {cell_type} cells found for patients {patient_ids}")

    counts, count_source = extract_counts(adata_ct, gene)
    expr = extract_normalized_expr(adata_ct, gene)
    total_counts_per_cell = extract_vector(adata_ct.layers["counts"].sum(axis=1))

    obs_ct["gene_count"] = counts
    obs_ct["positive"] = counts > 0
    obs_ct["expr"] = expr
    obs_ct["cell_total_counts"] = total_counts_per_cell

    sample_summary = (
        obs_ct.groupby(
            ["patient_id", "patient_label", "sample_id", "sample_name", "sample_base", "condition"],
            observed=True,
        )
        .agg(
            total_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            mean_expr=("expr", "mean"),
            mean_expr_positive=("expr", lambda x: float(np.mean(x[x > 0])) if np.any(x > 0) else np.nan),
            gene_count_sum=("gene_count", "sum"),
            total_counts_sum=("cell_total_counts", "sum"),
        )
        .reset_index()
        .sort_values(["patient_id", "condition", "sample_id"])
    )
    sample_summary["pct_positive"] = sample_summary["positive_cells"] / sample_summary["total_cells"] * 100
    sample_summary["cpm"] = sample_summary["gene_count_sum"] / sample_summary["total_counts_sum"] * 1e6
    sample_summary["log1p_cpm"] = np.log1p(sample_summary["cpm"])

    patient_summary = (
        obs_ct.groupby(["patient_id", "patient_label", "condition"], observed=True)
        .agg(
            total_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            mean_expr=("expr", "mean"),
            mean_expr_positive=("expr", lambda x: float(np.mean(x[x > 0])) if np.any(x > 0) else np.nan),
            gene_count_sum=("gene_count", "sum"),
            total_counts_sum=("cell_total_counts", "sum"),
        )
        .reset_index()
        .sort_values(["patient_id", "condition"])
    )
    patient_summary["pct_positive"] = patient_summary["positive_cells"] / patient_summary["total_cells"] * 100
    patient_summary["cpm"] = patient_summary["gene_count_sum"] / patient_summary["total_counts_sum"] * 1e6
    patient_summary["log1p_cpm"] = np.log1p(patient_summary["cpm"])

    pooled_summary = (
        obs_ct.groupby("condition", observed=True)
        .agg(
            total_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            mean_expr=("expr", "mean"),
            gene_count_sum=("gene_count", "sum"),
            total_counts_sum=("cell_total_counts", "sum"),
        )
        .reset_index()
        .sort_values("condition")
    )
    pooled_summary["pct_positive"] = pooled_summary["positive_cells"] / pooled_summary["total_cells"] * 100
    pooled_summary["cpm"] = pooled_summary["gene_count_sum"] / pooled_summary["total_counts_sum"] * 1e6
    pooled_summary["log1p_cpm"] = np.log1p(pooled_summary["cpm"])

    pair_pct = (
        patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="pct_positive")
        .reset_index()
        .dropna(subset=["Normal", "LE"])
        .sort_values("patient_id")
    )
    pair_pct["delta_LE_minus_Normal_pp"] = pair_pct["LE"] - pair_pct["Normal"]

    pair_me = (
        patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="mean_expr")
        .reset_index()
        .dropna(subset=["Normal", "LE"])
        .sort_values("patient_id")
    )
    pair_me["delta_LE_minus_Normal"] = pair_me["LE"] - pair_me["Normal"]

    pair_pb = (
        patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="log1p_cpm")
        .reset_index()
        .dropna(subset=["Normal", "LE"])
        .sort_values("patient_id")
    )
    pair_pb["delta_LE_minus_Normal"] = pair_pb["LE"] - pair_pb["Normal"]

    pair_counts = (
        patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="positive_cells")
        .reset_index()
        .rename(columns={"Normal": "Normal_positive_cells", "LE": "LE_positive_cells"})
    )
    pair_totals = (
        patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="total_cells")
        .reset_index()
        .rename(columns={"Normal": "Normal_total_cells", "LE": "LE_total_cells"})
    )
    pair_pct = pair_pct.merge(pair_counts, on=["patient_id", "patient_label"], how="left")
    pair_pct = pair_pct.merge(pair_totals, on=["patient_id", "patient_label"], how="left")

    thresholds = [0, 1, 2, 3]
    sensitivity_rows = []
    for threshold in thresholds:
        tmp = obs_ct.copy()
        tmp["positive_threshold"] = tmp["gene_count"] > threshold
        tmp_patient = (
            tmp.groupby(["patient_id", "patient_label", "condition"], observed=True)["positive_threshold"]
            .agg(total_cells="size", positive_cells="sum")
            .reset_index()
        )
        tmp_patient["pct_positive"] = tmp_patient["positive_cells"] / tmp_patient["total_cells"] * 100
        tmp_wide = (
            tmp_patient.pivot(index=["patient_id", "patient_label"], columns="condition", values="pct_positive")
            .reset_index()
            .dropna(subset=["Normal", "LE"])
            .sort_values("patient_id")
        )
        tmp_wide["delta_LE_minus_Normal_pp"] = tmp_wide["LE"] - tmp_wide["Normal"]
        for _, row in tmp_wide.iterrows():
            sensitivity_rows.append(
                {
                    "threshold_gt": threshold,
                    "threshold_label": f"counts > {threshold}",
                    "patient_id": int(row["patient_id"]),
                    "patient_label": row["patient_label"],
                    "Normal_pct": float(row["Normal"]),
                    "LE_pct": float(row["LE"]),
                    "delta_LE_minus_Normal_pp": float(row["delta_LE_minus_Normal_pp"]),
                }
            )
        sensitivity_rows.append(
            {
                "threshold_gt": threshold,
                "threshold_label": f"counts > {threshold}",
                "patient_id": -1,
                "patient_label": "All",
                "Normal_pct": float(tmp_wide["Normal"].mean()),
                "LE_pct": float(tmp_wide["LE"].mean()),
                "delta_LE_minus_Normal_pp": float(tmp_wide["delta_LE_minus_Normal_pp"].mean()),
                "paired_ttest_p": safe_paired_t(tmp_wide["LE"].to_numpy(float), tmp_wide["Normal"].to_numpy(float)),
            }
        )
    sensitivity_df = pd.DataFrame(sensitivity_rows)

    le_pct = pair_pct["LE"].to_numpy(dtype=float)
    normal_pct = pair_pct["Normal"].to_numpy(dtype=float)
    le_asin = np.arcsin(np.sqrt(np.clip(le_pct / 100.0, 1e-8, 1 - 1e-8)))
    normal_asin = np.arcsin(np.sqrt(np.clip(normal_pct / 100.0, 1e-8, 1 - 1e-8)))

    fisher_or, fisher_p = stats.fisher_exact(
        [
            [
                int(pooled_summary.loc[pooled_summary["condition"].eq("LE"), "positive_cells"].iloc[0]),
                int(
                    pooled_summary.loc[pooled_summary["condition"].eq("LE"), "total_cells"].iloc[0]
                    - pooled_summary.loc[pooled_summary["condition"].eq("LE"), "positive_cells"].iloc[0]
                ),
            ],
            [
                int(pooled_summary.loc[pooled_summary["condition"].eq("Normal"), "positive_cells"].iloc[0]),
                int(
                    pooled_summary.loc[pooled_summary["condition"].eq("Normal"), "total_cells"].iloc[0]
                    - pooled_summary.loc[pooled_summary["condition"].eq("Normal"), "positive_cells"].iloc[0]
                ),
            ],
        ],
        alternative="two-sided",
    )

    stats_df = pd.DataFrame(
        [
            {
                "gene": gene,
                "cell_type": cell_type,
                "patient_ids": ",".join(str(pid) for pid in patient_ids),
                "n_paired_donors": int(len(pair_pct)),
                "positive_rule": "raw_counts_gt_zero",
                "count_source": count_source,
                "mean_pct_LE": float(pair_pct["LE"].mean()),
                "mean_pct_Normal": float(pair_pct["Normal"].mean()),
                "mean_delta_pct_pp": float(pair_pct["delta_LE_minus_Normal_pp"].mean()),
                "paired_ttest_pct_p": safe_paired_t(le_pct, normal_pct),
                "paired_ttest_pct_arcsin_p": safe_paired_t(le_asin, normal_asin),
                "paired_wilcoxon_pct_p": safe_wilcoxon(le_pct, normal_pct),
                "sign_test_pct_p": safe_sign_test(pair_pct["delta_LE_minus_Normal_pp"].to_numpy(float)),
                "paired_ttest_mean_expr_p": safe_paired_t(
                    pair_me["LE"].to_numpy(float), pair_me["Normal"].to_numpy(float)
                ),
                "paired_ttest_log1p_cpm_p": safe_paired_t(
                    pair_pb["LE"].to_numpy(float), pair_pb["Normal"].to_numpy(float)
                ),
                "pooled_fisher_odds_ratio": float(fisher_or),
                "pooled_fisher_p": float(fisher_p),
            }
        ]
    )

    return {
        "obs_full": obs_full,
        "obs_ct": obs_ct,
        "sample_summary": sample_summary,
        "patient_summary": patient_summary,
        "pooled_summary": pooled_summary,
        "pair_pct": pair_pct,
        "pair_me": pair_me,
        "pair_pb": pair_pb,
        "sensitivity_df": sensitivity_df,
        "stats_df": stats_df,
        "adata_ct": adata_ct,
        "full_mask": full_mask,
        "ct_mask": ct_mask,
    }


def plot_paired_metric(
    ax,
    pair_df: pd.DataFrame,
    metric_col_le: str,
    metric_col_normal: str,
    title: str,
    y_label: str,
    p_value: float,
    cond_colors: dict[str, str],
    annotate_counts: bool = False,
    annotate_donors: bool = True,
    xticklabels: list[str] | None = None,
    show_p_value_in_title: bool = True,
) -> None:
    cond_order = DISPLAY_CONDITION_ORDER
    all_y = []
    for _, row in pair_df.iterrows():
        ys = [float(row[metric_col_normal]), float(row[metric_col_le])]
        all_y.extend(ys)
        ax.plot([0, 1], ys, color="#7A7A7A", lw=0.9, zorder=1)
        ax.scatter(0, ys[0], s=40, color=cond_colors["Normal"], edgecolor="white", linewidth=0.6, zorder=2)
        ax.scatter(1, ys[1], s=40, color=cond_colors["LE"], edgecolor="white", linewidth=0.6, zorder=2)
        if annotate_donors:
            ax.text(1.10, ys[1], str(row["patient_label"]), fontsize=viz.FONTSIZE["annotation"], va="center")
        if annotate_counts:
            ax.annotate(
                f"{int(row['Normal_positive_cells'])}/{int(row['Normal_total_cells'])}",
                (0, ys[0]),
                xytext=(-7, 8),
                textcoords="offset points",
                fontsize=viz.FONTSIZE["annotation"] - 0.3,
                ha="right",
                va="bottom",
                color="#4A4A4A",
            )
            ax.annotate(
                f"{int(row['LE_positive_cells'])}/{int(row['LE_total_cells'])}",
                (1, ys[1]),
                xytext=(7, -8),
                textcoords="offset points",
                fontsize=viz.FONTSIZE["annotation"] - 0.3,
                ha="left",
                va="top",
                color="#4A4A4A",
            )
    mean_normal = float(pair_df[metric_col_normal].mean())
    mean_le = float(pair_df[metric_col_le].mean())
    ax.scatter(0, mean_normal, s=62, marker="D", color="#111111", edgecolor="white", linewidth=0.6, zorder=3)
    ax.scatter(1, mean_le, s=62, marker="D", color="#111111", edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(cond_order if xticklabels is None else xticklabels)
    ax.set_xlim(-0.38, 1.65)
    ax.set_ylabel(y_label)
    title_text = f"{title}\npaired t p={p_value:.3g}" if show_p_value_in_title else title
    ax.set_title(title_text, fontsize=viz.FONTSIZE["title"], pad=2)
    if all_y:
        ymin = float(np.min(all_y))
        ymax = float(np.max(all_y))
        span = max(ymax - ymin, 1e-6)
        ax.set_ylim(ymin - 0.10 * span, ymax + 0.16 * span)
    ax.grid(axis="y", color="#ECECEC", linewidth=0.6)
    sns.despine(ax=ax)


def plot_expression_violin(
    ax,
    obs_ct: pd.DataFrame,
    gene: str,
    cell_type: str,
    cond_colors: dict[str, str],
    seed: int = 42,
    show_inline_counts: bool = True,
    note_text: str | None = "Cell-level distributions shown\nDonor-level inference in panel D",
    xticklabels: list[str] | None = None,
) -> None:
    order = DISPLAY_CONDITION_ORDER
    plot_df = obs_ct[["condition", "expr"]].copy()

    sns.violinplot(
        data=plot_df,
        x="condition",
        y="expr",
        hue="condition",
        order=order,
        hue_order=order,
        palette=cond_colors,
        inner=None,
        linewidth=0.75,
        dodge=False,
        cut=0,
        saturation=0.95,
        ax=ax,
    )

    sns.boxplot(
        data=plot_df,
        x="condition",
        y="expr",
        hue="condition",
        order=order,
        hue_order=order,
        width=0.22,
        showfliers=False,
        dodge=False,
        palette={"Normal": "#FFFFFF", "LE": "#FFFFFF"},
        boxprops={"facecolor": "white", "alpha": 0.85, "linewidth": 0.75},
        whiskerprops={"linewidth": 0.75},
        medianprops={"linewidth": 0.9, "color": "#1A1A1A"},
        capprops={"linewidth": 0.75},
        ax=ax,
    )

    rng = np.random.default_rng(seed)
    point_rows = []
    for cond in order:
        cond_df = plot_df.loc[plot_df["condition"].eq(cond)].copy()
        n_points = min(450, len(cond_df))
        if n_points == 0:
            continue
        take = rng.choice(cond_df.index.to_numpy(), size=n_points, replace=False)
        point_rows.append(cond_df.loc[take])
    if point_rows:
        point_df = pd.concat(point_rows, ignore_index=True)
        sns.stripplot(
            data=point_df,
            x="condition",
            y="expr",
            order=order,
            color="black",
            size=1.6,
            alpha=0.35,
            jitter=0.23,
            ax=ax,
        )
    if ax.legend_ is not None:
        ax.legend_.remove()

    if show_inline_counts:
        for xpos, cond in enumerate(order):
            n_cells = int(plot_df["condition"].eq(cond).sum())
            ax.text(
                xpos,
                0.98,
                f"n={n_cells:,}",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=viz.FONTSIZE["annotation"],
                color="#444444",
            )

    ax.set_title(f"{gene} expression in {cell_type}", fontsize=viz.FONTSIZE["title"], pad=4)
    ax.set_xlabel("")
    if xticklabels is not None:
        ax.set_xticks(range(len(xticklabels)))
        ax.set_xticklabels(xticklabels)
    ax.set_ylabel("Normalized expression")
    if note_text:
        ax.text(
            0.02,
            0.04,
            note_text,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=viz.FONTSIZE["annotation"],
            color="#555555",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
        )
    sns.despine(ax=ax)


def plot_count_table(
    ax,
    pair_df: pd.DataFrame,
    cond_colors: dict[str, str],
    show_title: bool = True,
) -> None:
    display = pair_df.sort_values("patient_id").copy()
    rows = []
    for _, row in display.iterrows():
        rows.append(
            [
                str(row["patient_label"]),
                f"{int(row['Normal_positive_cells'])}/{int(row['Normal_total_cells'])}",
                f"{float(row['Normal']):.1f}%",
                f"{int(row['LE_positive_cells'])}/{int(row['LE_total_cells'])}",
                f"{float(row['LE']):.1f}%",
                f"{float(row['delta_LE_minus_Normal_pp']):.1f}",
            ]
        )
    rows.append(
        [
            "Mean",
            "",
            f"{float(display['Normal'].mean()):.1f}%",
            "",
            f"{float(display['LE'].mean()):.1f}%",
            f"{float(display['delta_LE_minus_Normal_pp'].mean()):.1f}",
        ]
    )

    ax.axis("off")
    if show_title:
        ax.set_title("Donor-level F2RL1+ keratinocyte summary", fontsize=viz.FONTSIZE["title"], pad=3)
    table = ax.table(
        cellText=rows,
        colLabels=["Donor", "Normal +/Total", "Normal %", "LE +/Total", "LE %", "Delta pp"],
        cellLoc="center",
        colLoc="center",
        bbox=[0.02, 0.05, 0.96, 0.84],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.1)
    table.scale(1.0, 1.04)

    n_rows = len(rows)
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        cell.set_edgecolor("#D9DDE2")
        if row == 0:
            cell.set_facecolor("#F2F4F7")
            cell.set_text_props(weight="bold", color="#1A1A1A")
        else:
            is_mean_row = row == n_rows
            is_alt = row % 2 == 0
            base_color = "#FFFFFF" if not is_alt else "#FBFCFD"
            if is_mean_row:
                base_color = "#F3F6FA"
            cell.set_facecolor(base_color)
            if col == 0:
                cell.set_text_props(weight="bold")
            if is_mean_row:
                cell.set_text_props(weight="bold")
    ax.text(
        0.02,
        -0.02,
        "Positive rule: raw counts > 0. Delta pp = LE - Normal percentage points.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color="#555555",
    )


def plot_patient_stratified_violin(
    ax,
    obs_ct: pd.DataFrame,
    gene: str,
    cell_type: str,
    cond_colors: dict[str, str],
) -> None:
    order = sorted(obs_ct["patient_label"].dropna().unique(), key=lambda x: int(str(x).replace("P", "")))
    violin = sns.violinplot(
        data=obs_ct,
        x="patient_label",
        y="expr",
        hue="condition",
        order=order,
        hue_order=DISPLAY_CONDITION_ORDER,
        palette=cond_colors,
        inner=None,
        cut=0,
        linewidth=0.7,
        saturation=0.95,
        dodge=True,
        ax=ax,
    )
    for collection in violin.collections:
        if hasattr(collection, "set_alpha"):
            collection.set_alpha(0.72)

    sns.boxplot(
        data=obs_ct,
        x="patient_label",
        y="expr",
        hue="condition",
        order=order,
        hue_order=DISPLAY_CONDITION_ORDER,
        palette={"Normal": "#FFFFFF", "LE": "#FFFFFF"},
        width=0.22,
        showfliers=False,
        dodge=True,
        boxprops={"facecolor": "white", "alpha": 0.88, "linewidth": 0.7},
        whiskerprops={"linewidth": 0.7},
        medianprops={"linewidth": 0.9, "color": "#1A1A1A"},
        capprops={"linewidth": 0.7},
        ax=ax,
    )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        dedup: dict[str, object] = {}
        for handle, label in zip(handles, labels):
            if label in {"Normal", "LE"} and label not in dedup:
                dedup[label] = handle
        ax.legend(
            dedup.values(),
            dedup.keys(),
            title="Condition",
            loc="upper right",
            frameon=False,
            fontsize=viz.FONTSIZE["legend"],
            title_fontsize=viz.FONTSIZE["legend"],
        )

    ax.set_title(
        f"{gene} expression in {cell_type}\nstratified by paired donor",
        fontsize=viz.FONTSIZE["title"],
        pad=4,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Normalized expression")
    ax.grid(axis="y", color="#F0F0F0", linewidth=0.55)
    sns.despine(ax=ax)


def build_main_panel(
    prefix: str,
    gene: str,
    cell_type: str,
    adata: sc.AnnData,
    summary: dict[str, object],
) -> None:
    obs_ct = summary["obs_ct"]
    pair_pct = summary["pair_pct"]
    stats_df = summary["stats_df"]
    ct_mask = summary["ct_mask"]
    full_mask = summary["full_mask"]
    obs_full = summary["obs_full"].loc[full_mask].copy().reset_index(drop=True)
    adata_full = adata[full_mask, :].copy()

    counts_full, _ = extract_counts(adata_full, gene)
    expr_full = extract_normalized_expr(adata_full, gene)
    obs_full["positive"] = counts_full > 0
    obs_full["expr"] = expr_full

    celltype_summary = (
        obs_full.groupby("cell_type", observed=True)
        .agg(
            total_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            mean_expr=("expr", "mean"),
        )
        .reset_index()
    )
    celltype_summary["pct_positive"] = celltype_summary["positive_cells"] / celltype_summary["total_cells"] * 100
    celltype_summary["short_label"] = celltype_summary["cell_type"].map(short_cell_type_label)
    celltype_summary = celltype_summary.sort_values(["pct_positive", "mean_expr"], ascending=[True, True]).reset_index(drop=True)
    celltype_summary["ypos"] = np.arange(len(celltype_summary))

    cond_colors = viz.get_condition_colors()
    keratinocyte_color = viz.get_cell_type_colors([cell_type])[cell_type]

    fig = plt.figure(figsize=(6.9, 5.25))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.06, 1.0],
        height_ratios=[1.02, 0.92],
        wspace=0.40,
        hspace=0.40,
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    fig.text(0.017, 0.978, "G", fontsize=13, fontweight="bold", va="top", ha="left")
    fig.suptitle(f"{gene} / PAR2 expression", fontsize=9.7, y=0.976, fontweight="bold")
    fig.add_artist(
        Line2D(
            [0.04, 0.985],
            [0.928, 0.928],
            transform=fig.transFigure,
            color="#222222",
            lw=0.6,
        )
    )

    if "X_umap" in adata_full.obsm:
        coords_full = adata_full.obsm["X_umap"]
        vmax = float(np.quantile(expr_full[expr_full > 0], 0.995)) if np.any(expr_full > 0) else 1.0
        axA.scatter(coords_full[:, 0], coords_full[:, 1], s=3.0, color="#DADADA", alpha=0.42, linewidths=0, rasterized=True)
        pos = expr_full > 0
        handle = axA.scatter(
            coords_full[pos, 0],
            coords_full[pos, 1],
            s=5.5,
            c=np.clip(expr_full[pos], 0, vmax),
            cmap="magma",
            vmin=0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        coords_ct = adata.obsm["X_umap"][ct_mask]
        axA.scatter(
            coords_ct[:, 0],
            coords_ct[:, 1],
            s=4.0,
            color=keratinocyte_color,
            alpha=0.14,
            linewidths=0,
            rasterized=True,
        )
        axA.set_title("All-cell UMAP", fontsize=8.6, pad=5)
        axA.set_xlabel("UMAP1")
        axA.set_ylabel("UMAP2")
        axA.set_xticks([])
        axA.set_yticks([])
        sns.despine(ax=axA, left=True, bottom=True)
        cax = inset_axes(axA, width="5.5%", height="38%", loc="center right", borderpad=0.8)
        cbar = fig.colorbar(handle, cax=cax)
        cbar.ax.tick_params(labelsize=viz.FONTSIZE["tick"] - 0.3)
        cbar.ax.set_title(gene, fontsize=viz.FONTSIZE["axis_title"] - 0.8, pad=3)

        inset = inset_axes(axA, width="34%", height="34%", loc="lower left", borderpad=0.55)
        expr_ct = obs_ct["expr"].to_numpy(dtype=float)
        inset.scatter(coords_ct[:, 0], coords_ct[:, 1], s=6, color="#E8E8E8", alpha=0.65, linewidths=0, rasterized=True)
        inset_pos = expr_ct > 0
        inset.scatter(
            coords_ct[inset_pos, 0],
            coords_ct[inset_pos, 1],
            s=8,
            c=np.clip(expr_ct[inset_pos], 0, vmax),
            cmap="magma",
            vmin=0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        inset.set_xticks([])
        inset.set_yticks([])
        for spine in inset.spines.values():
            spine.set_color("#444444")
            spine.set_linewidth(0.55)
    else:
        axA.axis("off")

    axB.hlines(
        y=celltype_summary["ypos"],
        xmin=0.08,
        xmax=celltype_summary["pct_positive"].clip(lower=0.08),
        color="#D1D1D1",
        linewidth=1.1,
        zorder=1,
    )
    point_colors = [keratinocyte_color if ct == cell_type else "#B7B7B7" for ct in celltype_summary["cell_type"]]
    point_edges = ["#1A1A1A" if ct == cell_type else "white" for ct in celltype_summary["cell_type"]]
    point_sizes = [62 if ct == cell_type else 36 for ct in celltype_summary["cell_type"]]
    axB.scatter(
        celltype_summary["pct_positive"].clip(lower=0.08),
        celltype_summary["ypos"],
        s=point_sizes,
        c=point_colors,
        edgecolors=point_edges,
        linewidths=0.6,
        zorder=2,
    )
    axB.set_xscale("log")
    axB.set_xlim(0.08, 60)
    axB.set_xticks([0.1, 0.3, 1, 3, 10, 30])
    axB.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"])
    axB.set_yticks(celltype_summary["ypos"])
    axB.set_yticklabels(celltype_summary["short_label"])
    axB.tick_params(axis="y", labelsize=6.2)
    axB.tick_params(axis="x", labelsize=6.6)
    for tick, ct in zip(axB.get_yticklabels(), celltype_summary["cell_type"]):
        if ct == cell_type:
            tick.set_fontweight("bold")
            tick.set_color("#1A1A1A")
    axB.set_title("Cell-type specificity", fontsize=8.6, pad=5)
    axB.set_xlabel(f"Percent {gene}+ cells (log)")
    axB.grid(axis="x", color="#EFEFEF", linewidth=0.6)
    sns.despine(ax=axB, left=False, bottom=False)

    plot_expression_violin(
        ax=axC,
        obs_ct=obs_ct,
        gene=gene,
        cell_type=cell_type,
        cond_colors=cond_colors,
        show_inline_counts=False,
        note_text=None,
        xticklabels=DISPLAY_CONDITION_ORDER,
    )
    axC.set_title(f"{gene} in keratinocytes", fontsize=8.6, pad=5)
    axC.grid(axis="y", color="#F1F1F1", linewidth=0.55)
    axC.tick_params(axis="x", pad=2)

    plot_paired_metric(
        ax=axD,
        pair_df=pair_pct,
        metric_col_le="LE",
        metric_col_normal="Normal",
        title="Donor-paired quantification",
        y_label="Percent positive",
        p_value=float(stats_df.at[0, "paired_ttest_pct_p"]),
        cond_colors=cond_colors,
        annotate_counts=False,
        annotate_donors=False,
        xticklabels=DISPLAY_CONDITION_ORDER,
    )
    axD.set_title(f"% {gene}+ keratinocytes", fontsize=8.6, pad=5)
    axD.tick_params(axis="x", pad=2)

    viz.add_panel_label(axA, "A", x=-0.14, y=1.05)
    viz.add_panel_label(axB, "B", x=-0.14, y=1.05)
    viz.add_panel_label(axC, "C", x=-0.14, y=1.05)
    viz.add_panel_label(axD, "D", x=-0.14, y=1.05)

    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.09, top=0.875, wspace=0.40, hspace=0.42)
    save_figure_exports(
        fig,
        alias_stem="main_panel",
        export_stem=f"{analysis_slug(gene, cell_type)}_main_panel",
        export_dir=COMPOSITE_OUT,
        prefix=prefix,
        tight=False,
    )
    plt.close(fig)


def build_main_caption_pdf(
    prefix: str,
    gene: str,
    cell_type: str,
    summary: dict[str, object],
) -> None:
    pair_pct = summary["pair_pct"].sort_values("patient_id").copy()
    stats_df = summary["stats_df"]
    sample_summary = summary["sample_summary"]
    has_duplicate_runs = (
        sample_summary.groupby(["patient_id", "condition"]).size().max() > 1
        if len(sample_summary) > 0
        else False
    )
    duplicate_note = (
        f"Donor-aware quantification defined {gene}-positive keratinocytes as cells with raw {gene} counts > 0 after collapsing duplicate runs to the patient-condition level. "
        if has_duplicate_runs
        else f"Donor-aware quantification defined {gene}-positive keratinocytes as cells with raw {gene} counts > 0 using only the specified sample set. "
    )

    caption_lines = [
        f"Figure G. {gene} / PAR2 expression in keratinocytes from the paired scRNA-seq cohort used for the keratinocyte/PAR2 manuscript analysis.",
        "Panels show the all-cell UMAP with keratinocyte inset, cross-cell-type specificity of percent F2RL1-positive cells, keratinocyte expression violins, and donor-paired quantification.",
        "Normal is shown on the left and LE is shown on the right throughout.",
        (
            duplicate_note
            + 
            f"Across paired donors, the mean percentage of {gene}-positive keratinocytes was {stats_df.at[0, 'mean_pct_Normal']:.1f}% in Normal skin and "
            f"{stats_df.at[0, 'mean_pct_LE']:.1f}% in LE skin "
            f"(paired t test p = {stats_df.at[0, 'paired_ttest_pct_p']:.4f}; arcsin-square-root transformed paired t test p = {stats_df.at[0, 'paired_ttest_pct_arcsin_p']:.4f}). "
            "Donor-level counts are listed below."
        ),
    ]

    fig = plt.figure(figsize=(6.85, 3.15))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.15], hspace=0.02)
    ax_text = fig.add_subplot(gs[0, 0])
    ax_table = fig.add_subplot(gs[1, 0])

    ax_text.axis("off")
    y = 0.98
    for line in caption_lines:
        wrapped = textwrap.fill(line, width=125)
        ax_text.text(0.0, y, wrapped, ha="left", va="top", fontsize=6.5, color="#222222")
        y -= 0.15 * (wrapped.count("\n") + 1) + 0.035

    plot_count_table(ax_table, pair_pct, viz.get_condition_colors(), show_title=False)
    viz.save_figure(
        fig,
        NOTES_OUT / artifact_stem("main_panel_caption", prefix),
        formats=("pdf",),
        tight=False,
    )
    plt.close(fig)


def build_supplement_caption_pdf(
    prefix: str,
    gene: str,
    cell_type: str,
    summary: dict[str, object],
) -> None:
    stats_df = summary["stats_df"]
    sensitivity_df = summary["sensitivity_df"]
    pvals = (
        sensitivity_df[sensitivity_df["patient_label"].eq("All")]
        .sort_values("threshold_gt")[["threshold_gt", "paired_ttest_p"]]
        .to_records(index=False)
    )
    threshold_text = "; ".join(
        f"counts > {int(thr)}: p={float(p):.4f}" for thr, p in pvals
    )

    lines = [
        f"Supplementary figure. E, {gene} expression in {cell_type.lower()} stratified by paired donor and condition.",
        f"F, donor-paired mean normalized {gene} expression across all {cell_type.lower()} (paired t p={stats_df.at[0, 'paired_ttest_mean_expr_p']:.4f}).",
        f"G, donor-paired {cell_type.lower()} pseudobulk {gene} abundance shown as log1p CPM after summing raw counts within donor-condition and scaling to CPM (paired t p={stats_df.at[0, 'paired_ttest_log1p_cpm_p']:.4f}).",
        f"H, threshold sensitivity of the LE-Normal difference for defining {gene}-positive cells: {threshold_text}.",
    ]

    fig = plt.figure(figsize=(6.85, 1.85))
    ax = fig.add_subplot(111)
    ax.axis("off")
    y = 0.98
    for i, line in enumerate(lines):
        wrapped = textwrap.fill(line, width=132)
        ax.text(0.0, y, wrapped, ha="left", va="top", fontsize=6.3 if i == 0 else 6.1, color="#222222")
        y -= 0.20 * (wrapped.count("\n") + 1)
    viz.save_figure(
        fig,
        NOTES_OUT / artifact_stem("supplement_caption", prefix),
        formats=("pdf",),
        tight=False,
    )
    plt.close(fig)


def build_main_figure(
    prefix: str,
    gene: str,
    cell_type: str,
    adata: sc.AnnData,
    summary: dict[str, object],
) -> None:
    obs_ct = summary["obs_ct"]
    pair_pct = summary["pair_pct"]
    stats_df = summary["stats_df"]
    adata_ct = summary["adata_ct"]
    ct_mask = summary["ct_mask"]

    fig = plt.figure(figsize=(6.85, 5.45))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0], wspace=0.28, hspace=0.26)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    cond_colors = viz.get_condition_colors()

    if "X_umap" in adata.obsm:
        coords_all = adata.obsm["X_umap"]
        coords_ct = coords_all[ct_mask]
        positive_mask = obs_ct["positive"].to_numpy(dtype=bool)

        axA.scatter(coords_all[:, 0], coords_all[:, 1], s=3, color="#D8D8D8", alpha=0.35, linewidths=0, rasterized=True)
        axA.scatter(coords_ct[:, 0], coords_ct[:, 1], s=5, color=viz.get_cell_type_colors([cell_type])[cell_type], alpha=0.65, linewidths=0, rasterized=True)
        axA.scatter(coords_ct[positive_mask, 0], coords_ct[positive_mask, 1], s=6, color="#CC2442", alpha=0.75, linewidths=0, rasterized=True)
        centroid = np.median(coords_ct, axis=0)
        axA.annotate(
            cell_type,
            xy=(centroid[0], centroid[1]),
            xytext=(centroid[0] + 2.2, centroid[1] + 1.4),
            arrowprops=dict(arrowstyle="->", lw=0.8, color="#444444"),
            fontsize=viz.FONTSIZE["annotation"] + 0.2,
            color="#303030",
        )
        axA.set_title("All cells\nkeratinocytes highlighted", fontsize=viz.FONTSIZE["title"], pad=2)
        axA.set_xlabel("UMAP1")
        axA.set_ylabel("UMAP2")
        axA.set_xticks([])
        axA.set_yticks([])
        sns.despine(ax=axA, left=True, bottom=True)
        legend_handles = [
            Line2D([0], [0], marker="o", linestyle="", markersize=4, markerfacecolor="#D8D8D8", markeredgecolor="none", label="All cells"),
            Line2D([0], [0], marker="o", linestyle="", markersize=4, markerfacecolor=viz.get_cell_type_colors([cell_type])[cell_type], markeredgecolor="none", label=cell_type),
            Line2D([0], [0], marker="o", linestyle="", markersize=4, markerfacecolor="#CC2442", markeredgecolor="none", label=f"{gene}+ {cell_type}"),
        ]
        axA.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize=viz.FONTSIZE["annotation"])
    else:
        axA.axis("off")
        axA.text(0.5, 0.5, "UMAP not available", ha="center", va="center")

    if "X_umap" in adata_ct.obsm:
        coords_ct = adata_ct.obsm["X_umap"]
        expr = obs_ct["expr"].to_numpy(dtype=float)
        axB.scatter(coords_ct[:, 0], coords_ct[:, 1], s=4, color="#DDDDDD", alpha=0.45, linewidths=0, rasterized=True)
        nz = expr > 0
        handle = axB.scatter(
            coords_ct[nz, 0],
            coords_ct[nz, 1],
            s=6,
            c=expr[nz],
            cmap="magma",
            linewidths=0,
            rasterized=True,
        )
        cb = fig.colorbar(handle, ax=axB, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=viz.FONTSIZE["tick"] - 0.3)
        cb.ax.set_title(gene, fontsize=viz.FONTSIZE["axis_title"] - 0.7, pad=3)
        axB.set_title(f"{cell_type}-only UMAP\n{gene} normalized expression", fontsize=viz.FONTSIZE["title"], pad=2)
        axB.set_xlabel("UMAP1")
        axB.set_ylabel("UMAP2")
        axB.set_xticks([])
        axB.set_yticks([])
        sns.despine(ax=axB, left=True, bottom=True)
    else:
        axB.axis("off")
        axB.text(0.5, 0.5, "Keratinocyte UMAP not available", ha="center", va="center")

    plot_expression_violin(
        ax=axC,
        obs_ct=obs_ct,
        gene=gene,
        cell_type=cell_type,
        cond_colors=cond_colors,
    )

    plot_paired_metric(
        ax=axD,
        pair_df=pair_pct,
        metric_col_le="LE",
        metric_col_normal="Normal",
        title=f"Percent {gene}+ {cell_type}",
        y_label="Percent positive",
        p_value=float(stats_df.at[0, "paired_ttest_pct_p"]),
        cond_colors=cond_colors,
        annotate_counts=True,
    )
    axD.text(
        0.02,
        0.02,
        "Positive call: raw counts > 0\nBlack diamonds mark condition means",
        transform=axD.transAxes,
        ha="left",
        va="bottom",
        fontsize=viz.FONTSIZE["annotation"],
        color="#555555",
    )

    viz.add_panel_label(axA, "A")
    viz.add_panel_label(axB, "B")
    viz.add_panel_label(axC, "C")
    viz.add_panel_label(axD, "D")
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.08, top=0.94, wspace=0.30, hspace=0.30)
    viz.save_figure(fig, DOCS_OUT / artifact_stem("overview", prefix), tight=False)
    plt.close(fig)


def build_supplementary_figure(
    prefix: str,
    gene: str,
    cell_type: str,
    summary: dict[str, object],
) -> None:
    obs_ct = summary["obs_ct"]
    pair_me = summary["pair_me"]
    pair_pb = summary["pair_pb"]
    sensitivity_df = summary["sensitivity_df"]
    stats_df = summary["stats_df"]
    cond_colors = viz.get_condition_colors()

    fig = plt.figure(figsize=(6.85, 5.55))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0], wspace=0.32, hspace=0.34)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    plot_patient_stratified_violin(
        ax=ax1,
        obs_ct=obs_ct,
        gene=gene,
        cell_type=cell_type,
        cond_colors=cond_colors,
    )
    ax1.set_title(f"{gene} in keratinocytes\nby paired donor", fontsize=8.4, pad=4)
    if ax1.legend_ is not None:
        handles, labels = ax1.get_legend_handles_labels()
        dedup = {}
        for handle, label in zip(handles, labels):
            if label in DISPLAY_CONDITION_ORDER and label not in dedup:
                dedup[label] = handle
        ax1.legend(
            dedup.values(),
            dedup.keys(),
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
            frameon=False,
            fontsize=viz.FONTSIZE["legend"],
            title=None,
        )

    plot_paired_metric(
        ax=ax2,
        pair_df=pair_me,
        metric_col_le="LE",
        metric_col_normal="Normal",
        title=f"Mean {gene} per donor",
        y_label="Mean normalized expression",
        p_value=float(stats_df.at[0, "paired_ttest_mean_expr_p"]),
        cond_colors=cond_colors,
        annotate_counts=False,
        show_p_value_in_title=False,
    )

    plot_paired_metric(
        ax=ax3,
        pair_df=pair_pb,
        metric_col_le="LE",
        metric_col_normal="Normal",
        title=f"{gene} pseudobulk",
        y_label="log1p CPM",
        p_value=float(stats_df.at[0, "paired_ttest_log1p_cpm_p"]),
        cond_colors=cond_colors,
        annotate_counts=False,
        show_p_value_in_title=False,
    )

    heat_df = (
        sensitivity_df[sensitivity_df["patient_label"].ne("All")]
        .pivot(index="threshold_label", columns="patient_label", values="delta_LE_minus_Normal_pp")
        .reindex([f"counts > {k}" for k in [0, 1, 2, 3]])
    )
    sns.heatmap(
        heat_df,
        ax=ax4,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "LE - Normal (pp)"},
    )
    ax4.set_title("Threshold sensitivity", fontsize=8.4, pad=4)
    ax4.set_xlabel("Donor")
    ax4.set_ylabel("Positive-call rule")

    viz.add_panel_label(ax1, "E")
    viz.add_panel_label(ax2, "F")
    viz.add_panel_label(ax3, "G")
    viz.add_panel_label(ax4, "H")
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.11, top=0.92, wspace=0.34, hspace=0.36)
    save_figure_exports(
        fig,
        alias_stem="supplement",
        export_stem=f"{analysis_slug(gene, cell_type)}_supplement",
        export_dir=COMPOSITE_OUT,
        prefix=prefix,
        tight=False,
    )
    plt.close(fig)


def export_main_panel_figures(
    prefix: str,
    gene: str,
    cell_type: str,
    adata: sc.AnnData,
    summary: dict[str, object],
    panel_set: str = "all",
) -> None:
    obs_ct = summary["obs_ct"]
    pair_pct = summary["pair_pct"]
    stats_df = summary["stats_df"]
    ct_mask = summary["ct_mask"]
    full_mask = summary["full_mask"]
    obs_full = summary["obs_full"].loc[full_mask].copy().reset_index(drop=True)
    adata_full = adata[full_mask, :].copy()

    counts_full, _ = extract_counts(adata_full, gene)
    expr_full = extract_normalized_expr(adata_full, gene)
    obs_full["positive"] = counts_full > 0
    obs_full["expr"] = expr_full

    celltype_summary = (
        obs_full.groupby("cell_type", observed=True)
        .agg(
            total_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            mean_expr=("expr", "mean"),
        )
        .reset_index()
    )
    celltype_summary["pct_positive"] = celltype_summary["positive_cells"] / celltype_summary["total_cells"] * 100
    celltype_summary["short_label"] = celltype_summary["cell_type"].map(short_cell_type_label)
    celltype_summary = celltype_summary.sort_values(
        ["pct_positive", "mean_expr"], ascending=[True, True]
    ).reset_index(drop=True)
    celltype_summary["ypos"] = np.arange(len(celltype_summary))

    cond_colors = viz.get_condition_colors()
    keratinocyte_color = viz.get_cell_type_colors([cell_type])[cell_type]
    slug = analysis_slug(gene, cell_type)

    if "X_umap" in adata_full.obsm and panel_selected(panel_set, "A"):
        fig, ax = plt.subplots(figsize=(3.2, 3.9))
        coords_full = adata_full.obsm["X_umap"]
        vmax = float(np.quantile(expr_full[expr_full > 0], 0.995)) if np.any(expr_full > 0) else 1.0
        ax.scatter(coords_full[:, 0], coords_full[:, 1], s=3.0, color="#DADADA", alpha=0.42, linewidths=0, rasterized=True)
        pos = expr_full > 0
        handle = ax.scatter(
            coords_full[pos, 0],
            coords_full[pos, 1],
            s=5.5,
            c=np.clip(expr_full[pos], 0, vmax),
            cmap="magma",
            vmin=0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        coords_ct = adata.obsm["X_umap"][ct_mask]
        ax.scatter(
            coords_ct[:, 0],
            coords_ct[:, 1],
            s=4.0,
            color=keratinocyte_color,
            alpha=0.14,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title("All-cell UMAP", fontsize=8.8, pad=5)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
        sns.despine(ax=ax, left=True, bottom=True)
        cax = inset_axes(ax, width="6%", height="42%", loc="center right", borderpad=1.0)
        cbar = fig.colorbar(handle, cax=cax)
        cbar.ax.tick_params(labelsize=viz.FONTSIZE["tick"] - 0.3)
        cbar.ax.set_title(gene, fontsize=viz.FONTSIZE["axis_title"] - 0.8, pad=3)

        inset = inset_axes(ax, width="39%", height="39%", loc="lower left", borderpad=0.7)
        expr_ct = obs_ct["expr"].to_numpy(dtype=float)
        inset.scatter(coords_ct[:, 0], coords_ct[:, 1], s=6, color="#E8E8E8", alpha=0.65, linewidths=0, rasterized=True)
        inset_pos = expr_ct > 0
        inset.scatter(
            coords_ct[inset_pos, 0],
            coords_ct[inset_pos, 1],
            s=8,
            c=np.clip(expr_ct[inset_pos], 0, vmax),
            cmap="magma",
            vmin=0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        inset.set_xticks([])
        inset.set_yticks([])
        for spine in inset.spines.values():
            spine.set_color("#444444")
            spine.set_linewidth(0.55)
        fig.subplots_adjust(left=0.11, right=0.92, bottom=0.10, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_A_all_cell_umap",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "B"):
        fig, ax = plt.subplots(figsize=(3.15, 3.95))
        ax.hlines(
            y=celltype_summary["ypos"],
            xmin=0.08,
            xmax=celltype_summary["pct_positive"].clip(lower=0.08),
            color="#D1D1D1",
            linewidth=1.1,
            zorder=1,
        )
        point_colors = [keratinocyte_color if ct == cell_type else "#B7B7B7" for ct in celltype_summary["cell_type"]]
        point_edges = ["#1A1A1A" if ct == cell_type else "white" for ct in celltype_summary["cell_type"]]
        point_sizes = [62 if ct == cell_type else 36 for ct in celltype_summary["cell_type"]]
        ax.scatter(
            celltype_summary["pct_positive"].clip(lower=0.08),
            celltype_summary["ypos"],
            s=point_sizes,
            c=point_colors,
            edgecolors=point_edges,
            linewidths=0.6,
            zorder=2,
        )
        ax.set_xscale("log")
        ax.set_xlim(0.08, 60)
        ax.set_xticks([0.1, 0.3, 1, 3, 10, 30])
        ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"])
        ax.set_yticks(celltype_summary["ypos"])
        ax.set_yticklabels(celltype_summary["short_label"])
        ax.tick_params(axis="y", labelsize=6.2)
        ax.tick_params(axis="x", labelsize=6.6)
        for tick, ct in zip(ax.get_yticklabels(), celltype_summary["cell_type"]):
            if ct == cell_type:
                tick.set_fontweight("bold")
                tick.set_color("#1A1A1A")
        ax.set_title("Cell-type specificity", fontsize=8.8, pad=5)
        ax.set_xlabel(f"Percent {gene}+ cells (log scale)")
        ax.grid(axis="x", color="#EFEFEF", linewidth=0.6)
        sns.despine(ax=ax, left=False, bottom=False)
        fig.subplots_adjust(left=0.42, right=0.98, bottom=0.14, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_B_cell_type_specificity",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "C"):
        fig, ax = plt.subplots(figsize=(3.2, 3.4))
        plot_expression_violin(
            ax=ax,
            obs_ct=obs_ct,
            gene=gene,
            cell_type=cell_type,
            cond_colors=cond_colors,
            show_inline_counts=False,
            note_text=None,
            xticklabels=DISPLAY_CONDITION_ORDER,
        )
        ax.set_title(f"{gene} in keratinocytes", fontsize=8.8, pad=5)
        ax.grid(axis="y", color="#F1F1F1", linewidth=0.55)
        fig.subplots_adjust(left=0.16, right=0.98, bottom=0.15, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_C_keratinocyte_expression_violin",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "D"):
        fig, ax = plt.subplots(figsize=(3.2, 3.4))
        plot_paired_metric(
            ax=ax,
            pair_df=pair_pct,
            metric_col_le="LE",
            metric_col_normal="Normal",
            title=f"% {gene}+ keratinocytes",
            y_label="Percent positive",
            p_value=float(stats_df.at[0, "paired_ttest_pct_p"]),
            cond_colors=cond_colors,
            annotate_counts=False,
            annotate_donors=True,
            xticklabels=DISPLAY_CONDITION_ORDER,
        )
        fig.subplots_adjust(left=0.17, right=0.92, bottom=0.15, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_D_donor_paired_quantification",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "table"):
        fig, ax = plt.subplots(figsize=(6.0, 2.25))
        plot_count_table(ax, pair_pct, cond_colors, show_title=True)
        fig.subplots_adjust(left=0.03, right=0.99, bottom=0.12, top=0.93)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_donor_summary_table",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)


def export_supplement_panel_figures(
    prefix: str,
    gene: str,
    cell_type: str,
    summary: dict[str, object],
    panel_set: str = "all",
) -> None:
    obs_ct = summary["obs_ct"]
    pair_me = summary["pair_me"]
    pair_pb = summary["pair_pb"]
    sensitivity_df = summary["sensitivity_df"]
    stats_df = summary["stats_df"]
    cond_colors = viz.get_condition_colors()
    slug = analysis_slug(gene, cell_type)

    if panel_selected(panel_set, "E"):
        fig, ax = plt.subplots(figsize=(4.75, 3.35))
        plot_patient_stratified_violin(
            ax=ax,
            obs_ct=obs_ct,
            gene=gene,
            cell_type=cell_type,
            cond_colors=cond_colors,
        )
        ax.set_title(f"{gene} in keratinocytes\nby paired donor", fontsize=8.8, pad=4)
        if ax.legend_ is not None:
            handles, labels = ax.get_legend_handles_labels()
            dedup = {}
            for handle, label in zip(handles, labels):
                if label in DISPLAY_CONDITION_ORDER and label not in dedup:
                    dedup[label] = handle
            ax.legend(
                dedup.values(),
                dedup.keys(),
                loc="upper right",
                bbox_to_anchor=(0.98, 0.98),
                frameon=False,
                fontsize=viz.FONTSIZE["legend"],
                title=None,
            )
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.15, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_E_patient_stratified_violin",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "F"):
        fig, ax = plt.subplots(figsize=(3.2, 3.35))
        plot_paired_metric(
            ax=ax,
            pair_df=pair_me,
            metric_col_le="LE",
            metric_col_normal="Normal",
            title=f"Mean {gene} per donor",
            y_label="Mean normalized expression",
            p_value=float(stats_df.at[0, "paired_ttest_mean_expr_p"]),
            cond_colors=cond_colors,
            annotate_counts=False,
            annotate_donors=True,
        )
        fig.subplots_adjust(left=0.18, right=0.90, bottom=0.15, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_F_mean_expression_paired",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "G"):
        fig, ax = plt.subplots(figsize=(3.2, 3.35))
        plot_paired_metric(
            ax=ax,
            pair_df=pair_pb,
            metric_col_le="LE",
            metric_col_normal="Normal",
            title=f"{gene} pseudobulk",
            y_label="log1p CPM",
            p_value=float(stats_df.at[0, "paired_ttest_log1p_cpm_p"]),
            cond_colors=cond_colors,
            annotate_counts=False,
            annotate_donors=True,
        )
        fig.subplots_adjust(left=0.18, right=0.90, bottom=0.15, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_G_pseudobulk_paired",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)

    if panel_selected(panel_set, "H"):
        heat_df = (
            sensitivity_df[sensitivity_df["patient_label"].ne("All")]
            .pivot(index="threshold_label", columns="patient_label", values="delta_LE_minus_Normal_pp")
            .reindex([f"counts > {k}" for k in [0, 1, 2, 3]])
        )
        fig, ax = plt.subplots(figsize=(3.45, 3.35))
        sns.heatmap(
            heat_df,
            ax=ax,
            cmap="RdBu_r",
            center=0,
            annot=True,
            fmt=".1f",
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "LE - Normal (pp)"},
        )
        ax.set_title("Threshold sensitivity", fontsize=8.8, pad=4)
        ax.set_xlabel("Donor")
        ax.set_ylabel("Positive-call rule")
        fig.subplots_adjust(left=0.20, right=0.95, bottom=0.17, top=0.90)
        save_figure_exports(
            fig,
            alias_stem=None,
            export_stem=f"{slug}_H_threshold_sensitivity_heatmap",
            export_dir=PANELS_OUT,
            prefix=prefix,
            tight=False,
        )
        plt.close(fig)


def write_outputs(prefix: str, adata_path: Path, summary: dict[str, object], gene: str, cell_type: str, patient_ids: list[int]) -> None:
    sample_summary = summary["sample_summary"]
    patient_summary = summary["patient_summary"]
    pooled_summary = summary["pooled_summary"]
    pair_pct = summary["pair_pct"]
    pair_me = summary["pair_me"]
    pair_pb = summary["pair_pb"]
    sensitivity_df = summary["sensitivity_df"]
    stats_df = summary["stats_df"]
    has_duplicate_runs = (
        sample_summary.groupby(["patient_id", "condition"]).size().max() > 1
        if len(sample_summary) > 0
        else False
    )
    duplicate_note = (
        "after collapsing duplicate runs to the patient-condition level"
        if has_duplicate_runs
        else "using only the specified sample set"
    )

    donor_summary_table = patient_summary.copy()
    donor_summary_table["positive_fraction"] = (
        donor_summary_table["positive_cells"].astype(int).astype(str)
        + "/"
        + donor_summary_table["total_cells"].astype(int).astype(str)
    )
    donor_summary_table["positive_percent"] = donor_summary_table["pct_positive"].map(lambda x: f"{x:.1f}%")
    donor_summary_table = donor_summary_table[
        [
            "patient_label",
            "condition",
            "positive_cells",
            "total_cells",
            "positive_fraction",
            "positive_percent",
            "mean_expr",
            "cpm",
            "log1p_cpm",
        ]
    ]

    sample_summary.to_csv(TABLES_OUT / f"{artifact_stem('sample_summary', prefix)}.csv", index=False)
    patient_summary.to_csv(TABLES_OUT / f"{artifact_stem('patient_condition_summary', prefix)}.csv", index=False)
    donor_summary_table.to_csv(TABLES_OUT / f"{artifact_stem('donor_summary', prefix)}.csv", index=False)
    pooled_summary.to_csv(TABLES_OUT / f"{artifact_stem('pooled_summary', prefix)}.csv", index=False)
    pair_pct.to_csv(TABLES_OUT / f"{artifact_stem('patient_pairs_pct_positive', prefix)}.csv", index=False)
    pair_me.to_csv(TABLES_OUT / f"{artifact_stem('patient_pairs_mean_expression', prefix)}.csv", index=False)
    pair_pb.to_csv(TABLES_OUT / f"{artifact_stem('patient_pairs_pseudobulk', prefix)}.csv", index=False)
    sensitivity_df.to_csv(TABLES_OUT / f"{artifact_stem('threshold_sensitivity', prefix)}.csv", index=False)
    stats_df.to_csv(TABLES_OUT / f"{artifact_stem('stats', prefix)}.csv", index=False)

    main_caption_lines = [
        "Figure G. F2RL1/PAR2 expression in keratinocytes from the paired scRNA-seq cohort used for the keratinocyte/PAR2 manuscript analysis.",
        "Panels show the all-cell UMAP with keratinocyte inset, cross-cell-type specificity of percent F2RL1-positive cells, keratinocyte expression violins, and donor-paired quantification.",
        "Normal is shown on the left and LE is shown on the right throughout.",
        (
            f"Donor-aware quantification defined F2RL1-positive keratinocytes as cells with raw counts > 0 {duplicate_note}. "
            f"The mean percentage of F2RL1-positive keratinocytes was {stats_df.at[0, 'mean_pct_Normal']:.1f}% in normal skin and "
            f"{stats_df.at[0, 'mean_pct_LE']:.1f}% in lymphedema skin "
            f"(paired t test p = {stats_df.at[0, 'paired_ttest_pct_p']:.4f}; arcsin-square-root transformed paired t test p = {stats_df.at[0, 'paired_ttest_pct_arcsin_p']:.4f}). "
            "Donor-level counts are provided in the accompanying caption PDF."
        ),
    ]
    (NOTES_OUT / f"{artifact_stem('main_panel_caption', prefix)}.md").write_text(
        "\n".join(main_caption_lines) + "\n",
        encoding="utf-8",
    )

    supplement_caption_lines = [
        f"Supplementary figure. E, normalized {gene} expression in {cell_type.lower()} stratified by paired donor and condition.",
        f"F, donor-paired mean normalized {gene} expression across all {cell_type.lower()} (paired t p = {stats_df.at[0, 'paired_ttest_mean_expr_p']:.4f}).",
        f"G, donor-paired {cell_type.lower()} pseudobulk {gene} abundance shown as log1p CPM after summing raw counts within donor-condition and scaling to CPM (paired t p = {stats_df.at[0, 'paired_ttest_log1p_cpm_p']:.4f}).",
        (
            "H, sensitivity of the LE-Normal difference to stricter positive-cell thresholds. "
            f"Paired t p values were {sensitivity_df.loc[(sensitivity_df['patient_label'].eq('All')) & (sensitivity_df['threshold_gt'].eq(0)), 'paired_ttest_p'].iloc[0]:.4f} for counts > 0, "
            f"{sensitivity_df.loc[(sensitivity_df['patient_label'].eq('All')) & (sensitivity_df['threshold_gt'].eq(1)), 'paired_ttest_p'].iloc[0]:.4f} for counts > 1, "
            f"{sensitivity_df.loc[(sensitivity_df['patient_label'].eq('All')) & (sensitivity_df['threshold_gt'].eq(2)), 'paired_ttest_p'].iloc[0]:.4f} for counts > 2, and "
            f"{sensitivity_df.loc[(sensitivity_df['patient_label'].eq('All')) & (sensitivity_df['threshold_gt'].eq(3)), 'paired_ttest_p'].iloc[0]:.4f} for counts > 3."
        ),
    ]
    (NOTES_OUT / f"{artifact_stem('supplement_caption', prefix)}.md").write_text(
        "\n".join(supplement_caption_lines) + "\n",
        encoding="utf-8",
    )

    response_lines = [
        "We quantified F2RL1-expressing keratinocytes in the paired keratinocyte manuscript cohort.",
        (
            "Using the exact cohort rerun and defining F2RL1-positive keratinocytes as cells with raw F2RL1 counts greater than 0, duplicate sequencing runs were collapsed to the patient-condition level before statistical testing."
            if has_duplicate_runs
            else "Using the exact specified sample set and defining F2RL1-positive keratinocytes as cells with raw F2RL1 counts greater than 0, donor-level percentages were compared directly across the paired samples."
        ),
        (
            f"Across paired donors, the mean percentage of F2RL1-positive keratinocytes was {stats_df.at[0, 'mean_pct_Normal']:.2f}% in normal skin and "
            f"{stats_df.at[0, 'mean_pct_LE']:.2f}% in lymphedema skin "
            f"(delta {stats_df.at[0, 'mean_delta_pct_pp']:.2f} percentage points; paired t test p={stats_df.at[0, 'paired_ttest_pct_p']:.4f}; "
            f"arcsin-square-root transformed paired t test p={stats_df.at[0, 'paired_ttest_pct_arcsin_p']:.4f})."
        ),
        "The donor-level counts and percentages are shown in the revised figure and accompanying table.",
    ]
    (DOCS_OUT / f"{artifact_stem('analysis_brief', prefix)}.md").write_text(
        "\n".join(response_lines) + "\n",
        encoding="utf-8",
    )

    report_lines = [
        "# F2RL1/PAR2 Analysis in Keratinocytes",
        "",
        "## Inputs",
        f"- Annotated object: `{adata_path}`",
        f"- Cell type: `{cell_type}`",
        f"- Gene: `{gene}`",
        "- Donor pairing: retained in the tracked summary tables",
        "- Positive-call rule: raw counts > 0",
        (
            "- Duplicate runs were collapsed to the patient-condition level for donor-aware summaries"
            if has_duplicate_runs
            else "- Only the specified samples were used; no duplicate-run collapse was required"
        ),
        "",
        "## Primary Result",
        (
            f"- Mean percent positive: LE {stats_df.at[0, 'mean_pct_LE']:.2f}% vs "
            f"Normal {stats_df.at[0, 'mean_pct_Normal']:.2f}% "
            f"(delta {stats_df.at[0, 'mean_delta_pct_pp']:.2f} percentage points; "
            f"paired t p={stats_df.at[0, 'paired_ttest_pct_p']:.4f}, "
            f"arcsin-transformed paired t p={stats_df.at[0, 'paired_ttest_pct_arcsin_p']:.4f}, "
            f"paired Wilcoxon p={stats_df.at[0, 'paired_wilcoxon_pct_p']:.4f})."
        ),
        "",
        "## Secondary Checks",
        f"- Mean normalized expression per keratinocyte: paired t p={stats_df.at[0, 'paired_ttest_mean_expr_p']:.4f}",
        f"- Keratinocyte pseudobulk log1p CPM: paired t p={stats_df.at[0, 'paired_ttest_log1p_cpm_p']:.4f}",
        f"- Pooled Fisher exact test on cells: p={stats_df.at[0, 'pooled_fisher_p']:.4f}",
        "",
        "## Interpretation",
        "- The donor-level percent-positive metric shows a consistent decrease in LE across all three paired donors.",
        "- Secondary abundance-normalized expression metrics are not significant, so the strongest statement is about the fraction of keratinocytes with detectable F2RL1 transcripts rather than a robust increase in bulk transcript abundance.",
        "",
        "## Patient-Condition Summary",
        patient_summary.to_string(index=False),
        "",
        "## Threshold Sensitivity",
        sensitivity_df.to_string(index=False),
        "",
        "## Statistics",
        stats_df.iloc[0].to_frame().to_string(header=False),
        "",
    ]
    (DOCS_OUT / f"{artifact_stem('report', prefix)}.md").write_text("\n".join(report_lines), encoding="utf-8")

    summary_lines = [
        "Advanced F2RL1/PAR2 keratinocyte analysis",
        "",
        "Sample-level summary:",
        sample_summary.to_string(index=False),
        "",
        "Patient-condition summary:",
        patient_summary.to_string(index=False),
        "",
        "Pooled summary:",
        pooled_summary.to_string(index=False),
        "",
        "Paired percent-positive summary:",
        pair_pct.to_string(index=False),
        "",
        "Paired mean-expression summary:",
        pair_me.to_string(index=False),
        "",
        "Paired pseudobulk summary:",
        pair_pb.to_string(index=False),
        "",
        "Threshold sensitivity:",
        sensitivity_df.to_string(index=False),
        "",
        "Statistics:",
        stats_df.iloc[0].to_frame().to_string(header=False),
        "",
    ]
    (DOCS_OUT / f"{artifact_stem('summary', prefix)}.md").write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    adata_path = Path(args.adata).expanduser().resolve()

    print("=" * 70)
    print("ADVANCED F2RL1/PAR2 KERATINOCYTE ANALYSIS")
    print("=" * 70)
    print(f"Input adata: {adata_path}")

    adata = sc.read_h5ad(adata_path)
    summary = summarise_metrics(
        adata=adata,
        gene=args.gene,
        cell_type=args.cell_type,
        patient_ids=args.patient_ids,
    )

    if not args.figures_only:
        write_outputs(
            prefix=args.prefix,
            adata_path=adata_path,
            summary=summary,
            gene=args.gene,
            cell_type=args.cell_type,
            patient_ids=args.patient_ids,
        )

    if args.figure_set in {"all", "supplement"} and args.panel_set == "all":
        build_supplementary_figure(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            summary=summary,
        )

    if args.figure_set in {"all", "supplement"}:
        export_supplement_panel_figures(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            summary=summary,
            panel_set=args.panel_set,
        )

    if args.figure_set in {"all", "main_panel"} and args.panel_set == "all":
        build_main_panel(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            adata=adata,
            summary=summary,
        )

    if args.figure_set in {"all", "main_panel"}:
        export_main_panel_figures(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            adata=adata,
            summary=summary,
            panel_set=args.panel_set,
        )

    if not args.figures_only:
        build_main_caption_pdf(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            summary=summary,
        )
        build_supplement_caption_pdf(
            prefix=args.prefix,
            gene=args.gene,
            cell_type=args.cell_type,
            summary=summary,
        )

    print("Saved advanced F2RL1 analysis outputs to:")
    print(f"  {FIGURES_OUT}")
    print(f"  {COMPOSITE_OUT}")
    print(f"  {PANELS_OUT}")
    print(f"  {TABLES_OUT}")
    print(f"  {NOTES_OUT}")
    print(f"  {DOCS_OUT}")


if __name__ == "__main__":
    main()


