#!/usr/bin/env python3
"""State-level keratinocyte F2RL1 analysis for deeper manuscript support."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats


THIS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = THIS_DIR.parent
FIGURES_OUT = PROJECT_ROOT / "figures"
TABLES_OUT = PROJECT_ROOT / "tables"
NOTES_OUT = PROJECT_ROOT / "notes"
DEFAULT_ADATA = PROJECT_ROOT / "work" / "manuscript_cohort" / "results" / "adata_annotated.h5ad"
DEFAULT_PREFIX = ""

for directory in [FIGURES_OUT, TABLES_OUT, NOTES_OUT]:
    directory.mkdir(parents=True, exist_ok=True)

import sys

sys.path.insert(0, str(THIS_DIR.resolve()))
from support.sample_metadata import attach_patient_metadata
from support import viz_config as viz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Keratinocyte state-level F2RL1 analysis.")
    parser.add_argument("--adata", default=str(DEFAULT_ADATA), help="Annotated h5ad path.")
    parser.add_argument("--gene", default="F2RL1", help="Target gene.")
    parser.add_argument("--cell-type", default="Keratinocytes", help="Cell type to subcluster.")
    parser.add_argument("--resolution", type=float, default=0.4, help="Leiden resolution.")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX, help="Filename prefix.")
    parser.add_argument(
        "--figures-only",
        action="store_true",
        help="Render figure assets only and skip tracked table/note exports.",
    )
    parser.add_argument(
        "--figure-set",
        choices=["all", "overview", "sensitivity", "decomposition"],
        default="all",
        help="Subset of state figures to render. Default: all",
    )
    parser.add_argument(
        "--export-marker-table",
        action="store_true",
        help="Export the full marker table. Disabled by default to keep the repo surface thin.",
    )
    return parser.parse_args()


def artifact_stem(name: str, prefix: str = "") -> str:
    prefix = prefix.strip().strip("_")
    return f"{prefix}_{name}" if prefix else name


def compact_state_label(label: str, multiline: bool = False) -> str:
    if not isinstance(label, str):
        return label
    token, remainder = label.split(" ", 1)
    short = remainder
    replacements = {
        "Stress-response KC": "Stress",
        "Differentiated KC": "Diff",
        "Junctional KC": "Junc",
        "Mesenchymal-like KC": "Mesench",
        "Cycling KC": "Cycling",
        "Immune-like contaminant": "Immune",
        "Myeloid-like contaminant": "Myeloid",
        "KRT7/KRT18-like KC": "KRT7/18",
        "Unclassified KC": "Unclass",
        "Basal KC": "Basal",
    }
    for old, new in replacements.items():
        short = short.replace(old, new)
    token = token.replace("KC", "K")
    return f"{token}\n{short}" if multiline else f"{token} {short}"


def flatten_columns(df: pd.DataFrame) -> pd.DataFrame:
    flat = df.copy()
    flat.columns = [
        col[0] if isinstance(col, tuple) and (len(col) == 1 or col[1] == "") else "_".join(str(x) for x in col if x != "")
        if isinstance(col, tuple)
        else col
        for col in flat.columns
    ]
    return flat


def extract_vector(x) -> np.ndarray:
    if sparse.issparse(x):
        return x.toarray().ravel()
    return np.asarray(x).ravel()


def extract_counts(adata: sc.AnnData, gene: str) -> np.ndarray:
    if "counts" in adata.layers and gene in adata.var_names:
        return extract_vector(adata.layers["counts"][:, adata.var_names.get_loc(gene)])
    if adata.raw is not None and gene in adata.raw.var_names:
        return extract_vector(adata.raw[:, gene].X)
    return extract_vector(adata[:, gene].X)


def extract_expr(adata: sc.AnnData, gene: str) -> np.ndarray:
    return extract_vector(adata[:, gene].X)


def safe_paired_t(x: pd.Series, y: pd.Series) -> float:
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    if len(xv) < 2:
        return np.nan
    return float(stats.ttest_rel(xv, yv).pvalue)


def safe_one_sample_t(x: pd.Series | np.ndarray) -> float:
    xv = np.asarray(x, dtype=float)
    if len(xv) < 2:
        return np.nan
    return float(stats.ttest_1samp(xv, 0.0).pvalue)


def classify_state(top_genes: list[str]) -> str:
    gene_set = set(top_genes)
    if gene_set & {"CD74", "MRC1", "FCER1G", "CST3", "ZEB2", "PID1"}:
        return "Myeloid-like contaminant"
    if gene_set & {"PTPRC", "CXCR4", "IL32", "SAMSN1", "ARHGDIB", "FYN"}:
        return "Immune-like contaminant"
    if gene_set & {"MKI67", "TOP2A", "ASPM", "CENPF", "CDK1", "TPX2"}:
        return "Cycling KC"
    if gene_set & {"KRT7", "KRT18", "CHCHD10", "PPDPF"}:
        return "KRT7/KRT18-like KC"
    if gene_set & {"VIM", "DCN", "COL6A2", "CCDC80", "IGFBP7", "CFH"}:
        return "Mesenchymal-like KC"
    if gene_set & {"ATF3", "BTG2", "RHOV", "DEGS1"}:
        return "Stress-response KC"
    if gene_set & {"FTH1", "FTL", "LGALS3", "HPGD"}:
        return "Iron/stress KC"
    if gene_set & {"KRT14", "KRT5", "KRT15", "COL17A1", "DST", "ITGB1"}:
        return "Basal KC"
    if gene_set & {"NEAT1", "MACF1", "PARD3", "TIAM1", "ERC1", "ABLIM1"}:
        return "Junctional KC"
    if gene_set & {"KRT10", "SFN", "DMKN", "KRT1", "SBSN", "LYPD3", "PERP"}:
        return "Differentiated KC"
    return "Unclassified KC"


def representative_markers(top_genes: list[str], n: int = 2) -> str:
    cleaned = [
        gene
        for gene in top_genes
        if not gene.startswith(("RPL", "RPS", "MT-", "MTRNR")) and gene not in {"MALAT1", "NEAT1"}
    ]
    keep = cleaned[:n] if cleaned else top_genes[:n]
    return "/".join(keep)


def subcluster_states(adata: sc.AnnData, gene: str, cell_type: str, resolution: float) -> dict[str, object]:
    ker = adata[adata.obs["cell_type"].eq(cell_type)].copy()
    obs = attach_patient_metadata(ker.obs[["sample_id", "condition"]].copy())
    obs["gene_count"] = extract_counts(ker, gene)
    obs["expr"] = extract_expr(ker, gene)
    obs["positive"] = obs["gene_count"] > 0

    ker_sub = ker.copy()
    sc.pp.highly_variable_genes(
        ker_sub,
        n_top_genes=min(2000, ker_sub.n_vars),
        flavor="seurat_v3",
        layer="counts" if "counts" in ker_sub.layers else None,
    )
    ker_sub = ker_sub[:, ker_sub.var["highly_variable"]].copy()
    sc.pp.scale(ker_sub, max_value=10)
    sc.tl.pca(ker_sub, svd_solver="arpack")
    sc.pp.neighbors(ker_sub, n_neighbors=15, n_pcs=min(20, ker_sub.obsm["X_pca"].shape[1]))
    sc.tl.umap(ker_sub)
    sc.tl.leiden(ker_sub, resolution=resolution, key_added="state_cluster")

    obs["state_cluster"] = ker_sub.obs["state_cluster"].astype(str).values
    ker.obs["state_cluster"] = obs["state_cluster"].values
    ker.obsm["X_state_umap"] = ker_sub.obsm["X_umap"]

    sc.tl.rank_genes_groups(ker, "state_cluster", method="wilcoxon")
    marker_df = sc.get.rank_genes_groups_df(ker, None)

    state_rows = []
    for cluster in sorted(obs["state_cluster"].unique(), key=lambda x: int(x)):
        top = marker_df.loc[marker_df["group"].eq(cluster)].head(15).copy()
        top_genes = top["names"].tolist()
        category = classify_state(top_genes)
        marker_stub = representative_markers(top_genes)
        state_label = f"KC{cluster} {category}"
        state_rows.append(
            {
                "state_cluster": cluster,
                "state_label": state_label,
                "state_category": category,
                "marker_stub": marker_stub,
                "top_markers": ", ".join(top_genes[:10]),
            }
        )

    state_meta = pd.DataFrame(state_rows).sort_values("state_cluster")
    meta_map = state_meta.set_index("state_cluster")
    obs["state_label"] = obs["state_cluster"].map(meta_map["state_label"])
    obs["state_category"] = obs["state_cluster"].map(meta_map["state_category"])
    obs["marker_stub"] = obs["state_cluster"].map(meta_map["marker_stub"])

    marker_df = marker_df.merge(state_meta[["state_cluster", "state_label", "state_category"]], left_on="group", right_on="state_cluster", how="left")

    return {
        "keratinocytes": ker,
        "obs": obs,
        "marker_df": marker_df,
        "state_meta": state_meta,
    }


def summarise_states(obs: pd.DataFrame) -> dict[str, pd.DataFrame]:
    state_condition = (
        obs.groupby(["state_cluster", "state_label", "state_category", "condition"], observed=True)
        .agg(
            n_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            pct_positive=("positive", "mean"),
            mean_expr=("expr", "mean"),
            mean_expr_positive=("expr", lambda x: float(np.mean(x[x > 0])) if np.any(x > 0) else np.nan),
        )
        .reset_index()
    )
    state_condition["pct_positive"] *= 100
    state_condition["pct_of_condition_keratinocytes"] = (
        state_condition["n_cells"]
        / state_condition.groupby("condition")["n_cells"].transform("sum")
        * 100
    )

    patient_state = (
        obs.groupby(
            ["patient_id", "patient_label", "condition", "state_cluster", "state_label", "state_category"],
            observed=True,
        )
        .agg(
            n_cells=("positive", "size"),
            positive_cells=("positive", "sum"),
            pct_positive=("positive", "mean"),
            mean_expr=("expr", "mean"),
            mean_expr_positive=("expr", lambda x: float(np.mean(x[x > 0])) if np.any(x > 0) else np.nan),
        )
        .reset_index()
    )
    patient_state["pct_positive"] *= 100
    patient_state["pct_of_patient_condition_keratinocytes"] = (
        patient_state["n_cells"]
        / patient_state.groupby(["patient_id", "condition"])["n_cells"].transform("sum")
        * 100
    )

    abundance_delta = (
        patient_state.pivot_table(
            index=["state_cluster", "state_label", "state_category"],
            columns=["patient_label", "condition"],
            values="pct_of_patient_condition_keratinocytes",
        )
        .sort_index()
    )
    positive_delta = (
        patient_state.pivot_table(
            index=["state_cluster", "state_label", "state_category"],
            columns=["patient_label", "condition"],
            values="pct_positive",
        )
        .sort_index()
    )
    meanexpr_delta = (
        patient_state.pivot_table(
            index=["state_cluster", "state_label", "state_category"],
            columns=["patient_label", "condition"],
            values="mean_expr",
        )
        .sort_index()
    )

    patients = sorted(patient_state["patient_label"].unique(), key=lambda x: int(str(x).replace("P", "")))
    abundance_delta_df = flatten_columns(abundance_delta.reset_index())[["state_cluster", "state_label", "state_category"]].copy()
    positive_delta_df = abundance_delta_df.copy()
    meanexpr_delta_df = abundance_delta_df.copy()
    stats_rows = []

    for patient in patients:
        abundance_delta_df[patient] = (
            abundance_delta[(patient, "LE")] - abundance_delta[(patient, "Normal")]
        ).to_numpy()
        positive_delta_df[patient] = (
            positive_delta[(patient, "LE")] - positive_delta[(patient, "Normal")]
        ).to_numpy()
        meanexpr_delta_df[patient] = (
            meanexpr_delta[(patient, "LE")] - meanexpr_delta[(patient, "Normal")]
        ).to_numpy()

    abundance_delta_df["mean_delta_pp"] = abundance_delta_df[patients].mean(axis=1)
    positive_delta_df["mean_delta_pp"] = positive_delta_df[patients].mean(axis=1)
    meanexpr_delta_df["mean_delta"] = meanexpr_delta_df[patients].mean(axis=1)

    for idx, row in abundance_delta_df.iterrows():
        pos_row = positive_delta_df.iloc[idx]
        expr_row = meanexpr_delta_df.iloc[idx]
        stats_rows.append(
            {
                "state_label": row["state_label"],
                "state_category": row["state_category"],
                "abundance_mean_delta_pp": float(row["mean_delta_pp"]),
                "abundance_paired_t_p": safe_one_sample_t(row[patients].to_numpy()),
                "f2rl1_positive_mean_delta_pp": float(pos_row["mean_delta_pp"]),
                "f2rl1_positive_paired_t_p": safe_one_sample_t(pos_row[patients].to_numpy()),
                "mean_expr_delta": float(expr_row["mean_delta"]),
                "mean_expr_paired_t_p": safe_one_sample_t(expr_row[patients].to_numpy()),
            }
        )

    state_stats = pd.DataFrame(stats_rows)
    return {
        "state_condition": state_condition,
        "patient_state": patient_state,
        "abundance_delta": abundance_delta_df,
        "positive_delta": positive_delta_df,
        "meanexpr_delta": meanexpr_delta_df,
        "state_stats": state_stats,
    }


def contamination_sensitivity(obs: pd.DataFrame) -> pd.DataFrame:
    scenarios = {
        "all_keratinocytes": set(),
        "exclude_immune_like": {"Immune-like contaminant", "Myeloid-like contaminant"},
        "epithelial_core_only": {"Immune-like contaminant", "Myeloid-like contaminant", "Mesenchymal-like KC", "Iron/stress KC"},
    }
    rows = []
    for scenario, excluded_categories in scenarios.items():
        sub = obs.loc[~obs["state_category"].isin(excluded_categories)].copy()
        patient_summary = (
            sub.groupby(["patient_id", "patient_label", "condition"], observed=True)
            .agg(total_cells=("positive", "size"), positive_cells=("positive", "sum"))
            .reset_index()
        )
        patient_summary["pct_positive"] = patient_summary["positive_cells"] / patient_summary["total_cells"] * 100
        pairs = patient_summary.pivot(index=["patient_id", "patient_label"], columns="condition", values="pct_positive").reset_index()
        pvalue = safe_paired_t(pairs["Normal"], pairs["LE"])
        for _, row in patient_summary.iterrows():
            rows.append(
                {
                    "scenario": scenario,
                    "patient_id": row["patient_id"],
                    "patient_label": row["patient_label"],
                    "condition": row["condition"],
                    "total_cells": row["total_cells"],
                    "positive_cells": row["positive_cells"],
                    "pct_positive": row["pct_positive"],
                    "scenario_paired_t_p": pvalue,
                }
            )
        rows.append(
            {
                "scenario": scenario,
                "patient_id": -1,
                "patient_label": "All",
                "condition": "Normal",
                "total_cells": int(patient_summary.loc[patient_summary["condition"].eq("Normal"), "total_cells"].sum()),
                "positive_cells": int(patient_summary.loc[patient_summary["condition"].eq("Normal"), "positive_cells"].sum()),
                "pct_positive": float(pairs["Normal"].mean()),
                "scenario_paired_t_p": pvalue,
            }
        )
        rows.append(
            {
                "scenario": scenario,
                "patient_id": -1,
                "patient_label": "All",
                "condition": "LE",
                "total_cells": int(patient_summary.loc[patient_summary["condition"].eq("LE"), "total_cells"].sum()),
                "positive_cells": int(patient_summary.loc[patient_summary["condition"].eq("LE"), "positive_cells"].sum()),
                "pct_positive": float(pairs["LE"].mean()),
                "scenario_paired_t_p": pvalue,
            }
        )
    return pd.DataFrame(rows)


def decompose_state_signal(patient_state: pd.DataFrame) -> dict[str, pd.DataFrame]:
    rows = []
    for (patient_id, patient_label, state_cluster, state_label, state_category), group in patient_state.groupby(
        ["patient_id", "patient_label", "state_cluster", "state_label", "state_category"],
        observed=True,
    ):
        normal = group.loc[group["condition"].eq("Normal")].iloc[0]
        le = group.loc[group["condition"].eq("LE")].iloc[0]

        abundance_normal = normal["pct_of_patient_condition_keratinocytes"] / 100.0
        abundance_le = le["pct_of_patient_condition_keratinocytes"] / 100.0
        positive_normal = normal["pct_positive"] / 100.0
        positive_le = le["pct_positive"] / 100.0
        expr_normal = normal["mean_expr"]
        expr_le = le["mean_expr"]

        comp_pos = (abundance_le - abundance_normal) * ((positive_le + positive_normal) / 2.0)
        within_pos = (positive_le - positive_normal) * ((abundance_le + abundance_normal) / 2.0)
        comp_expr = (abundance_le - abundance_normal) * ((expr_le + expr_normal) / 2.0)
        within_expr = (expr_le - expr_normal) * ((abundance_le + abundance_normal) / 2.0)

        rows.append(
            {
                "patient_id": patient_id,
                "patient_label": patient_label,
                "state_cluster": state_cluster,
                "state_label": state_label,
                "state_category": state_category,
                "composition_pos_pp": comp_pos * 100.0,
                "within_state_pos_pp": within_pos * 100.0,
                "total_pos_pp": (comp_pos + within_pos) * 100.0,
                "composition_expr": comp_expr,
                "within_state_expr": within_expr,
                "total_expr": comp_expr + within_expr,
            }
        )

    donor_state = pd.DataFrame(rows)
    mean_state = (
        donor_state.groupby(["state_cluster", "state_label", "state_category"], observed=True)[
            ["composition_pos_pp", "within_state_pos_pp", "total_pos_pp", "composition_expr", "within_state_expr", "total_expr"]
        ]
        .mean()
        .reset_index()
    )
    donor_totals = (
        donor_state.groupby(["patient_id", "patient_label"], observed=True)[
            ["composition_pos_pp", "within_state_pos_pp", "total_pos_pp", "composition_expr", "within_state_expr", "total_expr"]
        ]
        .sum()
        .reset_index()
    )
    mean_totals = donor_totals[
        ["composition_pos_pp", "within_state_pos_pp", "total_pos_pp", "composition_expr", "within_state_expr", "total_expr"]
    ].mean().to_frame().T
    mean_totals.insert(0, "patient_id", -1)
    mean_totals.insert(1, "patient_label", "Mean")
    donor_totals = pd.concat([donor_totals, mean_totals], ignore_index=True)
    return {
        "donor_state": donor_state,
        "mean_state": mean_state,
        "donor_totals": donor_totals,
    }


def state_palette(state_labels: list[str]) -> dict[str, tuple[float, float, float]]:
    colors = sns.color_palette("tab10", n_colors=len(state_labels))
    return {label: colors[i] for i, label in enumerate(state_labels)}


def plot_state_umap(ax, ker: sc.AnnData, obs: pd.DataFrame, order: list[str], palette: dict[str, tuple[float, float, float]]) -> None:
    coords = ker.obsm["X_state_umap"]
    for label in order:
        mask = obs["state_label"].eq(label).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=6, color=palette[label], linewidths=0, alpha=0.80, rasterized=True)
        centroid = np.median(coords[mask], axis=0)
        ax.text(
            centroid[0],
            centroid[1],
            label.split()[0],
            fontsize=5.8,
            fontweight="bold",
            ha="center",
            va="center",
            color="#111111",
            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
        )
    ax.set_title("Keratinocyte states", fontsize=8.8, pad=4)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, left=True, bottom=True)


def plot_marker_heatmap(ax, ker: sc.AnnData, obs: pd.DataFrame, order: list[str]) -> None:
    marker_genes = [
        "KRT14", "KRT5", "COL17A1", "DST",
        "KRT10", "SFN", "DMKN", "KRT1",
        "MKI67", "TOP2A", "ATF3", "FTH1",
        "VIM", "DCN", "KRT7", "KRT18",
        "PTPRC", "CD74", "F2RL1",
    ]
    marker_genes = [gene for gene in marker_genes if gene in ker.var_names]
    expr_df = pd.DataFrame(ker[:, marker_genes].X.toarray() if sparse.issparse(ker[:, marker_genes].X) else ker[:, marker_genes].X, columns=marker_genes)
    expr_df["state_label"] = obs["state_label"].values
    mean_expr = expr_df.groupby("state_label").mean().reindex(order)
    z = mean_expr.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) + 1e-9), axis=0)
    z.index = [compact_state_label(label, multiline=True) for label in z.index]
    sns.heatmap(
        z.T,
        ax=ax,
        cmap="vlag",
        center=0,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "z-score", "shrink": 0.94},
    )
    ax.set_title("State markers", fontsize=8.8, pad=4)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=4.8, pad=0.6)
    ax.tick_params(axis="y", labelsize=6.0)
    for tick in ax.get_xticklabels():
        tick.set_ha("center")
        tick.set_linespacing(0.88)


def plot_delta_heatmap(ax, delta_df: pd.DataFrame, value_cols: list[str], title: str, cbar_label: str) -> None:
    plot_df = delta_df.set_index("state_label")[value_cols]
    plot_df.index = [compact_state_label(label) for label in plot_df.index]
    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".1f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": cbar_label},
    )
    ax.set_title(title, fontsize=8.8, pad=4)
    ax.set_xlabel("Donor")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=6.4)
    ax.tick_params(axis="y", labelsize=5.8)


def plot_sensitivity(ax, sensitivity_df: pd.DataFrame) -> None:
    plot_df = sensitivity_df.loc[sensitivity_df["patient_label"].ne("All")].copy()
    scenario_order = ["all_keratinocytes", "exclude_immune_like", "epithelial_core_only"]
    label_map = {
        "all_keratinocytes": "All",
        "exclude_immune_like": "Exclude immune-like",
        "epithelial_core_only": "Epithelial core",
    }
    plot_df["scenario_label"] = plot_df["scenario"].map(label_map)
    colors = viz.get_condition_colors()
    for patient in sorted(plot_df["patient_label"].unique(), key=lambda x: int(x.replace("P", ""))):
        patient_df = plot_df.loc[plot_df["patient_label"].eq(patient)].copy()
        normal = patient_df.loc[patient_df["condition"].eq("Normal"), "pct_positive"].to_numpy()
        le = patient_df.loc[patient_df["condition"].eq("LE"), "pct_positive"].to_numpy()
        xs = np.arange(len(scenario_order))
        ax.plot(xs, normal, color="#9A9A9A", lw=0.9, alpha=0.7)
        ax.plot(xs, le, color="#9A9A9A", lw=0.9, alpha=0.7)
        ax.scatter(xs - 0.07, normal, s=26, color=colors["Normal"], edgecolor="white", linewidth=0.5, zorder=3)
        ax.scatter(xs + 0.07, le, s=26, color=colors["LE"], edgecolor="white", linewidth=0.5, zorder=3)
    scenario_means = (
        sensitivity_df.loc[sensitivity_df["patient_label"].ne("All")]
        .groupby(["scenario", "condition"], observed=True)["pct_positive"]
        .mean()
        .reset_index()
    )
    for i, scenario in enumerate(scenario_order):
        for condition, shift in [("Normal", -0.07), ("LE", 0.07)]:
            val = scenario_means.loc[
                scenario_means["scenario"].eq(scenario) & scenario_means["condition"].eq(condition),
                "pct_positive",
            ].iloc[0]
            ax.scatter(i + shift, val, s=58, marker="D", color="#111111", edgecolor="white", linewidth=0.5, zorder=4)
    ax.set_xticks(np.arange(len(scenario_order)))
    ax.set_xticklabels([label_map[s] for s in scenario_order], rotation=12, ha="right")
    ax.set_ylabel("Percent F2RL1+")
    ax.set_title("Contaminant sensitivity", fontsize=8.8, pad=4)
    ax.grid(axis="y", color="#EFEFEF", linewidth=0.5)
    sns.despine(ax=ax)


def plot_state_delta_scatter(ax, state_stats: pd.DataFrame) -> None:
    plot_df = state_stats.copy()
    palette = {
        "Differentiated KC": "#4E9A8F",
        "Basal KC": "#3A7EC9",
        "Junctional KC": "#7B5EA7",
        "Stress-response KC": "#E8749A",
        "Iron/stress KC": "#E8A020",
        "Cycling KC": "#D94F3D",
        "Mesenchymal-like KC": "#8B6347",
        "KRT7/KRT18-like KC": "#9EB4C4",
        "Immune-like contaminant": "#999999",
        "Myeloid-like contaminant": "#666666",
        "Unclassified KC": "#BBBBBB",
    }
    for _, row in plot_df.iterrows():
        ax.scatter(
            row["abundance_mean_delta_pp"],
            row["f2rl1_positive_mean_delta_pp"],
            s=42,
            color=palette.get(row["state_category"], "#777777"),
            edgecolor="white",
            linewidth=0.5,
        )
        ax.text(
            row["abundance_mean_delta_pp"] + 0.15,
            row["f2rl1_positive_mean_delta_pp"] + 0.15,
            row["state_label"].replace("KC", ""),
            fontsize=5.6,
            color="#333333",
        )
    ax.axhline(0, color="#A9A9A9", lw=0.7, ls="--")
    ax.axvline(0, color="#A9A9A9", lw=0.7, ls="--")
    ax.set_xlabel("LE - Normal state abundance (pp)")
    ax.set_ylabel("LE - Normal F2RL1+ fraction (pp)")
    ax.set_title("State shifts vs F2RL1 shifts", fontsize=8.8, pad=4)
    ax.grid(color="#EFEFEF", linewidth=0.5)
    sns.despine(ax=ax)


def plot_decomposition(ax, mean_state_df: pd.DataFrame, value_comp: str, value_within: str, title: str, xlabel: str) -> None:
    plot_df = mean_state_df.copy().sort_values(f"{'total_pos_pp' if value_comp.endswith('pos_pp') else 'total_expr'}")
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df[value_comp], color="#6FA8DC", alpha=0.95, label="Composition")
    ax.barh(y, plot_df[value_within], left=plot_df[value_comp], color="#E69138", alpha=0.95, label="Within-state")
    ax.axvline(0, color="#A9A9A9", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["state_label"], fontsize=5.8)
    ax.set_title(title, fontsize=8.8, pad=4)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#EFEFEF", linewidth=0.5)
    sns.despine(ax=ax, left=False)


def write_outputs(
    prefix: str,
    gene: str,
    state_meta: pd.DataFrame,
    marker_df: pd.DataFrame,
    summary: dict[str, pd.DataFrame],
    sensitivity_df: pd.DataFrame,
    decomposition: dict[str, pd.DataFrame],
    export_marker_table: bool,
) -> None:
    epithelial_categories = {
        "Basal KC",
        "Stress-response KC",
        "Differentiated KC",
        "Junctional KC",
        "Cycling KC",
        "Mesenchymal-like KC",
        "KRT7/KRT18-like KC",
        "Unclassified KC",
    }
    state_meta.to_csv(TABLES_OUT / f"{artifact_stem('state_metadata', prefix)}.csv", index=False)
    if export_marker_table:
        marker_df.to_csv(TABLES_OUT / f"{artifact_stem('state_markers', prefix)}.csv", index=False)
    summary["state_condition"].to_csv(TABLES_OUT / f"{artifact_stem('state_condition_summary', prefix)}.csv", index=False)
    summary["patient_state"].to_csv(TABLES_OUT / f"{artifact_stem('state_patient_condition_summary', prefix)}.csv", index=False)
    summary["abundance_delta"].to_csv(TABLES_OUT / f"{artifact_stem('state_abundance_delta', prefix)}.csv", index=False)
    summary["positive_delta"].to_csv(TABLES_OUT / f"{artifact_stem('state_f2rl1_positive_delta', prefix)}.csv", index=False)
    summary["meanexpr_delta"].to_csv(TABLES_OUT / f"{artifact_stem('state_meanexpr_delta', prefix)}.csv", index=False)
    summary["state_stats"].to_csv(TABLES_OUT / f"{artifact_stem('state_stats', prefix)}.csv", index=False)
    summary["state_stats"].loc[summary["state_stats"]["state_category"].isin(epithelial_categories)].to_csv(
        TABLES_OUT / f"{artifact_stem('state_stats_epithelial_core', prefix)}.csv",
        index=False,
    )
    sensitivity_df.to_csv(TABLES_OUT / f"{artifact_stem('contaminant_sensitivity', prefix)}.csv", index=False)
    decomposition["donor_state"].to_csv(TABLES_OUT / f"{artifact_stem('state_signal_decomposition_by_donor', prefix)}.csv", index=False)
    decomposition["mean_state"].to_csv(TABLES_OUT / f"{artifact_stem('state_signal_decomposition_mean', prefix)}.csv", index=False)
    decomposition["donor_totals"].to_csv(TABLES_OUT / f"{artifact_stem('state_signal_decomposition_totals', prefix)}.csv", index=False)

    state_condition = summary["state_condition"]
    state_stats = summary["state_stats"].sort_values("f2rl1_positive_mean_delta_pp")
    epithelial_stats = state_stats.loc[state_stats["state_category"].isin(epithelial_categories)].copy()
    decomposition_mean = decomposition["mean_state"].loc[
        decomposition["mean_state"]["state_category"].isin(epithelial_categories)
    ].copy()
    pos_driver = decomposition_mean.reindex(
        decomposition_mean["total_pos_pp"].abs().sort_values(ascending=False).index
    ).iloc[0]
    expr_driver = decomposition_mean.reindex(
        decomposition_mean["total_expr"].abs().sort_values(ascending=False).index
    ).iloc[0]
    top_depleted = epithelial_stats.sort_values("f2rl1_positive_mean_delta_pp").iloc[0]
    top_enriched = epithelial_stats.sort_values("f2rl1_positive_mean_delta_pp").iloc[-1]
    abundance_depleted = epithelial_stats.sort_values("abundance_mean_delta_pp").iloc[0]
    abundance_enriched = epithelial_stats.sort_values("abundance_mean_delta_pp").iloc[-1]
    sensitivity_pivot = sensitivity_df.loc[sensitivity_df["patient_label"].eq("All")].pivot(index="scenario", columns="condition", values="pct_positive")
    scenario_p = sensitivity_df.loc[sensitivity_df["patient_label"].eq("All"), ["scenario", "scenario_paired_t_p"]].drop_duplicates()
    p_map = scenario_p.set_index("scenario")["scenario_paired_t_p"].to_dict()
    total_decomp = decomposition["donor_totals"].loc[decomposition["donor_totals"]["patient_label"].eq("Mean")].iloc[0]

    lines = [
        "Keratinocyte state-level F2RL1 analysis",
        "",
        f"Subclustering resolved {len(state_meta)} keratinocyte states from the paired keratinocyte manuscript cohort.",
        "Marker-based labels indicate a mixture of differentiated, basal, junctional, stress-response, cycling, mesenchymal-like, and likely contaminant states.",
        (
            f"Among epithelial keratinocyte states, the most LE-enriched state was {abundance_enriched['state_label']} "
            f"(mean abundance delta {abundance_enriched['abundance_mean_delta_pp']:.2f} pp)."
        ),
        (
            f"The most LE-depleted epithelial state was {abundance_depleted['state_label']} "
            f"(mean abundance delta {abundance_depleted['abundance_mean_delta_pp']:.2f} pp)."
        ),
        (
            f"The epithelial state with the largest increase in within-state {gene}-positive fraction in LE was {top_enriched['state_label']} "
            f"({top_enriched['f2rl1_positive_mean_delta_pp']:.2f} pp)."
        ),
        (
            f"The epithelial state with the largest decrease in within-state {gene}-positive fraction in LE was {top_depleted['state_label']} "
            f"({top_depleted['f2rl1_positive_mean_delta_pp']:.2f} pp)."
        ),
        "",
        "State-signal decomposition of the overall LE - Normal shift:",
        (
            f"- Overall F2RL1-positive fraction delta: {total_decomp['total_pos_pp']:.2f} pp, composed of "
            f"{total_decomp['composition_pos_pp']:.2f} pp from state abundance changes and "
            f"{total_decomp['within_state_pos_pp']:.2f} pp from within-state F2RL1 shifts."
        ),
        (
            f"- Overall mean-expression delta: {total_decomp['total_expr']:.4f}, composed of "
            f"{total_decomp['composition_expr']:.4f} from state abundance changes and "
            f"{total_decomp['within_state_expr']:.4f} from within-state expression shifts."
        ),
        (
            f"- Largest epithelial contributor to the percent-positive difference: {pos_driver['state_label']} "
            f"({pos_driver['total_pos_pp']:.2f} pp total; composition {pos_driver['composition_pos_pp']:.2f} pp, "
            f"within-state {pos_driver['within_state_pos_pp']:.2f} pp)."
        ),
        (
            f"- Largest epithelial contributor to the mean-expression difference: {expr_driver['state_label']} "
            f"({expr_driver['total_expr']:.4f} total; composition {expr_driver['composition_expr']:.4f}, "
            f"within-state {expr_driver['within_state_expr']:.4f})."
        ),
        "",
        "Contaminant sensitivity of the overall F2RL1-positive fraction:",
        (
            f"- All keratinocytes: Normal {sensitivity_pivot.loc['all_keratinocytes', 'Normal']:.2f}% vs "
            f"LE {sensitivity_pivot.loc['all_keratinocytes', 'LE']:.2f}% "
            f"(paired t p={p_map['all_keratinocytes']:.4f})"
        ),
        (
            f"- Excluding immune-like contaminants: Normal {sensitivity_pivot.loc['exclude_immune_like', 'Normal']:.2f}% vs "
            f"LE {sensitivity_pivot.loc['exclude_immune_like', 'LE']:.2f}% "
            f"(paired t p={p_map['exclude_immune_like']:.4f})"
        ),
        (
            f"- Epithelial core only: Normal {sensitivity_pivot.loc['epithelial_core_only', 'Normal']:.2f}% vs "
            f"LE {sensitivity_pivot.loc['epithelial_core_only', 'LE']:.2f}% "
            f"(paired t p={p_map['epithelial_core_only']:.4f})"
        ),
        "",
        "Interpretation:",
        "The reverse scRNA trend is not explained by obvious contaminant clusters alone.",
        "Instead, F2RL1 appears heterogeneous across keratinocyte states, so tissue-level bulk RNA or protein assays may diverge from a simple single-cell percent-positive metric if state composition changes between Normal and LE.",
    ]
    (NOTES_OUT / f"{artifact_stem('state_analysis_summary', prefix)}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(str(Path(args.adata).resolve()))
    result = subcluster_states(adata, args.gene, args.cell_type, args.resolution)
    summary = summarise_states(result["obs"])
    sensitivity_df = contamination_sensitivity(result["obs"])
    decomposition = decompose_state_signal(summary["patient_state"])

    state_order = result["state_meta"]["state_label"].tolist()
    palette = state_palette(state_order)
    donor_cols = sorted(summary["abundance_delta"].filter(regex=r"^P").columns, key=lambda x: int(x.replace("P", "")))

    if args.figure_set in {"all", "overview"}:
        fig = plt.figure(figsize=(8.8, 6.35))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.04, 1.20], height_ratios=[1.0, 1.0], wspace=0.46, hspace=0.40)
        axA = fig.add_subplot(gs[0, 0])
        axB = fig.add_subplot(gs[0, 1])
        axC = fig.add_subplot(gs[1, 0])
        axD = fig.add_subplot(gs[1, 1])

        plot_state_umap(axA, result["keratinocytes"], result["obs"], state_order, palette)
        plot_marker_heatmap(axB, result["keratinocytes"], result["obs"], state_order)
        plot_delta_heatmap(axC, summary["abundance_delta"], donor_cols, "State abundance shift\nLE - Normal (pp)", "LE - Normal (pp)")
        plot_delta_heatmap(axD, summary["positive_delta"], donor_cols, f"{args.gene}+ fraction shift\nwithin state", "LE - Normal (pp)")
        viz.add_panel_label(axA, "A")
        viz.add_panel_label(axB, "B")
        viz.add_panel_label(axC, "C")
        viz.add_panel_label(axD, "D")
        fig.subplots_adjust(left=0.06, right=0.97, bottom=0.08, top=0.95)
        viz.save_figure(fig, FIGURES_OUT / artifact_stem("state_overview", args.prefix), tight=False)
        plt.close(fig)

    if args.figure_set in {"all", "sensitivity"}:
        fig = plt.figure(figsize=(7.0, 3.2))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.36)
        axE = fig.add_subplot(gs[0, 0])
        axF = fig.add_subplot(gs[0, 1])
        plot_sensitivity(axE, sensitivity_df)
        plot_state_delta_scatter(axF, summary["state_stats"])
        viz.add_panel_label(axE, "E")
        viz.add_panel_label(axF, "F")
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.18, top=0.92)
        viz.save_figure(fig, FIGURES_OUT / artifact_stem("state_sensitivity", args.prefix), tight=False)
        plt.close(fig)

    if args.figure_set in {"all", "decomposition"}:
        fig = plt.figure(figsize=(8.2, 3.9))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.42)
        axG = fig.add_subplot(gs[0, 0])
        axH = fig.add_subplot(gs[0, 1])
        epi_decomp = decomposition["mean_state"].loc[
            decomposition["mean_state"]["state_category"].isin(
                {
                    "Basal KC",
                    "Stress-response KC",
                    "Differentiated KC",
                    "Junctional KC",
                    "Cycling KC",
                    "Mesenchymal-like KC",
                    "KRT7/KRT18-like KC",
                    "Unclassified KC",
                }
            )
        ].copy()
        plot_decomposition(
            axG,
            epi_decomp,
            "composition_pos_pp",
            "within_state_pos_pp",
            "Decomposition of F2RL1+ fraction shift",
            "LE - Normal (pp)",
        )
        plot_decomposition(
            axH,
            epi_decomp,
            "composition_expr",
            "within_state_expr",
            "Decomposition of mean-expression shift",
            "LE - Normal",
        )
        viz.add_panel_label(axG, "G")
        viz.add_panel_label(axH, "H")
        handles, labels = axG.get_legend_handles_labels()
        axH.legend(handles, labels, loc="lower right", frameon=False, fontsize=6.4)
        fig.subplots_adjust(left=0.17, right=0.98, bottom=0.14, top=0.92)
        viz.save_figure(fig, FIGURES_OUT / artifact_stem("state_decomposition", args.prefix), tight=False)
        plt.close(fig)

    if not args.figures_only:
        write_outputs(
            args.prefix,
            args.gene,
            result["state_meta"],
            result["marker_df"],
            summary,
            sensitivity_df,
            decomposition,
            args.export_marker_table,
        )

    print("Saved state-level keratinocyte analysis to:")
    print(f"  {FIGURES_OUT}")
    print(f"  {TABLES_OUT}")
    print(f"  {NOTES_OUT}")


if __name__ == "__main__":
    main()
