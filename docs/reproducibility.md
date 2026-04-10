# Reproducibility

This repository is the curated keratinocyte/PAR2 analysis companion. Generic
paired-tissue scRNA-seq workflow logic is provided by the shared
`paired_sc` package.

## Inputs

- curated cohort definitions in [../data/manifests](../data/manifests)
- local 10x H5 files supplied at rerun time via `--matrix-root`

## Cohort selection

The tracked manifest definitions are versioned under
[../data/manifests](../data/manifests). The public default for the manuscript
rerun is `manuscript_cohort`. Additional paired rerun definitions are preserved
alongside it for comparison or validation runs.

## Runtime layout

Package-backed reruns are written under the ignored `work/` directory:

- `work/<cohort>/results/`
- `work/<cohort>/figures/`
- `work/<cohort>/reports/`
- `work/<cohort>/logs/`

Tracked figure, table, and note assets live at the repo root and are refreshed
by [../scripts/regenerate_outputs.sh](../scripts/regenerate_outputs.sh).

For paper-figure refreshes without touching the tracked tables or notes, use:

- [../scripts/figure_g.sh](../scripts/figure_g.sh)
- [../scripts/figure_s2.sh](../scripts/figure_s2.sh)

## Key commands

```bash
python ./scripts/build_runtime_assets.py --cohort manuscript_cohort --matrix-root /path/to/h5_files --validate
./scripts/run_core_rerun.sh --cohort manuscript_cohort --matrix-root /path/to/h5_files
./scripts/regenerate_outputs.sh --cohort manuscript_cohort
```

Figure-only rerenders:

```bash
./scripts/figure_g.sh --cohort manuscript_cohort
./scripts/figure_s2.sh --cohort manuscript_cohort
```

Panel-level Figure G rerenders:

```bash
./scripts/figure_g_panel_a.sh --cohort manuscript_cohort
./scripts/figure_g_panel_b.sh --cohort manuscript_cohort
./scripts/figure_g_panel_c.sh --cohort manuscript_cohort
./scripts/figure_g_panel_d.sh --cohort manuscript_cohort
./scripts/figure_g_donor_table.sh --cohort manuscript_cohort
```

Supplementary figure rerenders:

```bash
./scripts/figure_s2_panel_e.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_f.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_g.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_h.sh --cohort manuscript_cohort
./scripts/figure_s2_state_overview.sh --cohort manuscript_cohort
./scripts/figure_s2_state_sensitivity.sh --cohort manuscript_cohort
./scripts/figure_s2_state_decomposition.sh --cohort manuscript_cohort
```

