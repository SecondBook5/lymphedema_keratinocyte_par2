# Quickstart

This repository packages the keratinocyte/PAR2-specific cohort definitions,
tracked outputs, and rerun wrappers used for this analysis. The shared
`paired_sc` workflow provides the general single-cell processing
machinery underneath that study-specific layer.

## Expected local inputs

- a local installation of `paired_sc`
- a directory containing the 10x H5 files referenced by the cohort manifest

For internal use, the matrix root can point at the existing
`Lymphadema_Mast_Cell_scrna-seq_Analysis/data/h5_files` directory.

## Typical rerun

```bash
pip install -e ../paired-single-cell-pipeline
pip install -r requirements.txt

python ./scripts/build_runtime_assets.py \
  --cohort manuscript_cohort \
  --matrix-root /path/to/h5_files \
  --validate

./scripts/run_core_rerun.sh \
  --cohort manuscript_cohort \
  --matrix-root /path/to/h5_files

./scripts/regenerate_outputs.sh --cohort manuscript_cohort
```

This sequence writes ignored runtime outputs under `work/` and refreshes the
tracked figure/table assets at the repo root.

If you only want the paper figures, run the figure-specific scripts instead of
the full output refresh:

```bash
./scripts/figure_g.sh --cohort manuscript_cohort
./scripts/figure_s2.sh --cohort manuscript_cohort
```

Panel-level rerenders for Figure G:

```bash
./scripts/figure_g_panel_a.sh --cohort manuscript_cohort
./scripts/figure_g_panel_b.sh --cohort manuscript_cohort
./scripts/figure_g_panel_c.sh --cohort manuscript_cohort
./scripts/figure_g_panel_d.sh --cohort manuscript_cohort
./scripts/figure_g_donor_table.sh --cohort manuscript_cohort
```

Supplementary panel and state-figure rerenders:

```bash
./scripts/figure_s2_panel_e.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_f.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_g.sh --cohort manuscript_cohort
./scripts/figure_s2_panel_h.sh --cohort manuscript_cohort
./scripts/figure_s2_state_overview.sh --cohort manuscript_cohort
./scripts/figure_s2_state_sensitivity.sh --cohort manuscript_cohort
./scripts/figure_s2_state_decomposition.sh --cohort manuscript_cohort
```

To rerun a different tracked cohort, replace `manuscript_cohort` with the
manifest stem you want to use from `data/manifests/`.

## What gets refreshed

- `figures/main_panel.*`
- `figures/supplement.*`
- `figures/state_overview.*`
- `figures/state_sensitivity.*`
- `figures/state_decomposition.*`
- `tables/*.csv`
- `notes/*.md`

