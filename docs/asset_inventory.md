# Asset Inventory

This repository keeps only curated keratinocyte/PAR2 assets.

## Figures

- `figures/main_panel.*`
  - primary multi-panel figure for the keratinocyte `F2RL1 / PAR2` analysis
- `figures/supplement.*`
  - supplementary supporting panels for expression, pseudobulk, and threshold
    sensitivity
- `figures/state_overview.*`
  - keratinocyte state map, marker overview, abundance shift, and within-state
    `F2RL1` shift
- `figures/state_sensitivity.*`
  - sensitivity analyses excluding likely contaminant states
- `figures/state_decomposition.*`
  - decomposition of total LE-Normal shift into state-composition and
    within-state signal components
- `figures/panels/`
  - standalone panels suitable for figure rearrangement
- `figures/composite/`
  - alternate composite exports retained for convenience

## Tables

- `tables/donor_summary.csv`
  - donor-level `F2RL1+ / total` keratinocyte counts and proportions
- `tables/stats.csv`
  - paired statistical summaries for the main keratinocyte comparison
- `tables/sample_summary.csv`
  - sample-level keratinocyte summaries
- `tables/pooled_summary.csv`
  - pooled condition summary
- `tables/threshold_sensitivity.csv`
  - threshold-dependent sensitivity analysis
- `tables/state_*.csv`
  - state-level abundance, positivity, mean-expression, and decomposition
    summaries

## Notes

- `notes/main_panel_caption.md`
- `notes/supplement_caption.md`
- `notes/state_analysis.md`
- `notes/state_analysis_summary.md`

## Regeneration

- `scripts/figure_g.sh`
  - refreshes the main manuscript figure and its standalone panel exports
- `scripts/figure_g_panel_a.sh`, `scripts/figure_g_panel_b.sh`,
  `scripts/figure_g_panel_c.sh`, `scripts/figure_g_panel_d.sh`
  - refresh specific standalone Figure G panels
- `scripts/figure_g_donor_table.sh`
  - refreshes the donor summary table panel paired with Figure G
- `scripts/figure_s2.sh`
  - refreshes the supplementary figure and the state-analysis figure set
- `scripts/figure_s2_panel_e.sh`, `scripts/figure_s2_panel_f.sh`,
  `scripts/figure_s2_panel_g.sh`, `scripts/figure_s2_panel_h.sh`
  - refresh specific standalone supplementary panels
- `scripts/figure_s2_state_overview.sh`
  - refreshes the tracked state overview figure
- `scripts/figure_s2_state_sensitivity.sh`
  - refreshes the tracked state sensitivity figure
- `scripts/figure_s2_state_decomposition.sh`
  - refreshes the tracked state decomposition figure
- `scripts/regenerate_outputs.sh`
  - refreshes the tracked figures, tables, notes, and generated report summaries
