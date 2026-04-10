# Keratinocyte F2RL1/PAR2 State Analysis

This note summarizes the deeper keratinocyte state analysis supporting the
single-cell component of *Critical Role of Keratinocytes and PAR2 in Secondary
Lymphedema Development* (Park, Pal, Chen, Shin, Garcia Nores, Baik,
Stull-Lane, Book, et al.; update citation details as needed). Public data
accession: GEO accession pending.

## Overview

- Cell type: annotated keratinocytes from the paired manuscript cohort
- Primary endpoint: percent of keratinocytes with raw `F2RL1` counts `> 0`

## Primary donor-aware result

- Mean across paired donors: Normal `41.65%` vs LE `36.55%`
- Mean LE - Normal difference: `-5.10 pp`
- Paired `t` test on donor-level percentages: `p = 0.0838`
- Arcsin-square-root transformed paired `t` test: `p = 0.0924`
- Paired Wilcoxon: `p = 0.25`
- Sign test: `p = 0.25`
- Pooled Fisher exact test: odds ratio `0.880`, `p = 0.0556`

## Alternative keratinocyte-level readouts

- Mean normalized expression paired `t`: `p = 0.6782`
- Mean expression among `F2RL1+` keratinocytes: heterogeneous across donors
- Keratinocyte pseudobulk (`log1p CPM`) paired `t`: `p = 0.9234`

These alternate expression-level summaries do not support a consistent LE-wide
increase.

## Contaminant sensitivity

- All keratinocytes: Normal `41.65%` vs LE `36.55%`, paired `t p = 0.0838`
- Excluding immune-like and myeloid-like contaminants: Normal `41.79%` vs LE `36.32%`, paired `t p = 0.0443`
- Epithelial core only: Normal `43.47%` vs LE `37.32%`, paired `t p = 0.0889`

The reverse scRNA trend is therefore not explained by obvious contaminant
clusters alone.

## Keratinocyte state analysis

Subclustering resolved `11` keratinocyte states with marker-based labels
including basal, stress-response, differentiated, junctional, cycling,
mesenchymal-like, and likely contaminant states.

- Most LE-enriched epithelial state: `KC3 Differentiated KC`, `+5.04 pp`
- Next LE-enriched epithelial state: `KC5 Mesenchymal-like KC`, `+4.15 pp`
- Most LE-depleted epithelial state: `KC0 Basal KC`, `-7.85 pp`
- Next LE-depleted epithelial state: `KC1 Stress-response KC`, `-4.15 pp`
- Largest epithelial increase in within-state `F2RL1+` fraction: `KC3 Differentiated KC`, `+5.62 pp`
- Next epithelial increase: `KC5 Mesenchymal-like KC`, `+3.37 pp`
- Largest epithelial decrease: `KC6 Cycling KC`, `-10.11 pp`
- Next epithelial decrease: `KC1 Stress-response KC`, `-8.80 pp`

## Decomposition of the overall LE-Normal difference

Using the identity `overall signal = sum(state abundance x within-state
signal)`, the overall LE - Normal difference was decomposed into a composition
effect and a within-state effect.

For percent `F2RL1+`:

- Total LE - Normal difference: `-5.10 pp`
- Composition contribution: `-0.27 pp`
- Within-state contribution: `-4.83 pp`

For mean normalized expression:

- Total LE - Normal difference: `-0.0228`
- Composition contribution: `+0.0001`
- Within-state contribution: `-0.0229`

Largest epithelial contributor to the percent-positive difference:

- `KC1 Stress-response KC`: `-5.46 pp total`
- Breakdown: `-2.49 pp` composition and `-2.97 pp` within-state

Largest epithelial contributor to the mean-expression difference:

- `KC1 Stress-response KC`: `-0.0356 total`
- Breakdown: `-0.0191` composition and `-0.0165` within-state

## Bottom line

`F2RL1/PAR2` is heterogeneously expressed across keratinocyte states. In the
paired keratinocyte manuscript cohort, the overall fraction of `F2RL1+`
keratinocytes is lower in LE than in paired Normal skin, and this pattern is
not explained by obvious contaminant clusters alone. However, LE is enriched
for differentiated and mesenchymal-like keratinocyte states, some of which
show increased within-state `F2RL1` positivity, offering a plausible
explanation for discrepancies with tissue-level assays.
