# 05 — Temporal Modelling of Atrophy Trajectories

Non-linear least-squares (NLS) modelling of expression trajectories across
0, 24, 72, 168 h to define an early–intermediate–late biomarker tier protocol.

## Goal

Fit candidate kinetic models (exponential decay, sigmoidal logistic, biphasic)
to each panel gene and select the best by AICc. Use t1/2 to assign genes to:

- **Tier 1 — Early (≤24 h)**: rapid responders for acute monitoring
- **Tier 2 — Late (≥72 h)**: chronic responders for sustained-mission tracking

## Inputs

- `data/processed/timecourse_normalised.csv` (gene × timepoint matrix)
- `data/processed/panel_18.csv` (final panel)

## Outputs

- `data/processed/nls_fits.csv` (parameters + AICc per gene per model)
- `data/processed/two_tier_assignment.csv`
- `figures/figureS7_temporal_trajectories.pdf`

## Run order

```bash
Rscript nls_model_selection.R
Rscript two_tier_protocol.R
```
