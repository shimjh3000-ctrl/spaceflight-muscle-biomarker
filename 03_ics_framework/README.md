# 03 — Integration Confidence Score (ICS) Framework

The ICS framework combines discovery-level evidence (S4522) with cross-study consistency (S726) to rank the 103 intersection genes and nominate the **18-gene validated biomarker panel**.

## Formulation

```
ICS  = 0.6 × S4522 + 0.4 × S726
S4522 = 0.6 × Log2FC_norm + 0.4 × (-log10 FDR_norm)
S726  = 0.5 × N_studies_norm + 0.25 × DisGeNET_norm + 0.25 × Category_norm
```

All sub-scores are min–max normalised to [0, 10] across the 103 intersection genes.

## Tier definitions

| Tier | ICS range |
|---|---|
| **Gold** | ≥ 8.0 |
| **High** | 7.0 – 7.9 |
| **Moderate** | 6.0 – 6.9 |

## Scripts

| Script | Purpose |
|---|---|
| `ics_compute.R`         | Compute ICS for 103 intersection genes |
| `bootstrap_bca.R`       | n = 10,000 BCa bootstrap confidence intervals |
| `permutation_null.R`    | n = 1,000 permutation null distribution; empirical p-values |
| `sensitivity_configs.R` | Four alternative weightings (A–D) + Robustness Index |

## Sensitivity configurations

| Config | S4522 weight | S726 weight | Notes |
|---|---|---|---|
| Primary | 0.6 | 0.4 | Reported in main text |
| A | 0.5 | 0.5 | Equal weighting |
| B | 0.7 | 0.3 | Discovery-heavy |
| C | 0.4 | 0.6 | Meta-analysis-heavy |
| D | 0.6 | 0.4 with reweighted S726 sub-scores | DisGeNET emphasised |

**Robustness Index** = fraction of genes whose tier classification is identical to Primary across all configurations. All Gold-tier genes: Robustness Index = 1.00 (no reclassification).

## Outputs

- `data/ics_scores.csv` (matches Supplementary Table S2)
- `data/18_gene_panel.csv`
- `data/sensitivity_robustness.csv`
- `data/bootstrap_ci.csv` — TNNT3 9.4 (9.1–9.7), TRIM63 9.2 (8.9–9.5), ACTC1 8.7 (8.4–9.0), IGF2 8.5 (8.2–8.8)
