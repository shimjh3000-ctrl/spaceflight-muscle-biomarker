# 02 — PRISMA 2020 Meta-Analysis (412 studies → 726-gene consensus)

PRISMA 2020-compliant systematic review and random-effects meta-analysis producing the **S726 consensus set**.

## Search strategy

- Databases: PubMed, EMBASE, Web of Science
- Time window: 2000–2024
- Records identified: 4,127 → screened → 412 included
- Full PRISMA flow: see Supplementary Figure S1 and Supplementary Checklist 1
- Search terms: see `prisma_search_strategy.md`

## Inclusion criteria

1. FDR < 0.05 reported
2. |Log2FC| > 0.5 reported or computable
3. Gene appears in ≥ 3 independent studies
4. Concordant fold-change direction in ≥ 66% of studies
5. DisGeNET disease-association score ≥ 0.3

## Statistical framework

- Effect-size pooling: **REML random-effects** model (`metafor::rma` v3.8.1)
- Heterogeneity reported: I², τ², Q-statistic
- All genes: I² < 75% (range 8–67%)
- Egger's regression for publication bias

## Scripts

- `prisma_search_strategy.md`     — full search strings and screening protocol
- `inclusion_exclusion.R`         — apply 4 criteria, build inclusion table
- `effect_size_pooling.R`         — REML random-effects pooling for each gene

## Output

`data/726_consensus_genes.csv` — meta-analysis summary statistics for 726 genes.
