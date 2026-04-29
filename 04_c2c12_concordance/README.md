# 04 — C2C12 Microgravity Concordance Analysis

Cross-validation of the in vivo–derived 18-gene panel against an independent
in vitro C2C12 myotube simulated-microgravity dataset (RWV bioreactor; 72 h).

## Goal

Test whether (i) k-means co-expression modules in the C2C12 dataset recapitulate
in vivo functional clusters (atrogenes, contractile, IGF/PI3K-AKT, ROS) and
(ii) per-gene log2FC values correlate between systems (target Pearson r ≥ 0.70
for the 18-gene panel).

## Inputs

- Raw FASTQ from independent C2C12 RWV experiment (3 vs 3, 72 h)
- 18-gene panel (`data/processed/panel_18.csv` from stage 03)

## Outputs

- `data/processed/c2c12_log2fc.csv`
- `data/processed/c2c12_kmeans_modules.csv`
- `figures/figureS6_c2c12_concordance.pdf`

## Run order

```bash
bash reanalyze_pipeline.sh             # STAR + featureCounts + DESeq2
Rscript kmeans_modules.R               # k-means clustering of co-expression
Rscript pearson_correlation.R          # in vivo vs in vitro log2FC correlation
```

## Expected results

- 4 co-expression modules matching atrogene / contractile / IGF-AKT / ROS
- Pearson r = 0.78 across the 18-gene panel (manuscript Fig. 4)
