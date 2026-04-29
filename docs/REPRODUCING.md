# Reproducing the Analysis — Step-by-Step Guide

This document walks through reproducing every result in the manuscript from raw publicly available data. Skip to **Section 5** if you only want to reproduce ICS scores from the included summary tables.

---

## 1. Prerequisites

- Linux (Ubuntu 20.04+ recommended) or macOS
- ≥ 32 GB RAM, ≥ 200 GB free disk space (raw data + intermediates)
- Conda / Mambaforge installed
- R ≥ 4.2.1 with `renv` installed
- Stable internet connection (~100 GB downloads)

```bash
# Set up environment
conda env create -f environment.yml
conda activate spaceflight-muscle

R -e "renv::restore()"
```

---

## 2. Download raw data

See [`docs/DATA_SOURCES.md`](DATA_SOURCES.md) for individual download commands.

```bash
# Example: ISS THP-1 transcriptome (primary discovery)
mkdir -p data/raw/GSE245789 && cd data/raw/GSE245789
# Download FASTQ via SRA Toolkit (pre-installed in environment.yml)
prefetch GSE245789
fasterq-dump --split-files SRR_accession_list.txt
cd ../../..
```

> **Note**: NASA Twins Study (dbGaP phs001453) requires an approved data access request.

---

## 3. Run the RNA-seq pipeline (ISS THP-1)

```bash
# Step 1: STAR alignment to GRCh38.p13 / GENCODE v38
bash 01_rnaseq_pipeline/01_star_align.sh \
     data/raw/GSE245789/ \
     references/star_index_GRCh38_v38/ \
     data/aligned/

# Step 2: featureCounts gene-level quantification
bash 01_rnaseq_pipeline/02_featurecounts.sh \
     data/aligned/ \
     references/gencode.v38.annotation.gtf \
     data/counts/

# Step 3: DESeq2 differential expression
Rscript 01_rnaseq_pipeline/03_deseq2.R \
     data/counts/ \
     data/deseq2_results/
```

**Output**: `data/deseq2_results/THP1_SF_vs_GC.csv` — 4,522 DEGs at FDR < 0.05 and |Log2FC| > 0.5.

---

## 4. PRISMA meta-analysis (412 studies → 726-gene consensus)

```bash
Rscript 02_meta_analysis/inclusion_exclusion.R       # Apply PRISMA filters
Rscript 02_meta_analysis/effect_size_pooling.R       # REML random-effects via metafor
```

**Output**: `data/726_consensus_genes.csv` — 726 genes meeting all four inclusion criteria.

---

## 5. Compute ICS scores and 18-gene panel (no raw data needed)

The repository ships with the 103-gene intersection table, so this step is independent of upstream pipelines:

```bash
Rscript 03_ics_framework/ics_compute.R               # Computes ICS for 103 genes
Rscript 03_ics_framework/bootstrap_bca.R             # 10,000 BCa bootstrap iterations
Rscript 03_ics_framework/permutation_null.R          # 1,000 shuffled-null permutations
Rscript 03_ics_framework/sensitivity_configs.R       # Configs A–D, Robustness Index
```

**Outputs**:
- `data/ics_scores.csv` — Matches Supplementary Table S2
- `data/18_gene_panel.csv` — Final ranked panel (Gold / High / Moderate tiers)
- `data/sensitivity_robustness.csv` — Tier stability across 5 configurations

---

## 6. C2C12 cross-tissue concordance

```bash
bash 04_c2c12_concordance/reanalyze_pipeline.sh      # STAR + DESeq2 on 6 datasets
Rscript 04_c2c12_concordance/kmeans_modules.R         # k=2 module assignment
Rscript 04_c2c12_concordance/pearson_correlation.R    # r values + 95% CIs
```

**Output**: `data/c2c12_concordance.csv` — 18 genes × 6 datasets concordance matrix.

---

## 7. Temporal kinetic modelling

```bash
Rscript 05_temporal_modeling/nls_model_selection.R    # Linear / exp / sigmoid / biphasic; AIC
Rscript 05_temporal_modeling/two_tier_protocol.R      # Defines fortnightly/monthly schedule
```

---

## 8. Drug-Gene Interaction Network and Bliss synergy

```bash
Rscript 06_dgi_network/dgidb_query.R                  # 15 FDA-approved compounds
python3 06_dgi_network/string_louvain.py              # STRING v11.5 + Louvain
Rscript 06_dgi_network/bliss_synergy.R                # Bliss independence scores
```

---

## 9. Reproduce all figures

```bash
# Main text figures
python3 07_figures/figure1_overview.py
Rscript  07_figures/figure2_ics.R
# ... (Figs 3–8)

# Supplementary figures
python3 07_figures/figureS1-S8/figS1_prisma.py
# ... (Figs S2–S8)
```

All figures are written to `07_figures/output/` as PNG (300 dpi) and PDF.

---

## Expected runtimes (16-core / 64 GB workstation)

| Stage | Approximate runtime |
|---|---|
| FASTQ download (full) | 3–6 hours |
| STAR alignment (16 ISS samples) | 2 hours |
| featureCounts | 15 minutes |
| DESeq2 | 5 minutes |
| PRISMA meta-analysis | 30 minutes |
| ICS framework + bootstrap | 20 minutes |
| C2C12 reanalysis (6 datasets) | 4 hours |
| Temporal modelling | 5 minutes |
| DGI / Bliss | 10 minutes |
| Figure regeneration | 15 minutes |

---

## Troubleshooting

- **`STAR: not found`** → re-activate the conda environment.
- **`metafor` not found** → run `R -e "renv::restore()"` from the repo root.
- **dbGaP access denied** → submit a data access request at https://dbgap.ncbi.nlm.nih.gov/ for accession phs001453.
- **STRING API rate limit** → the script auto-retries with exponential backoff; if it still fails, set `STRINGDB_DELAY=2` before running.

If you encounter other issues, please open a GitHub issue or contact the corresponding authors.
