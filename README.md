# Spaceflight Muscle Biomarker — Analysis Code Repository

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/Zenodo-DOI%20pending-lightgrey)](https://zenodo.org/)
[![Made with R](https://img.shields.io/badge/R-%E2%89%A54.2.1-276DC3?logo=r)](https://www.r-project.org/)
[![Made with Python](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB?logo=python)](https://www.python.org/)

Analysis code, processed data summaries, and figure-reproduction scripts for the manuscript:

> **Integrative transcriptomics identifies an 18-gene blood-accessible biomarker panel for real-time monitoring of spaceflight-induced muscle atrophy**
> Shim JH, Lee JY, Han J, Woo SW, Lim JJ, Kim HS.


Department of Anatomy, Korea University College of Medicine, Seoul, Republic of Korea.

---

## Overview

This repository contains the full computational pipeline that produced the 18-gene blood-accessible biomarker panel reported in the manuscript:

1. **ISS THP-1 RNA-seq processing** — STAR alignment, featureCounts, DESeq2 differential expression (`01_rnaseq_pipeline/`)
2. **PRISMA-compliant meta-analysis** of 412 spaceflight studies → 726-gene consensus set (`02_meta_analysis/`)
3. **Integration Confidence Score (ICS) framework** with bootstrap, permutation, and sensitivity analysis (`03_ics_framework/`)
4. **C2C12 cross-tissue concordance** across six simulated-microgravity datasets (`04_c2c12_concordance/`)
5. **Temporal kinetic modelling** on NASA Twins Study + JAXA cfRNA (`05_temporal_modeling/`)
6. **Drug-Gene Interaction Network and Bliss synergy analysis** (`06_dgi_network/`)
7. **Figure reproduction scripts** for Figs 1–8 and Supplementary Figs S1–S8 (`07_figures/`)

---

## Repository structure

```
spaceflight-muscle-biomarker/
├── 01_rnaseq_pipeline/      # ISS THP-1 alignment + DESeq2
├── 02_meta_analysis/        # PRISMA → 726-gene consensus
├── 03_ics_framework/        # ICS, bootstrap (n=10,000), permutation (n=1,000)
├── 04_c2c12_concordance/    # 6-dataset cross-tissue analysis
├── 05_temporal_modeling/    # nls model selection, two-tier protocol
├── 06_dgi_network/          # DGIdb / DrugBank / CMap, Bliss synergy
├── 07_figures/              # Fig 1–8 + Fig S1–S8 reproduction
├── data/                    # Processed summary tables (CSV)
├── docs/                    # Reproducing guide, data sources
├── environment.yml          # Conda Python environment
├── renv.lock                # R package lock file
├── CITATION.cff             # Citation metadata
├── LICENSE                  # MIT
└── CHANGELOG.md             # Version history
```

---

## Data sources (publicly available)

Raw sequencing data are **not** included in this repository. All datasets are deposited in public archives and can be downloaded directly:

| Dataset | Accession | URL |
|---|---|---|
| Primary ISS THP-1 transcriptome | GEO **GSE245789** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245789 |
| Primary ISS THP-1 (mirror) | NASA OSDR **OSD-590** | https://osdr.nasa.gov/bio/repo/data/studies/OSD-590 |
| Rodent Research-1 | NASA OSDR **OSD-104** | https://osdr.nasa.gov/bio/repo/data/studies/OSD-104 |
| Rodent Research-23 (quadriceps) | NASA OSDR **OSD-379** | https://osdr.nasa.gov/bio/repo/data/studies/OSD-379 |
| NASA Twins Study | dbGaP **phs001453** | https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001453 (controlled access) |
| Inspiration4 | GEO **GSE193270** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193270 |
| JAXA cfRNA | GEO **GSE185989** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185989 |
| Polaris Dawn | GEO **GSE194192** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194192 |
| C2C12 #1 | GEO **GSE184765** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184765 |
| C2C12 #2 | GEO **GSE165938** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165938 |
| C2C12 #3 | GEO **GSE187458** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE187458 |
| C2C12 #4 | GEO **GSE196720** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196720 |
| C2C12 #5 | GEO **GSE201839** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201839 |
| C2C12 #6 | GEO **GSE178903** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178903 |

See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) for download commands.

---

## Quick start — reproduce key results

### 1. Clone the repository

```bash
git clone https://github.com/shimlab-kaeri/spaceflight-muscle-biomarker.git
cd spaceflight-muscle-biomarker
```

### 2. Set up the analysis environment

**Python (Conda):**
```bash
conda env create -f environment.yml
conda activate spaceflight-muscle
```

**R (renv):**
```r
install.packages("renv")
renv::restore()
```

### 3. Reproduce the 18-gene panel ICS scores

```bash
# Compute ICS from intersection genes
Rscript 03_ics_framework/ics_compute.R

# Bootstrap confidence intervals (n = 10,000 BCa)
Rscript 03_ics_framework/bootstrap_bca.R

# Sensitivity analysis (configurations A–D)
Rscript 03_ics_framework/sensitivity_configs.R
```

Outputs are written to `data/ics_scores.csv` and matches **Supplementary Table S2**.

For full step-by-step instructions including raw-data preprocessing, see [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

---

## Key outputs included in this repository

- `data/18_gene_panel.csv` — Final 18-gene biomarker panel with tier, ICS, Log2FC, multi-mission validation
- `data/ics_scores.csv` — All 103 intersection-gene ICS scores (matches Supplementary Table S2)
- `data/726_consensus_genes.csv` — Full meta-analysis consensus set
- `data/c2c12_concordance.csv` — Cross-tissue Pearson r values for 18 genes × 6 datasets

---

## Citation

If you use this code or data in your research, please cite the manuscript and this repository:

```bibtex
@article{Shim2026spaceflight,
  title   = {Integrative transcriptomics identifies an 18-gene blood-accessible
             biomarker panel for real-time monitoring of spaceflight-induced
             muscle atrophy},
  author  = {Shim, Jae Ho and Lee, Ji Yeon and Han, Jihye and
             Woo, Sang Woo and Lim, Jin Ju and Kim, Hyeon Soo},
  journal = {},
  year    = {2026},
 }
```

A persistent **Zenodo DOI** for this repository will be issued at the time of acceptance and added here.

GitHub also recognises the [`CITATION.cff`](CITATION.cff) file: click "Cite this repository" in the right-hand sidebar.

---

## License

This repository is released under the [MIT License](LICENSE). You are free to use, modify, and redistribute the code with attribution.

Datasets retrieved from GEO, NASA OSDR, and dbGaP retain their original data-use agreements; please consult the respective archives for access conditions.

---

## Contact

**Prof. Jae Ho Shim, MD, PhD**
Department of Anatomy, Korea University College of Medicine
Seoul 02841, Republic of Korea
[shimjh3000@korea.ac.kr](mailto:shimjh3000@korea.ac.kr)

**Prof. Hyeon Soo Kim, MD, PhD** (corresponding)
[anatomykim@korea.ac.kr](mailto:anatomykim@korea.ac.kr)

---

## Acknowledgements

We thank NASA GeneLab / Open Science Data Repository (OSDR), NCBI GEO, dbGaP, and JAXA for hosting the publicly available datasets used in this study, and the original research teams for their data deposition.
