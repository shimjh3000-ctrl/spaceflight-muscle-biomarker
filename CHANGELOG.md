# Changelog

All notable changes to this repository are documented here. The project adheres to [Semantic Versioning](https://semver.org/) and the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## [1.0.0] — 2026-04-30

Initial public release accompanying the manuscript submission to *Scientific Reports*.

### Added
- ISS THP-1 RNA-seq pipeline (STAR v2.7.10a, featureCounts v2.0.1, DESeq2 v1.34.0).
- PRISMA 2020-compliant meta-analysis pipeline (412 studies, REML random-effects via metafor v3.8.1).
- Integration Confidence Score (ICS) framework with:
  - Min–max normalisation across 103 intersection genes
  - Bootstrap BCa confidence intervals (n = 10,000)
  - Permutation null distribution (n = 1,000)
  - Sensitivity analysis across four weighting configurations (A–D)
- C2C12 cross-tissue concordance pipeline for six simulated-microgravity datasets.
- Temporal kinetic modelling (linear / exponential / sigmoid / biphasic; AIC selection) on NASA Twins Study and JAXA cfRNA.
- Drug-Gene Interaction Network with Bliss independence synergy scoring (DGIdb v4.2.0, DrugBank v5.1.9, CMap L1000).
- Figure reproduction scripts for Figures 1–8 and Supplementary Figures S1–S8.
- Processed summary tables: 18-gene panel, ICS scores, C2C12 concordance, 726-gene consensus.
- MIT licence, CITATION.cff, README, REPRODUCING guide, DATA_SOURCES guide.
- Conda Python environment (`environment.yml`) and R `renv.lock`.

### Notes
- Raw FASTQ files are not redistributed here; all datasets are publicly available at GEO, NASA OSDR, and dbGaP — see [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md).
- A persistent Zenodo DOI will be issued at journal acceptance.
