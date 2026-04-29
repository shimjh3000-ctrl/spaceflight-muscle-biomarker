# 07 — Figure Generation

Scripts that render every main and supplementary figure published in the paper.
Each figure has its own subdirectory and entry script; all read from
`data/processed/` and write PDF/PNG to `figures/`.

## Main figures

| Figure  | Description                                              | Script                                          |
|---------|----------------------------------------------------------|-------------------------------------------------|
| Fig. 1  | Study design and analytical workflow                     | `figure1_workflow.R`                            |
| Fig. 2  | Volcano + heatmap of 4,522 DEGs and 103 intersection     | `figure2_deg_overview.R`                        |
| Fig. 3  | ICS scoring and 18-gene panel tier ranking               | `figure3_ics_panel.R`                           |
| Fig. 4  | C2C12 concordance scatter and module recovery            | `figure4_c2c12.R`                               |
| Fig. 5  | DGI network and Bliss synergy heatmap                    | `figure5_dgi_network.R`                         |

## Supplementary figures

| Figure | Description                                              |
|--------|----------------------------------------------------------|
| S1     | PRISMA flow diagram                                      |
| S2     | Per-dataset QC (mapping rate, library complexity)        |
| S3     | Bootstrap and permutation null distributions             |
| S4     | ICS sensitivity heatmap across weight configurations     |
| S5     | GO/KEGG enrichment dot plots                             |
| S6     | C2C12 in vitro concordance (extended)                    |
| S7     | Temporal trajectories with NLS fits                      |
| S8     | STRING community structure                               |

## Requirements

- R 4.2+, Bioconductor 3.16+
- `ggplot2`, `pheatmap`, `EnhancedVolcano`, `igraph`, `ggraph`, `cowplot`

Run all figures:

```bash
for s in figureS*/*.R figure*.R; do Rscript "$s"; done
```
