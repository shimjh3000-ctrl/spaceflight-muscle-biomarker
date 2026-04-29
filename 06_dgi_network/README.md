# 06 — Drug-Gene Interaction Network and Combinatorial Synergy

Bridge from the 18-gene panel to actionable countermeasure candidates by
integrating DGIdb interactions, STRING-Louvain network communities, and
Bliss-independence synergy scoring.

## Goal

Identify drugs / drug classes whose combined target footprint covers
multiple network communities and shows positive synergy (Bliss > 0.1).

## Inputs

- `data/processed/panel_18.csv`
- DGIdb v5.0.6 API (queried at runtime)
- STRING v12 protein interactions (high-confidence, score ≥ 0.7)

## Outputs

- `data/processed/dgi_interactions.csv`
- `data/processed/string_communities.csv` (Louvain modularity)
- `data/processed/bliss_synergy_matrix.csv`
- `figures/figure5_dgi_network.pdf`

## Run order

```bash
Rscript dgidb_query.R                 # DGIdb v5.0.6 -> drug-gene table
python string_louvain.py              # STRING + Louvain communities
Rscript bliss_synergy.R               # Bliss independence over community pairs
```
