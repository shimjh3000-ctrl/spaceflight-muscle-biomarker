#!/usr/bin/env python3
"""string_louvain.py
Build STRING v12 high-confidence subnetwork around the 18-gene panel and
detect communities by Louvain modularity. Saves community assignments to
data/processed/string_communities.csv.
"""
from __future__ import annotations

import io
import sys
import pandas as pd
import requests
import networkx as nx
import community as community_louvain  # python-louvain

PANEL_CSV = "data/processed/panel_18.csv"
STRING_VERSION = "12.0"
SPECIES = 9606  # human; rerun with 10090 for mouse if needed
SCORE_THRESHOLD = 700  # 0.7 confidence

panel = pd.read_csv(PANEL_CSV)
genes = panel["symbol"].tolist()

# 1) STRING API — get interactions
url = f"https://string-db.org/api/tsv/network"
params = {
    "identifiers": "%0d".join(genes),
    "species": SPECIES,
    "required_score": SCORE_THRESHOLD,
    "network_type": "physical",
    "caller_identity": "spaceflight-muscle-biomarker",
}
r = requests.post(url, data=params, timeout=60)
r.raise_for_status()
edges = pd.read_csv(io.StringIO(r.text), sep="\t")

# 2) Build graph
G = nx.Graph()
for _, row in edges.iterrows():
    G.add_edge(row["preferredName_A"], row["preferredName_B"],
               weight=float(row["score"]))

# 3) Louvain community detection
partition = community_louvain.best_partition(G, weight="weight", random_state=42)
mod = community_louvain.modularity(partition, G, weight="weight")

out = pd.DataFrame({
    "symbol": list(partition.keys()),
    "community": [f"C{v}" for v in partition.values()],
})
out.to_csv("data/processed/string_communities.csv", index=False)

print(f"Louvain modularity: {mod:.3f}")
print(out["community"].value_counts())
