# =============================================================================
# bliss_synergy.R
# Bliss-independence synergy scoring across drug pairs that target distinct
# Louvain communities. Bliss > 0.1 -> synergistic pair (manuscript threshold).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

dgi  <- read_csv("data/processed/dgi_interactions.csv",  show_col_types = FALSE)
comm <- read_csv("data/processed/string_communities.csv", show_col_types = FALSE)

# Map drug -> set of communities it covers (via target gene -> community)
drug_comm <- dgi %>%
  inner_join(comm, by = c("gene" = "symbol")) %>%
  group_by(drug) %>%
  summarise(communities = list(unique(community)),
            n_targets   = n_distinct(gene),
            approved    = any(approved), .groups = "drop")

# Pairwise Bliss: assume single-drug effects e_i (proxy = mean |interaction score|)
single_eff <- dgi %>%
  group_by(drug) %>%
  summarise(e = mean(abs(score), na.rm = TRUE), .groups = "drop")

drugs <- intersect(drug_comm$drug, single_eff$drug)
pairs <- expand.grid(d1 = drugs, d2 = drugs, stringsAsFactors = FALSE) %>%
  filter(d1 < d2) %>%
  inner_join(single_eff, by = c("d1" = "drug")) %>% rename(e1 = e) %>%
  inner_join(single_eff, by = c("d2" = "drug")) %>% rename(e2 = e) %>%
  mutate(bliss_expected = e1 + e2 - e1*e2)

# Observed combination effect placeholder: union of community coverage as proxy
cov <- setNames(drug_comm$communities, drug_comm$drug)
pairs <- pairs %>%
  mutate(coverage_union = mapply(function(a,b) length(union(cov[[a]], cov[[b]])),
                                 d1, d2),
         bliss_observed = pmin(1, coverage_union / max(coverage_union)),
         bliss_score    = bliss_observed - bliss_expected,
         synergistic    = bliss_score > 0.1)

write_csv(pairs, "data/processed/bliss_synergy_matrix.csv")
cat("Synergistic pairs (Bliss > 0.1):", sum(pairs$synergistic), "\n")
