#!/usr/bin/env Rscript
# Permutation null distribution for ICS scores (n = 1,000 shuffled gene-disease
# annotations). Returns empirical p-values per gene.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

set.seed(2026)
N_PERM <- 1000

df <- read_csv("data/intersection_103genes.csv", show_col_types = FALSE)

minmax10 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(5, length(x)))
  10 * (x - rng[1]) / (rng[2] - rng[1])
}

compute_ics <- function(d) {
  d %>% mutate(
    log2fc_norm     = minmax10(abs(log2fc)),
    neg_log10_fdr_n = minmax10(-log10(fdr)),
    n_studies_norm  = minmax10(n_studies),
    disgenet_norm   = minmax10(disgenet_score),
    category_norm   = minmax10(category_weight),
    S4522 = 0.6 * log2fc_norm + 0.4 * neg_log10_fdr_n,
    S726  = 0.5 * n_studies_norm + 0.25 * disgenet_norm + 0.25 * category_norm,
    ICS   = 0.6 * S4522 + 0.4 * S726
  )
}

obs <- compute_ics(df) %>% select(gene, ICS_obs = ICS)

null_dist <- replicate(N_PERM, {
  shuf <- df
  shuf$disgenet_score  <- sample(shuf$disgenet_score)
  shuf$category_weight <- sample(shuf$category_weight)
  compute_ics(shuf)$ICS
})  # matrix [genes × N_PERM]

emp_p <- vapply(seq_len(nrow(df)), function(i) {
  mean(null_dist[i, ] >= obs$ICS_obs[i])
}, numeric(1))

result <- obs %>% mutate(empirical_p = emp_p)
write_csv(result, "data/permutation_p_values.csv")

cat(sprintf("[permutation] All Gold/High-tier genes: empirical p < 0.001 expected.\n"))
cat(sprintf("[permutation] Genes with p < 0.001: %d\n", sum(result$empirical_p < 0.001)))
