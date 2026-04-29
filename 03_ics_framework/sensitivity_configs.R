#!/usr/bin/env Rscript
# Sensitivity analysis across four alternative ICS weighting configurations (A–D).
# Computes Robustness Index = fraction of tier-stable genes vs Primary.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

df <- read_csv("data/intersection_103genes.csv", show_col_types = FALSE)

minmax10 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(5, length(x)))
  10 * (x - rng[1]) / (rng[2] - rng[1])
}

base <- df %>% mutate(
  l = minmax10(abs(log2fc)),
  f = minmax10(-log10(fdr)),
  s = minmax10(n_studies),
  d = minmax10(disgenet_score),
  c = minmax10(category_weight)
)

configs <- list(
  Primary = c(w_S4522 = 0.6, w_S726 = 0.4, w_l = 0.6, w_f = 0.4, w_s = 0.5,  w_d = 0.25, w_c = 0.25),
  A       = c(w_S4522 = 0.5, w_S726 = 0.5, w_l = 0.6, w_f = 0.4, w_s = 0.5,  w_d = 0.25, w_c = 0.25),
  B       = c(w_S4522 = 0.7, w_S726 = 0.3, w_l = 0.6, w_f = 0.4, w_s = 0.5,  w_d = 0.25, w_c = 0.25),
  C       = c(w_S4522 = 0.4, w_S726 = 0.6, w_l = 0.6, w_f = 0.4, w_s = 0.5,  w_d = 0.25, w_c = 0.25),
  D       = c(w_S4522 = 0.6, w_S726 = 0.4, w_l = 0.6, w_f = 0.4, w_s = 0.4,  w_d = 0.4,  w_c = 0.2)
)

tier_of <- function(ics) ifelse(ics >= 8, "Gold",
                          ifelse(ics >= 7, "High",
                            ifelse(ics >= 6, "Moderate", "Low")))

apply_cfg <- function(cfg) {
  with(cfg, {
    s4522 <- w_l * base$l + w_f * base$f
    s726  <- w_s * base$s + w_d * base$d + w_c * base$c
    ics   <- w_S4522 * s4522 + w_S726 * s726
    tier_of(ics)
  })
}

mat <- sapply(configs, apply_cfg)
colnames(mat) <- names(configs)
out <- bind_cols(gene = base$gene, as.data.frame(mat))

# Robustness Index
ri <- mean(apply(mat, 1, function(r) length(unique(r)) == 1))
gold_stable <- sum(apply(mat[mat[, "Primary"] == "Gold", , drop = FALSE], 1,
                         function(r) all(r == "Gold")))
gold_total <- sum(mat[, "Primary"] == "Gold")

write_csv(out, "data/sensitivity_robustness.csv")
cat(sprintf("[sensitivity] Robustness Index (overall): %.3f\n", ri))
cat(sprintf("[sensitivity] Gold-tier stable across all 5 configs: %d / %d\n",
            gold_stable, gold_total))
cat(sprintf("[sensitivity] Expected: Gold-tier Robustness Index = 1.00\n"))
