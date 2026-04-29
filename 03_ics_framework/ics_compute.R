#!/usr/bin/env Rscript
# Integration Confidence Score (ICS) computation
# ICS = 0.6 × S4522 + 0.4 × S726 (all sub-scores min-max normalised to [0,10])
#
# Inputs:
#   data/intersection_103genes.csv  — must contain columns:
#     gene, log2fc, fdr, n_studies, disgenet_score, category_weight
#
# Outputs:
#   data/ics_scores.csv      (all 103 intersection genes ranked)
#   data/18_gene_panel.csv   (Gold + High tier genes)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

input  <- "data/intersection_103genes.csv"
out_a  <- "data/ics_scores.csv"
out_b  <- "data/18_gene_panel.csv"

if (!file.exists(input)) {
  stop("Place intersection_103genes.csv in the data/ directory before running.\n",
       "Source columns required: gene, log2fc, fdr, n_studies, disgenet_score, category_weight")
}

minmax10 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(5, length(x)))
  10 * (x - rng[1]) / (rng[2] - rng[1])
}

df <- read_csv(input, show_col_types = FALSE) %>%
  mutate(
    log2fc_abs       = abs(log2fc),
    neg_log10_fdr    = -log10(fdr),
    log2fc_norm      = minmax10(log2fc_abs),
    neg_log10_fdr_n  = minmax10(neg_log10_fdr),
    n_studies_norm   = minmax10(n_studies),
    disgenet_norm    = minmax10(disgenet_score),
    category_norm    = minmax10(category_weight),
    S4522            = 0.6 * log2fc_norm + 0.4 * neg_log10_fdr_n,
    S726             = 0.5 * n_studies_norm + 0.25 * disgenet_norm + 0.25 * category_norm,
    ICS              = 0.6 * S4522 + 0.4 * S726,
    tier = case_when(
      ICS >= 8.0 ~ "Gold",
      ICS >= 7.0 ~ "High",
      ICS >= 6.0 ~ "Moderate",
      TRUE       ~ "Low"
    )
  ) %>%
  arrange(desc(ICS))

write_csv(df, out_a)

panel18 <- df %>% filter(tier %in% c("Gold", "High"))
write_csv(panel18, out_b)

cat(sprintf("[ICS] %d intersection genes scored.\n", nrow(df)))
cat(sprintf("[ICS] Gold tier: %d | High tier: %d | Moderate tier: %d\n",
            sum(df$tier == "Gold"), sum(df$tier == "High"), sum(df$tier == "Moderate")))
cat(sprintf("[ICS] 18-gene panel written to: %s\n", out_b))
