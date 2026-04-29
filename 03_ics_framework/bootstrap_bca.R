#!/usr/bin/env Rscript
# BCa bootstrap confidence intervals for ICS scores (n = 10,000 iterations)
#
# Manuscript values to reproduce:
#   TNNT3 9.4 (95% CI 9.1–9.7); TRIM63 9.2 (8.9–9.5);
#   ACTC1 8.7 (8.4–9.0); IGF2 8.5 (8.2–8.8)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(boot)
})

set.seed(2026)
N_BOOT <- 10000

df <- read_csv("data/intersection_103genes.csv", show_col_types = FALSE)

minmax10 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(5, length(x)))
  10 * (x - rng[1]) / (rng[2] - rng[1])
}

ics_for <- function(gene_row, ref) {
  log2fc_n <- minmax10(c(ref$log2fc_abs, abs(gene_row$log2fc)))[length(ref$log2fc_abs) + 1]
  fdr_n    <- minmax10(c(ref$neg_log10_fdr, -log10(gene_row$fdr)))[length(ref$neg_log10_fdr) + 1]
  ns_n     <- minmax10(c(ref$n_studies, gene_row$n_studies))[length(ref$n_studies) + 1]
  dg_n     <- minmax10(c(ref$disgenet_score, gene_row$disgenet_score))[length(ref$disgenet_score) + 1]
  cat_n    <- minmax10(c(ref$category_weight, gene_row$category_weight))[length(ref$category_weight) + 1]
  S4522 <- 0.6 * log2fc_n + 0.4 * fdr_n
  S726  <- 0.5 * ns_n + 0.25 * dg_n + 0.25 * cat_n
  0.6 * S4522 + 0.4 * S726
}

ref <- df %>% transmute(log2fc_abs = abs(log2fc),
                        neg_log10_fdr = -log10(fdr),
                        n_studies, disgenet_score, category_weight)

ics_stat <- function(data, idx) {
  resampled <- data[idx, ]
  vapply(seq_len(nrow(df)), function(i) ics_for(df[i, ], resampled), numeric(1))[seq_len(nrow(df))]
}

# Per-gene bootstrap (BCa)
result <- list()
for (i in seq_len(nrow(df))) {
  gene <- df$gene[i]
  b <- boot(ref, statistic = function(d, idx) ics_for(df[i, ], d[idx, ]), R = N_BOOT)
  ci <- tryCatch(boot.ci(b, type = "bca", conf = 0.95), error = function(e) NULL)
  if (!is.null(ci)) {
    result[[gene]] <- tibble(gene = gene,
                             ICS_point = ci$t0,
                             CI_low    = ci$bca[4],
                             CI_high   = ci$bca[5])
  }
}
out <- bind_rows(result)
write_csv(out, "data/bootstrap_ci.csv")
cat(sprintf("[bootstrap] %d genes × %d iterations done.\n", nrow(out), N_BOOT))
