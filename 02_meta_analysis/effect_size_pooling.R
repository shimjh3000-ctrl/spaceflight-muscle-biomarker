#!/usr/bin/env Rscript
# Random-effects meta-analysis (REML) across 412 spaceflight studies
# Output: 726-gene consensus set with pooled effect sizes
#
# Usage: Rscript effect_size_pooling.R <input_long_table> <output_csv>

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
in_path  <- args[[1]]   # long-format: gene, study, log2fc, se, n
out_path <- args[[2]]

# Long-format table extracted from the 412 included studies
df <- read_csv(in_path, show_col_types = FALSE)

# Pool per gene
pool_one <- function(g) {
  gd <- df %>% filter(gene == g)
  if (nrow(gd) < 3) return(NULL)
  fit <- tryCatch(
    rma(yi = gd$log2fc, sei = gd$se, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  tibble(
    gene = g,
    n_studies = nrow(gd),
    pooled_log2fc = as.numeric(fit$beta),
    pooled_se = fit$se,
    p_value = fit$pval,
    I2 = fit$I2,
    tau2 = fit$tau2,
    egger_p = tryCatch(regtest(fit, model = "lm")$pval, error = function(e) NA_real_)
  )
}

genes <- unique(df$gene)
out <- bind_rows(lapply(genes, pool_one))
out <- out %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  filter(I2 < 75, n_studies >= 3, abs(pooled_log2fc) > 0.5, fdr < 0.05)

write_csv(out, out_path)
cat(sprintf("[meta-analysis] 726-gene consensus expected; obtained %d genes.\n", nrow(out)))
