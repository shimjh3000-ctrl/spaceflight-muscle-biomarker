# =============================================================================
# pearson_correlation.R
# Pearson correlation between in vivo (GSE245789/OSD-590 meta) and
# in vitro C2C12 RWV log2 fold changes for the 18-gene panel.
# Manuscript reports r = 0.78 (P < 0.001).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

panel <- read_csv("data/processed/panel_18.csv", show_col_types = FALSE)
invivo <- read_csv("data/processed/meta_pooled_log2fc.csv", show_col_types = FALSE) %>%
  select(symbol, log2fc_invivo = pooled_log2fc)
invitro <- read_csv("data/processed/c2c12_log2fc.csv", show_col_types = FALSE) %>%
  rename(symbol = 1, log2fc_invitro = log2FoldChange) %>%
  select(symbol, log2fc_invitro)

merged <- panel %>%
  inner_join(invivo,  by = "symbol") %>%
  inner_join(invitro, by = "symbol")

ct <- cor.test(merged$log2fc_invivo, merged$log2fc_invitro, method = "pearson")
cat(sprintf("Pearson r = %.3f (95%% CI %.3f-%.3f), P = %.2e, n = %d\n",
            ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value, nrow(merged)))

write_csv(merged, "data/processed/c2c12_panel_concordance.csv")

p <- ggplot(merged, aes(log2fc_invivo, log2fc_invitro, label = symbol)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE) +
  ggrepel::geom_text_repel(size = 3) +
  labs(x = "In vivo log2FC (pooled)", y = "In vitro C2C12 RWV log2FC",
       title = sprintf("18-gene panel concordance (r = %.2f)", ct$estimate)) +
  theme_classic()
ggsave("figures/figureS6_c2c12_concordance.pdf", p, width = 6, height = 5)
