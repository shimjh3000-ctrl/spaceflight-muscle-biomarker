# =============================================================================
# figure2_deg_overview.R
# Volcano plot of DESeq2 results + heatmap of 103-gene intersection.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(EnhancedVolcano); library(pheatmap)
  library(readr); library(dplyr)
})

deg <- read_csv("data/processed/deseq2_results.csv", show_col_types = FALSE)
inter <- read_csv("data/processed/intersection_103.csv", show_col_types = FALSE)
expr <- as.matrix(read.csv("data/processed/vsd_expression.csv",
                           row.names = 1, check.names = FALSE))

# Panel A — volcano
p1 <- EnhancedVolcano(deg,
                      lab = deg$symbol,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 0.5,
                      title = "Spaceflight vs Ground Control",
                      subtitle = "DESeq2 (n = 8 + 8)")
ggsave("figures/figure2A_volcano.pdf", p1, width = 7, height = 6)

# Panel B — 103-gene heatmap
mat <- expr[rownames(expr) %in% inter$symbol, ]
ann <- data.frame(condition = ifelse(grepl("SF", colnames(mat)), "SF", "GC"),
                  row.names = colnames(mat))
pheatmap(mat, scale = "row",
         show_rownames = FALSE, annotation_col = ann,
         filename = "figures/figure2B_heatmap.pdf",
         width = 6, height = 8)
