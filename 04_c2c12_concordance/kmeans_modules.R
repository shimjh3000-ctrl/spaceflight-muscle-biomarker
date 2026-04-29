# =============================================================================
# kmeans_modules.R
# k-means co-expression modules on C2C12 RWV variance-stabilised expression.
# Tests whether 4 modules (atrogenes / contractile / IGF-AKT / ROS) emerge.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(readr)
  library(dplyr)
})

set.seed(42)

# Load DESeq2 object saved in reanalyze_pipeline.sh (re-build vsd here)
cts <- as.matrix(read.table("data/processed/c2c12/counts/c2c12_counts.tsv",
                            header=TRUE, row.names=1, skip=1, check.names=FALSE))
cts <- cts[, grep("Aligned", colnames(cts))]
colnames(cts) <- sub("_Aligned.*", "", basename(colnames(cts)))
coldata <- data.frame(row.names = colnames(cts),
                      condition = factor(ifelse(grepl("RWV", colnames(cts)),
                                                "RWV", "static")))
dds  <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
vsd  <- vst(dds, blind = FALSE)
expr <- assay(vsd)

# Use top 2000 most variable genes for clustering
vars <- rowVars(expr)
top  <- order(vars, decreasing = TRUE)[1:2000]
expr_top <- t(scale(t(expr[top, ])))

# k-means with k = 4 (4 hypothesised modules)
km <- kmeans(expr_top, centers = 4, nstart = 50, iter.max = 100)

modules <- tibble(
  gene_id = rownames(expr_top),
  module  = paste0("M", km$cluster)
)
write_csv(modules, "data/processed/c2c12_kmeans_modules.csv")

# Module composition summary
cat("Module sizes:\n"); print(table(modules$module))
