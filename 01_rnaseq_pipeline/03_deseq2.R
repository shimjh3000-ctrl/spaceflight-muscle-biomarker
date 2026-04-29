#!/usr/bin/env Rscript
# DESeq2 v1.34.0 differential expression: ISS THP-1 (n=8 SF) vs ground control (n=8 GC)
# FDR < 0.05, |Log2FC| > 0.5; apeglm LFC shrinkage
#
# Usage: Rscript 03_deseq2.R <counts_dir> <output_dir>

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
counts_dir <- args[[1]]
out_dir    <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load count matrix ------------------------------------------------
counts_path <- file.path(counts_dir, "gene_counts.tsv")
raw <- read.table(counts_path, header = TRUE, sep = "\t",
                  comment.char = "#", check.names = FALSE)
cnt <- as.matrix(raw[, -(1:6)])
rownames(cnt) <- raw$Geneid
colnames(cnt) <- sub("_Aligned.*", "", basename(colnames(cnt)))

# Drop SF08b technical replicate from statistical analysis (kept for vis only)
cnt <- cnt[, !grepl("SF08b", colnames(cnt))]

# ---- Sample metadata --------------------------------------------------
condition <- ifelse(grepl("^SF", colnames(cnt)), "spaceflight", "ground_control")
coldata <- data.frame(row.names = colnames(cnt),
                      condition = factor(condition, levels = c("ground_control", "spaceflight")))

stopifnot(table(coldata$condition)[["spaceflight"]] == 8,
          table(coldata$condition)[["ground_control"]] == 8)

# ---- DESeq2 -----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = coldata,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# apeglm LFC shrinkage
res <- lfcShrink(dds, coef = "condition_spaceflight_vs_ground_control",
                 type = "apeglm")

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

# ---- Filter to DEGs (FDR < 0.05, |Log2FC| > 0.5) ----------------------
degs <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.5)

write_csv(res_df, file.path(out_dir, "THP1_SF_vs_GC_full.csv"))
write_csv(degs,  file.path(out_dir, "THP1_SF_vs_GC.csv"))

cat(sprintf("[DESeq2] %d DEGs at FDR<0.05 and |Log2FC|>0.5\n", nrow(degs)))
cat(sprintf("[DESeq2] Expected from manuscript: 4,522\n"))
