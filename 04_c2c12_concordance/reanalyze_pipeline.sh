#!/usr/bin/env bash
# =============================================================================
# reanalyze_pipeline.sh
# End-to-end reanalysis of C2C12 RWV simulated-microgravity RNA-seq
# Same parameters as 01_rnaseq_pipeline (STAR 2.7.10a, featureCounts 2.0.3,
# DESeq2 1.34.0) for direct comparability.
# =============================================================================
set -euo pipefail

FASTQ_DIR="data/raw/c2c12_rwv"
OUT_DIR="data/processed/c2c12"
GENOME_DIR="reference/STAR_GRCm39"
GTF="reference/gencode.vM30.primary_assembly.annotation.gtf"
THREADS=8

mkdir -p "${OUT_DIR}/bam" "${OUT_DIR}/counts"

# 1) STAR alignment
for fq1 in "${FASTQ_DIR}"/*_R1.fastq.gz; do
  sample=$(basename "${fq1}" _R1.fastq.gz)
  fq2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
  STAR --runThreadN ${THREADS} \
       --genomeDir "${GENOME_DIR}" \
       --readFilesIn "${fq1}" "${fq2}" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix "${OUT_DIR}/bam/${sample}_"
done

# 2) featureCounts (paired-end, stranded reverse)
featureCounts -T ${THREADS} -p --countReadPairs -s 2 \
              -a "${GTF}" \
              -o "${OUT_DIR}/counts/c2c12_counts.tsv" \
              "${OUT_DIR}"/bam/*_Aligned.sortedByCoord.out.bam

# 3) DESeq2 (RWV vs static control)
Rscript - <<'RS'
suppressPackageStartupMessages({ library(DESeq2); library(readr) })
cts <- as.matrix(read.table("data/processed/c2c12/counts/c2c12_counts.tsv",
                            header=TRUE, row.names=1, skip=1, check.names=FALSE))
cts <- cts[, grep("Aligned", colnames(cts))]
colnames(cts) <- sub("_Aligned.*", "", basename(colnames(cts)))
coldata <- data.frame(row.names = colnames(cts),
                      condition = factor(ifelse(grepl("RWV", colnames(cts)),
                                                "RWV", "static"),
                                         levels = c("static", "RWV")))
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "RWV", "static"))
write.csv(as.data.frame(res), "data/processed/c2c12_log2fc.csv")
RS

echo "C2C12 reanalysis pipeline complete."
