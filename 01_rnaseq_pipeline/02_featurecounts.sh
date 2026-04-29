#!/usr/bin/env bash
# featureCounts v2.0.1 stranded paired-end gene-level quantification
# Usage: bash 02_featurecounts.sh <bam_dir> <gtf> <output_dir>
set -euo pipefail

BAM_DIR="${1:?bam_dir required}"
GTF="${2:?gtf required}"
OUT_DIR="${3:?output_dir required}"
THREADS="${THREADS:-8}"

mkdir -p "${OUT_DIR}"

featureCounts \
  -T "${THREADS}" \
  -p --countReadPairs \
  -s 2 \
  -t exon \
  -g gene_id \
  -a "${GTF}" \
  -o "${OUT_DIR}/gene_counts.tsv" \
  "${BAM_DIR}"/*_Aligned.sortedByCoord.out.bam

echo "[featureCounts] Done: ${OUT_DIR}/gene_counts.tsv"
