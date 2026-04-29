#!/usr/bin/env bash
# STAR v2.7.10a paired-end alignment to GRCh38.p13 / GENCODE v38
# Usage: bash 01_star_align.sh <fastq_dir> <star_index> <output_dir>
set -euo pipefail

FASTQ_DIR="${1:?fastq_dir required}"
STAR_INDEX="${2:?star_index required}"
OUT_DIR="${3:?output_dir required}"
THREADS="${THREADS:-16}"

mkdir -p "${OUT_DIR}"

for R1 in "${FASTQ_DIR}"/*_R1*.fastq.gz; do
    SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
    R2="${R1/_R1/_R2}"
    echo "[STAR] ${SAMPLE}"

    STAR \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${R1}" "${R2}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --quantMode GeneCounts
done

echo "[STAR] Alignment complete. Outputs in: ${OUT_DIR}"
