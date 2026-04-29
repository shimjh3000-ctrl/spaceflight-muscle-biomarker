# 01 — RNA-seq Pipeline (ISS THP-1)

Processing pipeline for the primary discovery dataset: **GSE245789 / OSD-590** (n = 8 spaceflight + 8 ground control THP-1 samples; 16 biological + 1 technical replicate = 17 lanes).

## Workflow

```
FASTQ → trim-galore → STAR (GRCh38.p13/GENCODE v38) → featureCounts → DESeq2 → DEGs
```

## Scripts

| Script | Purpose |
|---|---|
| `01_star_align.sh`   | STAR v2.7.10a paired-end alignment, 150 bp reads |
| `02_featurecounts.sh` | featureCounts v2.0.1 gene-level quantification (stranded) |
| `03_deseq2.R`        | DESeq2 v1.34.0 differential expression with apeglm log2FC shrinkage |

## Key parameters

- Alignment: STAR `--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFilterType BySJout`
- Quantification: featureCounts `-s 2 -p --countReadPairs -t exon -g gene_id`
- DE thresholds: **FDR < 0.05**, **|Log2FC| > 0.5**
- LFC shrinkage: apeglm
- Sample exclusion: SF08b (technical replicate of SF08; Pearson r = 0.997) included for visualisation only; statistical analyses use n = 8 biological replicates per group.

## Quality control acceptance criteria

| Metric | Threshold | All 16 biological samples |
|---|---|---|
| Mapped reads | ≥ 25.8 M | ✓ |
| Mapping rate | 92.6–96.9% | ✓ |
| RIN | ≥ 8.5 | ✓ (≥ 8.8 observed) |

See **Supplementary Table S12** for per-sample QC metrics.

## Output

`THP1_SF_vs_GC.csv` — 4,522 DEGs (Set S4522 in the manuscript).
