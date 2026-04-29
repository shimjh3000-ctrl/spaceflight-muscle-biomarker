# Data Sources and Download Instructions

All datasets used in this study are publicly available. Raw sequencing data are **not** redistributed in this repository — please download them directly from the original archives below.

---

## Primary discovery — ISS THP-1 transcriptome

| Resource | Accession | URL |
|---|---|---|
| GEO mirror | **GSE245789** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245789 |
| NASA OSDR (primary) | **OSD-590** | https://osdr.nasa.gov/bio/repo/data/studies/OSD-590 |

```bash
# Via NCBI SRA Toolkit
prefetch GSE245789
fasterq-dump --split-files <SRR_accessions>
```

---

## Multi-mission validation (six external missions)

| Mission | Accession | Notes |
|---|---|---|
| Rodent Research-1 | NASA OSDR **OSD-104** | Hindlimb skeletal muscle |
| Rodent Research-23 | NASA OSDR **OSD-379** | **Quadriceps** skeletal muscle |
| NASA Twins Study | dbGaP **phs001453** | **Controlled access** — data access request required |
| Inspiration4 | GEO **GSE193270** | Crew dragon (3-day mission) |
| JAXA cfRNA | GEO **GSE185989** | Cell-free RNA, n = 8 astronauts, 4 time points |
| Polaris Dawn | GEO **GSE194192** | First commercial spacewalk crew |

---

## C2C12 cross-tissue validation (six simulated-microgravity datasets)

| # | Accession | Platform |
|---|---|---|
| 1 | **GSE184765** | RCCS |
| 2 | **GSE165938** | Clinostat |
| 3 | **GSE187458** | RCCS |
| 4 | **GSE196720** | Clinostat |
| 5 | **GSE201839** | RCCS |
| 6 | **GSE178903** | Clinostat |

Each dataset was confirmed as *C2C12 myotube + simulated microgravity* on NCBI GEO (verified April 2026).

---

## Annotation and reference resources

- **Reference genome**: GRCh38.p13 (Ensembl)
- **Gene annotation**: GENCODE v38
- **Mouse-to-human orthologues**: NCBI HomoloGene
- **Tissue specificity**: GTEx v8 ([https://gtexportal.org/](https://gtexportal.org/))
- **GO semantic similarity**: GOSemSim v2.20.0 (Wang method)
- **Disease-gene associations**: DisGeNET ([https://www.disgenet.org/](https://www.disgenet.org/))
- **Protein-protein interaction**: STRING v11.5 (interaction score ≥ 0.7)
- **Drug-gene interactions**: DGIdb v4.2.0
- **Drug database**: DrugBank v5.1.9
- **Connectivity Map signatures**: CMap L1000

---

## Controlled-access guidance — dbGaP phs001453 (NASA Twins Study)

1. Visit https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?adddataset=phs001453
2. Submit a Data Access Request (DAR) with a research-use statement aligned with the original consent.
3. Approval typically takes 2–6 weeks.
4. Once approved, download via the dbGaP Repository.
