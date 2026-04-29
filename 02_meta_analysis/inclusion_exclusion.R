# =============================================================================
# inclusion_exclusion.R
# Apply PRISMA inclusion/exclusion criteria to candidate datasets
# =============================================================================
# Inclusion criteria:
#   - Spaceflight or simulated microgravity exposure (HU, RWV, clinostat)
#   - Skeletal muscle, myogenic cell line, or whole-body muscle-relevant tissue
#   - Mammalian (Mus musculus, Rattus norvegicus, Homo sapiens)
#   - Public raw or processed expression data (RNA-seq or microarray)
#   - n >= 3 per condition
# Exclusion criteria:
#   - Cardiac/smooth muscle only
#   - In vitro non-myogenic cell types
#   - No control group
#   - Withdrawn or superseded records
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Input: candidate dataset list from search strategy
candidates <- read_csv("data/raw/prisma_candidates.csv", show_col_types = FALSE)

# Apply criteria sequentially and log counts at each PRISMA stage
n0 <- nrow(candidates)

# Stage 1: deduplication
dedup <- candidates %>% distinct(accession, .keep_all = TRUE)
n1 <- nrow(dedup)

# Stage 2: title/abstract screening (manually curated flag)
screened <- dedup %>% filter(passes_title_abstract == TRUE)
n2 <- nrow(screened)

# Stage 3: full eligibility check
eligible <- screened %>%
  filter(
    exposure %in% c("spaceflight", "hindlimb_unloading", "RWV", "clinostat"),
    tissue_class %in% c("skeletal_muscle", "myogenic_cell"),
    organism %in% c("Mus musculus", "Rattus norvegicus", "Homo sapiens"),
    n_per_group >= 3,
    has_control == TRUE,
    !is_withdrawn
  )
n3 <- nrow(eligible)

# Stage 4: data availability confirmed (raw or processed downloadable)
included <- eligible %>% filter(data_available == TRUE)
n4 <- nrow(included)

# Write PRISMA flow numbers
prisma_flow <- tibble::tribble(
  ~stage,                          ~n,
  "Records identified",            n0,
  "After deduplication",           n1,
  "After title/abstract screen",   n2,
  "Eligible (full criteria)",      n3,
  "Included in meta-analysis",     n4
)
write_csv(prisma_flow, "data/processed/prisma_flow.csv")
write_csv(included,    "data/processed/included_datasets.csv")

cat("PRISMA flow saved. Final included datasets: n =", n4, "\n")
