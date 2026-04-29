# =============================================================================
# two_tier_protocol.R
# Assign panel genes to Early (Tier 1) or Late (Tier 2) tiers based on t1/2
# (time to half-maximal change) derived from best-fit NLS model.
# Threshold: t1/2 <= 24 h -> Tier 1; t1/2 >= 72 h -> Tier 2; else mixed
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })

fits  <- read_csv("data/processed/nls_fits.csv", show_col_types = FALSE)
panel <- read_csv("data/processed/panel_18.csv", show_col_types = FALSE)

# Half-times (precomputed from best-fit parameters; see nls_model_selection.R)
halftimes <- read_csv("data/processed/nls_halftimes.csv", show_col_types = FALSE)

assigned <- panel %>%
  inner_join(halftimes, by = "symbol") %>%
  mutate(tier = case_when(
    t_half_h <= 24 ~ "Tier1_Early",
    t_half_h >= 72 ~ "Tier2_Late",
    TRUE           ~ "Mixed"
  ))

write_csv(assigned, "data/processed/two_tier_assignment.csv")

cat("Two-tier assignment summary:\n"); print(table(assigned$tier))
