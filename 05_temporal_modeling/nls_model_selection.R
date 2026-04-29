# =============================================================================
# nls_model_selection.R
# Fit competing kinetic models to each panel gene's timecourse and select by AICc
# Models:
#   M1 exponential : y = a * (1 - exp(-k*t)) + b
#   M2 logistic    : y = L / (1 + exp(-k*(t - t0))) + b
#   M3 biphasic    : y = a1*(1-exp(-k1*t)) + a2*(1-exp(-k2*t)) + b
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(MuMIn); library(minpack.lm)
})

tc    <- read_csv("data/processed/timecourse_normalised.csv", show_col_types = FALSE)
panel <- read_csv("data/processed/panel_18.csv",              show_col_types = FALSE)

fit_models <- function(df) {
  out <- list()
  out$M1 <- try(nlsLM(value ~ a*(1-exp(-k*time)) + b, data = df,
                      start = list(a = 1, k = 0.05, b = 0)), silent = TRUE)
  out$M2 <- try(nlsLM(value ~ L/(1+exp(-k*(time-t0))) + b, data = df,
                      start = list(L = 1, k = 0.05, t0 = 48, b = 0)), silent = TRUE)
  out$M3 <- try(nlsLM(value ~ a1*(1-exp(-k1*time)) + a2*(1-exp(-k2*time)) + b,
                      data = df,
                      start = list(a1 = 0.5, k1 = 0.1, a2 = 0.5, k2 = 0.01, b = 0)),
                silent = TRUE)
  out
}

results <- panel$symbol |> lapply(function(g) {
  df <- tc |> filter(symbol == g) |> select(time, value)
  if (nrow(df) < 4) return(NULL)
  fits <- fit_models(df)
  aicc <- sapply(fits, function(f) if (inherits(f, "try-error")) NA else AICc(f))
  best <- names(which.min(aicc))
  tibble(symbol = g, best_model = best,
         aicc_M1 = aicc["M1"], aicc_M2 = aicc["M2"], aicc_M3 = aicc["M3"])
}) |> bind_rows()

write_csv(results, "data/processed/nls_fits.csv")
cat("NLS fits saved. Model selection summary:\n"); print(table(results$best_model))
