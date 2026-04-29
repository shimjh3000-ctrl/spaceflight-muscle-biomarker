# =============================================================================
# figure3_ics_panel.R
# ICS bar plot for the 18-gene panel with tier annotation (Gold/Silver/Bronze).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(dplyr)
})

panel <- read_csv("data/processed/panel_18.csv", show_col_types = FALSE) %>%
  arrange(desc(ics))

panel$symbol <- factor(panel$symbol, levels = panel$symbol)

p <- ggplot(panel, aes(symbol, ics, fill = tier)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f", ics)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c(Gold = "#D4A017",
                               Silver = "#A8A8A8",
                               Bronze = "#B87333")) +
  labs(x = NULL, y = "Integrated Concordance Score",
       title = "18-gene blood-accessible biomarker panel") +
  theme_classic()

ggsave("figures/figure3_ics_panel.pdf", p, width = 7, height = 6)
