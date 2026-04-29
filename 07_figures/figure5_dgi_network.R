# =============================================================================
# figure5_dgi_network.R
# Render the DGI + STRING-Louvain network coloured by community, with edge
# weight proportional to Bliss synergy score.
# =============================================================================

suppressPackageStartupMessages({
  library(igraph); library(ggraph); library(ggplot2)
  library(readr); library(dplyr); library(tidyr)
})

dgi   <- read_csv("data/processed/dgi_interactions.csv",   show_col_types = FALSE)
comm  <- read_csv("data/processed/string_communities.csv", show_col_types = FALSE)
bliss <- read_csv("data/processed/bliss_synergy_matrix.csv", show_col_types = FALSE)

# Build bipartite gene-drug graph
edges <- dgi %>%
  transmute(from = gene, to = drug, weight = abs(score))
nodes <- bind_rows(
  tibble(name = unique(edges$from), type = "gene") %>%
    left_join(comm, by = c("name" = "symbol")),
  tibble(name = unique(edges$to),   type = "drug", community = NA_character_)
)
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.4, colour = "grey60") +
  geom_node_point(aes(colour = community, shape = type), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_width(range = c(0.3, 1.5)) +
  theme_void() +
  labs(title = "DGI + STRING network around 18-gene panel")

ggsave("figures/figure5_dgi_network.pdf", p, width = 9, height = 7)
