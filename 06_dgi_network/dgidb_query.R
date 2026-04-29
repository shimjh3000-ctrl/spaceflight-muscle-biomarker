# =============================================================================
# dgidb_query.R
# Query DGIdb v5.0.6 GraphQL API for drug-gene interactions of panel genes.
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(readr)
})

panel <- read_csv("data/processed/panel_18.csv", show_col_types = FALSE)

dgidb_url <- "https://dgidb.org/api/graphql"

query_template <- '
{
  genes(names: [%s]) {
    nodes {
      name
      interactions {
        drug { name conceptId approved }
        interactionScore
        interactionTypes { type directionality }
        sources { sourceDbName }
      }
    }
  }
}'

gene_list <- paste(sprintf("\"%s\"", panel$symbol), collapse = ", ")
qbody <- list(query = sprintf(query_template, gene_list))

resp <- POST(dgidb_url,
             body = toJSON(qbody, auto_unbox = TRUE),
             content_type_json())
stop_for_status(resp)
parsed <- content(resp, as = "parsed", simplifyVector = FALSE)

rows <- list()
for (g in parsed$data$genes$nodes) {
  for (i in g$interactions) {
    rows[[length(rows)+1]] <- tibble(
      gene  = g$name,
      drug  = i$drug$name,
      approved = i$drug$approved %||% FALSE,
      score = i$interactionScore %||% NA_real_,
      type  = paste(sapply(i$interactionTypes, `[[`, "type"), collapse = ";")
    )
  }
}
dgi <- bind_rows(rows)

write_csv(dgi, "data/processed/dgi_interactions.csv")
cat("DGIdb interactions saved: n =", nrow(dgi), "\n")
