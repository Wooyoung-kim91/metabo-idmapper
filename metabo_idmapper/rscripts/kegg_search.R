#!/usr/bin/env Rscript
# kegg_search — KEGGREST keggFind("compound", term) candidate search (verify later!).
# stdin  JSON: {"terms":["N-methylarginine","monomethylarginine"], "max_per_term":8}
# stdout JSON: [{"term":..,"kegg":"C00123","kegg_name":"..."}]  (flat candidate list)
suppressMessages({library(KEGGREST); library(jsonlite)})

req <- fromJSON(paste(readLines("stdin", warn = FALSE), collapse = "\n"),
                simplifyVector = FALSE)
terms <- unlist(req$terms)
maxn <- if (is.null(req$max_per_term)) 8L else as.integer(req$max_per_term)

rows <- list()
for (term in terms) {
  r <- tryCatch(keggFind("compound", term), error = function(e) character(0))
  if (length(r) > 0) {
    n <- min(length(r), maxn)
    for (k in seq_len(n)) {
      rows[[length(rows) + 1]] <- list(
        term = term,
        kegg = sub("^cpd:", "", names(r)[k]),
        kegg_name = unname(r[k])
      )
    }
  }
}
writeLines(toJSON(rows, auto_unbox = TRUE, null = "null"), req$out)
