#!/usr/bin/env Rscript
# kegg_entry — KEGGREST keggGet(): fetch a compound entry's OWN record by id.
# This is the back-check that catches a plausible-but-wrong id: the entry's NAME field is
# the DB's own statement of what the id is (C05472 = "Urocortisol", not a bile acid), and
# the linkage/anomer wording in NAME is the isomer evidence a mass check cannot give.
# stdin  JSON: {"ids":["C05465","C04922"]}
# stdout JSON: [{"kegg":..,"names":[..],"formula":..,"exact_mass":..,"mol_weight":..,
#                "chebi":..,"pubchem":..,"lipidmaps":..,"found":true}]
suppressMessages({library(KEGGREST); library(jsonlite)})

req <- fromJSON(paste(readLines("stdin", warn = FALSE), collapse = "\n"),
                simplifyVector = FALSE)
ids <- unique(unlist(req$ids))
ids <- ids[!is.na(ids) & nzchar(ids)]

dblink <- function(entry, key) {
  dl <- entry$DBLINKS
  if (is.null(dl)) return(NULL)
  hit <- grep(paste0("^", key, ":"), dl, value = TRUE)
  if (length(hit) == 0) return(NULL)
  trimws(sub(paste0("^", key, ":"), "", hit[1]))
}

rows <- list()
# keggGet accepts at most 10 ids per request.
for (start in seq(1, max(length(ids), 1), by = 10)) {
  chunk <- ids[start:min(start + 9, length(ids))]
  chunk <- chunk[!is.na(chunk)]
  if (length(chunk) == 0) next
  res <- tryCatch(keggGet(paste0("cpd:", chunk)), error = function(e) NULL)
  got <- character(0)
  if (!is.null(res)) {
    for (entry in res) {
      cid <- sub("^cpd:", "", entry$ENTRY[[1]])
      if (is.null(cid) || is.na(cid)) cid <- entry$ENTRY[[1]]
      got <- c(got, cid)
      nm <- unname(entry$NAME)
      rows[[length(rows) + 1]] <- list(
        kegg = cid, found = TRUE,
        names = as.list(trimws(sub(";$", "", nm))),
        formula = if (is.null(entry$FORMULA)) NULL else entry$FORMULA[[1]],
        exact_mass = if (is.null(entry$EXACT_MASS)) NULL else as.numeric(entry$EXACT_MASS[[1]]),
        mol_weight = if (is.null(entry$MOL_WEIGHT)) NULL else as.numeric(entry$MOL_WEIGHT[[1]]),
        chebi = dblink(entry, "ChEBI"), pubchem = dblink(entry, "PubChem"),
        lipidmaps = dblink(entry, "LIPIDMAPS")
      )
    }
  }
  for (miss in setdiff(chunk, got)) {
    rows[[length(rows) + 1]] <- list(kegg = miss, found = FALSE, names = list(),
                                     formula = NULL, exact_mass = NULL, mol_weight = NULL,
                                     chebi = NULL, pubchem = NULL, lipidmaps = NULL)
  }
}
writeLines(toJSON(rows, auto_unbox = TRUE, null = "null", na = "null"), req$out)
