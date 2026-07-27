#!/usr/bin/env Rscript
# bridge_xref — BridgeDb metabolite xref bridging (batch, one DB load per call).
# stdin  JSON: {"db":"<path>", "queries":[{"id":"..","source":"InChIKey","targets":["KEGG Compound","HMDB"]}]}
# stdout JSON: [{"id":..,"source":..,"mappings":{"KEGG Compound":[..],"HMDB":[..]}}]
suppressMessages({library(BridgeDbR); library(jsonlite)})

req <- fromJSON(paste(readLines("stdin", warn = FALSE), collapse = "\n"),
                simplifyVector = FALSE)
db <- req$db
mapper <- loadDatabase(db)
message("loadDatabase OK: ", db)

# map one id (source system name) to a target system name -> character vector
map1 <- function(src_name, id, tgt_name) {
  sc <- tryCatch(getSystemCode(src_name), error = function(e) "")
  tc <- tryCatch(getSystemCode(tgt_name), error = function(e) "")
  if (sc == "" || tc == "" || is.null(id) || is.na(id) || id == "") return(character(0))
  r <- tryCatch(map(mapper, source = sc, identifier = as.character(id), target = tc),
                error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(character(0))
  unique(r$mapping)
}

out <- lapply(req$queries, function(q) {
  targets <- unlist(q$targets)
  mp <- list()
  for (t in targets) {
    # ChEBI may be stored bare or prefixed; try both
    ids <- as.character(q$id)
    if (identical(q$source, "ChEBI")) ids <- unique(c(ids, sub("^CHEBI:", "", ids)))
    hits <- unique(unlist(lapply(ids, function(x) map1(q$source, x, t))))
    mp[[t]] <- as.list(hits)
  }
  list(id = q$id, source = q$source, mappings = mp)
})

writeLines(toJSON(out, auto_unbox = TRUE, null = "null"), req$out)
