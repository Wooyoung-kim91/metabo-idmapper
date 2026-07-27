#!/usr/bin/env Rscript
# metaboanalyst_map — MetaboAnalystR batch name -> HMDB/KEGG/PubChem/ChEBI.
# Verified invocation (skill 02): InitDataObjects(..., default.dpi = 72) is REQUIRED
# (works around a reproducible recursion bug — do not remove).
# stdin  JSON: {"names":["Taurine","Glucose"], "workdir":"/abs/scratch"}
# stdout JSON: [{"query":..,"match":..,"kegg":..,"hmdb":..,"pubchem":..,"chebi":..,"smiles":..,"status":..}]
suppressMessages({library(MetaboAnalystR); library(jsonlite)})

req <- fromJSON(paste(readLines("stdin", warn = FALSE), collapse = "\n"),
                simplifyVector = FALSE)
cmpds <- unlist(req$names)
workdir <- if (is.null(req$workdir)) file.path(tempdir(), "ma_work") else req$workdir
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
# Robustness: if NONE of the queried names match, MetaboAnalystR's CrossReferencing hits an
# all-fail branch and dies ("object 'current.msg' not found"). Inject guaranteed-match anchors
# so it always proceeds; anchor rows are dropped from the output below.
ANCHORS <- c("Taurine", "L-Glutamic acid", "D-Glucose")
uniq <- unique(c(ANCHORS, cmpds))

old <- getwd(); on.exit(setwd(old), add = TRUE)
setwd(workdir)
mSet <- InitDataObjects("conc", "msetora", FALSE, default.dpi = 72)
mSet <- Setup.MapData(mSet, uniq)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)
res <- as.data.frame(mSet$dataSet$map.table, stringsAsFactors = FALSE)
setwd(old)

nz <- function(x) !(is.na(x) | x %in% c("", "NA"))
col <- function(df, n) if (n %in% colnames(df)) df[[n]] else rep(NA, nrow(df))

rows <- lapply(seq_len(nrow(res)), function(i) {
  q <- as.character(col(res, "Query")[i])
  if (q %in% ANCHORS) return(NULL)   # drop injected anchors
  list(
    query   = q,
    match   = as.character(col(res, "Match")[i]),
    kegg    = as.character(col(res, "KEGG")[i]),
    hmdb    = as.character(col(res, "HMDB")[i]),
    pubchem = as.character(col(res, "PubChem")[i]),
    chebi   = as.character(col(res, "ChEBI")[i]),
    smiles  = as.character(col(res, "SMILES")[i]),
    status  = if (nz(col(res, "Match")[i])) "matched" else "no_match"
  )
})
rows <- rows[!vapply(rows, is.null, logical(1))]
writeLines(toJSON(rows, auto_unbox = TRUE, null = "null", na = "null"), req$out)
