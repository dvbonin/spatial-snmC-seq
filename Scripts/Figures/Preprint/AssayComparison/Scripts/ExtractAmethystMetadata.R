#!/usr/bin/env Rscript
# Per-cell metadata of a published sci-MET v3 amethyst object, as a flat table.
#
# GEO ships the objects as gzip around either an RDS (NovaSeq) or an .RData holding a
# whole workspace (pilot), so both are tried; the caller gunzips first. The amethyst
# package is not needed: the object's slots are read as plain attributes, and only the
# metadata slot is kept. cg_cov is covered CpG cytosines per cell, per strand.

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--input", type = "character"),
  make_option("--output", type = "character")
)))

obj <- tryCatch(readRDS(opt$input), error = function(e) NULL)
if (is.null(obj)) {
  env <- new.env()
  load(opt$input, envir = env)
  hits <- Filter(function(n) inherits(get(n, env), "amethyst"), ls(env, all.names = TRUE))
  if (length(hits) != 1) stop("expected one amethyst object in ", opt$input, ", found ", length(hits))
  obj <- get(hits, env)
}
md <- attributes(obj)$metadata
if (is.null(md)) stop("no metadata slot in ", opt$input)

out <- data.frame(cell = rownames(md), md, row.names = NULL, check.names = FALSE)
write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %s (%d cells; median cg_cov %.0f)\n", opt$output, nrow(out), median(out$cg_cov)))
