library(optparse)
library(data.table)

# Pools the bulk-plate wells of one cell type into a per-CpG profile.
#
# Wells are combined coverage-weighted, which is the same thing as summing counts:
#   sum(meth) / sum(cov)  ==  sum( beta_w * cov_w / sum(cov) )
# so a well with coverage 2 at a site contributes a fifth as much as one with
# coverage 10. Counts are summed rather than the weights applied directly because
# the bismark .cov file already stores them exactly.
#
# The bulk plate was processed with bismark, not biscuit, so the input here is
# .bismark.cov.gz (chrom, pos, pos, pct, meth, unmeth) with 1-based positions --
# converted to 0-based BED on the way out so it matches the single-cell and public
# files. The plain .txt.gz bedGraph beside it carries only a percentage, no depth,
# and so cannot be pooled this way.
#
# The two cytosines of each CpG dyad are summed into one site, matching the merged
# single-cell beds -- CorrHeatmaps joins the two on coordinates, so they have to
# agree. Bismark .cov carries no strand column and adjacency cannot settle it (in a
# CGCG stretch positions 2 and 3 are both CpG cytosines but sit in different
# dyads), so the dyad is looked up against the reference CpG set from
# BuildCpGSites.py: a position is the plus-strand C if it is in that set, otherwise
# it is the minus-strand C of the dyad one base earlier.
#
# Wells are aggregated in batches: pooling every well of a type at once is hundreds
# of millions of CpG rows, which will not fit in memory.

option_list <- list(
  make_option("--annotation", type = "character",
              help = "Well stats CSV; columns 1-3 are read positionally as well, cell_type, cells"),
  make_option("--cov-dir", type = "character", help = "Directory holding *.bismark.cov.gz"),
  make_option("--cpg-sites", type = "character",
              help = "CpGSites_hg38.bed.gz from BuildCpGSites.py: plus-strand C of every dyad"),
  make_option("--cell-type", type = "character", help = "cell_type value to pool"),
  make_option("--min-cells", type = "integer", default = 5L,
              help = "Wells with at least this many cells are pooled [default: %default]"),
  make_option("--batch-size", type = "integer", default = 8L),
  make_option("--output", type = "character", help = "Per-CpG bed.gz: chrom,start,end,beta,coverage")
)
opt <- parse_args(OptionParser(option_list = option_list))

ann <- fread(opt$annotation)
setnames(ann, 1:3, c("well", "cell_type", "cells"))
sel <- ann[cell_type == opt[["cell-type"]] & cells >= opt[["min-cells"]]]
if (nrow(sel) == 0)
  stop("No wells matched cell_type '", opt[["cell-type"]], "' with >= ", opt[["min-cells"]],
       " cells. Present: ", paste(unique(ann$cell_type), collapse = ", "))
message("Pooling ", nrow(sel), " wells (", opt[["cell-type"]], ", >= ", opt[["min-cells"]],
        " cells), cell counts: ", paste(sort(unique(sel$cells)), collapse = "/"))

# The well token is always followed by "_S" in the filename, so "A1" cannot match "A10".
cov_files <- list.files(opt[["cov-dir"]], pattern = "\\.bismark\\.cov\\.gz$", full.names = TRUE)
paths <- vapply(sel$well, function(w) {
  hit <- cov_files[grepl(paste0("CpG_context_", w, "_S"), basename(cov_files), fixed = TRUE)]
  if (length(hit) != 1) stop("Expected exactly 1 cov file for well ", w, ", found ", length(hit))
  hit
}, character(1))

cpg <- fread(opt[["cpg-sites"]], header = FALSE, select = c(1, 2),
             col.names = c("chrom", "start"), showProgress = FALSE)
setkey(cpg, chrom, start)
message("Reference CpG dyads: ", nrow(cpg))

acc <- NULL
unassigned <- 0L
collapse <- function(d) d[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chrom, pos)]
batches <- split(paths, ceiling(seq_along(paths) / opt[["batch-size"]]))
for (i in seq_along(batches)) {
  b <- rbindlist(lapply(batches[[i]], function(p) {
    d <- fread(p, header = FALSE, showProgress = FALSE,
               col.names = c("chrom", "start", "end", "pct", "meth", "unmeth"))
    # 1-based cytosine -> 0-based; then fold onto the dyad's plus-strand C.
    d[, p0 := start - 1L]
    d[, self := !is.na(cpg[.(d$chrom, d$p0), which = TRUE])]
    d[, prev := !is.na(cpg[.(d$chrom, d$p0 - 1L), which = TRUE])]
    d[, pos := fifelse(self, p0, p0 - 1L)]
    # A handful of rows per well sit at positions that are not a CpG in the reference
    # at all -- neither (p0,p0+1) nor (p0-1,p0) reads CG. They are ~0.4% of rows and
    # almost all coverage 1. There is no dyad to assign them to, so they are dropped
    # rather than given a coordinate that means nothing.
    unassigned <<- unassigned + d[!(self | prev), .N]
    d[(self | prev), .(chrom, pos, meth, unmeth)]
  }))
  acc <- if (is.null(acc)) collapse(b) else collapse(rbind(acc, collapse(b)))
  rm(b); invisible(gc())
  message("  batch ", i, "/", length(batches), ": ", nrow(acc), " CpGs pooled")
}

acc[, `:=`(coverage = meth + unmeth, beta = meth / (meth + unmeth))]
acc[, `:=`(start = pos, end = pos + 2L)]
setorder(acc, chrom, start)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(acc[, .(chrom, start, end, beta, coverage)], opt$output,
       sep = "\t", col.names = FALSE, compress = "gzip")
message("Dropped ", unassigned, " rows at positions that are not a reference CpG")
message("Wrote: ", opt$output, " (", nrow(acc), " CpG dyads, ", sum(acc$coverage), " total reads)")
