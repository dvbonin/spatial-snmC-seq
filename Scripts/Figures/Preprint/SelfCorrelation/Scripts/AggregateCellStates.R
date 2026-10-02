library(optparse)
library(data.table)

# Per CpG site, how many QC-passed single cells of one cell type cover it (k) and how
# many of those call it methylated (M). Those two numbers are all the pairwise
# discordance downstream needs, which is why no site x cell matrix is built -- that
# would be ~30M x 130 entries for no gain.
#
# Each cell contributes one binary call per site. Nearly always it has a single read
# there, so the call is that read; where it has several and they disagree the majority
# wins, with ties broken at random (reported, and rare). Reads are not pooled across
# cells here: the question is whether cells agree with each other, so each cell gets
# one vote regardless of its depth.
#
# Cells are aggregated in batches -- the union across cells runs to tens of millions
# of sites, and the per-cell tables cannot all be held at once.

option_list <- list(
  make_option("--well-stats", type = "character",
              help = "Comma-separated beds_dir=well_stats_csv pairs, one per plate"),
  make_option("--cell-type", type = "character",
              help = "well_class value to pool (holds the type only for CellCount == 1 wells)"),
  make_option("--mixed-cells", type = "character", default = "majority",
              help = "A cell whose reads disagree at a site (beta not 0 or 1): 'majority' takes the majority read with a coin flip on ties, 'exclude' drops that cell from that site and counts it [default: %default]"),
  make_option("--batch-size", type = "integer", default = 25L),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output", type = "character", help = "chrom, start, k, M")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

wells <- rbindlist(lapply(strsplit(opt[["well-stats"]], ",")[[1]], function(spec) {
  parts <- strsplit(spec, "=")[[1]]
  ws <- fread(parts[2])
  sel <- ws[well_class == opt[["cell-type"]] & passes_cpg_qc == TRUE, WellPosition]
  message(basename(parts[2]), ": ", length(sel), " QC-passed ", opt[["cell-type"]], " wells")
  data.table(path = file.path(parts[1], paste0(sel, "_MethCPG.bed.gz")))
}))
stopifnot(nrow(wells) > 0, all(file.exists(wells$path)))

n_mixed <- 0L
collapse <- function(d) d[, .(k = sum(k), M = sum(M)), by = .(chrom, start)]

acc <- NULL
batches <- split(wells$path, ceiling(seq_len(nrow(wells)) / opt[["batch-size"]]))
for (i in seq_along(batches)) {
  b <- rbindlist(lapply(batches[[i]], function(p) {
    # Merged-CpG bed: beta is column 4 (rows are CpG dyads, not single cytosines).
    d <- fread(p, header = FALSE, select = c(1, 2, 4), showProgress = FALSE,
               col.names = c("chrom", "start", "beta"))
    if (opt[["mixed-cells"]] == "exclude") {
      # A cell only votes where its reads agree outright; anything in between is a
      # cell that disagrees with itself and carries no direction.
      mixed <- d$beta > 0 & d$beta < 1
      n_mixed <<- n_mixed + sum(mixed)
      d <- d[!mixed]
      d[, M := as.integer(beta == 1)]
    } else {
      ties <- d$beta == 0.5
      n_mixed <<- n_mixed + sum(ties)
      d[, M := as.integer(beta > 0.5)]          # majority read
      if (any(ties)) d[ties, M := as.integer(runif(.N) < 0.5)]   # coin flip on an even split
    }
    d[, .(chrom, start, k = 1L, M)]
  }))
  acc <- if (is.null(acc)) collapse(b) else collapse(rbind(acc, collapse(b)))
  rm(b); invisible(gc())
  message("  batch ", i, "/", length(batches), ": ", nrow(acc), " sites")
}

setorder(acc, chrom, start)
fwrite(acc, opt$output, sep = "\t", col.names = FALSE, compress = "gzip")
message("Wrote: ", opt$output, " (", nrow(acc), " sites)")
message("  mixed cell-site observations (", opt[["mixed-cells"]], "): ", n_mixed)
message("  sites with k >= 2: ", sum(acc$k >= 2))
