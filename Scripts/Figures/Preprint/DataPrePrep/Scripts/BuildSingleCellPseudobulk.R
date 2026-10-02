library(optparse)
library(data.table)

# Pools the QC-passed single-cell wells of one cell type into a per-CpG pseudobulk.
#
# Two poolings, chosen with --weighting, because the figures want different things:
#
#   reads  methylated and total READS are summed across wells and divided once, so a
#          deeply covered cell contributes proportionally more. Column 5 is the read
#          count. This is what the correlation figures compare against bulk profiles.
#
#   wells  each well contributes its own beta once, whatever its depth, and column 5 is
#          the NUMBER OF WELLS covering the CpG. At ~1 read per well per CpG the two are
#          nearly the same thing; this variant states the equal-per-cell intent outright
#          rather than relying on that.
#
# well_class is the selector rather than CellType because it already encodes both:
# it holds the cell type only for CellCount == 1 wells ("Empty" for 0, NA for
# multi-cell), so this pools single cells without a second condition.
#
# Wells are aggregated in batches: the full pool is ~300 million CpG rows across
# both cell types, which will not fit in memory at once. Each batch is collapsed to
# per-CpG sums and folded into the accumulator, keeping peak memory near the size of
# the covered CpG set rather than the size of the input.

option_list <- list(
  make_option("--well-stats", type = "character",
              help = "Comma-separated beds_dir=well_stats_csv pairs, one per plate"),
  make_option("--cell-type", type = "character",
              help = "well_class value to pool (e.g. HealthyKeratinocytes)"),
  make_option("--batch-size", type = "integer", default = 25L,
              help = "Wells aggregated per batch [default: %default]"),
  make_option("--weighting", type = "character", default = "reads",
              help = "'reads': read-weighted mean, column 5 is reads; 'wells': equal weight per well, column 5 is wells [default: %default]"),
  make_option("--output", type = "character",
              help = "Per-CpG bed.gz: chrom, start, end, beta, coverage")
)
opt <- parse_args(OptionParser(option_list = option_list))

wells <- rbindlist(lapply(strsplit(opt[["well-stats"]], ",")[[1]], function(spec) {
  parts <- strsplit(spec, "=")[[1]]
  ws <- fread(parts[2])
  sel <- ws[well_class == opt[["cell-type"]] & passes_cpg_qc == TRUE, WellPosition]
  message(basename(parts[2]), ": ", length(sel), " QC-passed ", opt[["cell-type"]], " wells")
  data.table(path = file.path(parts[1], paste0(sel, "_MethCPG.bed.gz")))
}))
stopifnot(nrow(wells) > 0, all(file.exists(wells$path)))

acc <- NULL
collapse <- function(d) d[, .(methylated = sum(methylated), coverage = sum(coverage)),
                          by = .(chrom, start, end)]

batches <- split(wells$path, ceiling(seq_len(nrow(wells)) / opt[["batch-size"]]))
for (i in seq_along(batches)) {
  # Merged-CpG bed columns (vcf2bed | mergecg): chrom, start, end, beta, coverage, detail.
  # Rows are CpG dyads spanning two bases; there is no context column to filter on and
  # none is needed, since vcf2bed emits CG context only.
  # Under 'reads' the methylated count is not stored, so it is reconstructed as
  # beta * coverage and summed unrounded: rounding each well first would inject a
  # per-well error for no benefit. Under 'wells' the well's beta is the contribution and
  # the weight is 1, so the accumulator sums betas and counts wells.
  b <- rbindlist(lapply(batches[[i]], function(p) {
    d <- fread(p, sep = "\t", header = FALSE, showProgress = FALSE)
    if (opt$weighting == "wells")
      d[, .(chrom = V1, start = V2, end = V3, methylated = V4, coverage = 1L)]
    else
      d[, .(chrom = V1, start = V2, end = V3, methylated = V4 * V5, coverage = V5)]
  }))
  acc <- if (is.null(acc)) collapse(b) else collapse(rbind(acc, collapse(b)))
  rm(b); invisible(gc())
  message("  batch ", i, "/", length(batches), ": ", nrow(acc), " CpGs pooled")
}

acc[, beta := methylated / coverage]
setorder(acc, chrom, start)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(acc[, .(chrom, start, end, beta, coverage)], opt$output,
       sep = "\t", col.names = FALSE, compress = "gzip")
message("Wrote: ", opt$output, " (", nrow(acc), " CpGs, ", sum(acc$coverage), " total ",
        if (opt$weighting == "wells") "well observations)" else "reads)")
