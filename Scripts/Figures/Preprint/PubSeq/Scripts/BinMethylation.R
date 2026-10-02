library(optparse)
library(data.table)

# Collapses a per-CpG bed into fixed-width bins.
#
# mean_beta is the plain (unweighted) mean of the per-CpG betas in the bin, not a
# coverage-weighted one. That is the only choice that compares like with like here:
# the public bulk file carries no read depths (its coverage column is a constant 1),
# so a weighted mean would silently reduce to the unweighted one on that side while
# staying weighted on the single-cell side.
#
# total_coverage therefore means different things per input -- summed reads for the
# single-cell pseudobulks, but simply the number of CpG rows for the public file.
# Whatever threshold is applied to it downstream is read depth on one axis and CpG
# count on the other.

option_list <- list(
  make_option("--input", type = "character", help = "Per-CpG bed.gz: chrom,start,end,beta,coverage"),
  make_option("--bin-size", type = "integer", default = 10000L),
  make_option("--primary-only", action = "store_true", default = TRUE,
              help = "Keep chr1-22,X,Y; drops the alt/random contigs the liftOver of the public file left behind, which the single-cell side has no coverage on anyway"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

d <- fread(opt$input, header = FALSE,
           col.names = c("chrom", "start", "end", "beta", "coverage"))
n0 <- nrow(d)
d <- d[!is.na(beta) & !is.na(coverage)]
if (isTRUE(opt[["primary-only"]])) d <- d[chrom %in% paste0("chr", c(1:22, "X", "Y"))]
message("Read ", n0, " CpGs, kept ", nrow(d))

d[, bin_start := floor(start / opt[["bin-size"]]) * opt[["bin-size"]]]
agg <- d[, .(mean_beta = mean(beta), total_coverage = sum(coverage), n_cpgs = .N),
         by = .(chrom, bin_start)]
agg[, bin_end := bin_start + opt[["bin-size"]]]
setcolorder(agg, c("chrom", "bin_start", "bin_end", "mean_beta", "total_coverage", "n_cpgs"))
setorder(agg, chrom, bin_start)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(agg, opt$output, sep = "\t", col.names = FALSE, compress = "gzip")
message("Wrote: ", opt$output, " (", nrow(agg), " bins of ", opt[["bin-size"]], " bp)")
