library(optparse)
library(data.table)

# Top differentially methylated CpG sites between the two bulk profiles, at single-site
# resolution.
#
# Bulk is the selector because it is independent of the single-cell data, so choosing
# sites on it does not bias the pseudobulk comparisons the figure is about. Bulk_KZ vs
# Bulk_MCC is consequently self-defining in the resulting plot -- r near -1 by
# construction, not a result.
#
# Ties are the thing to watch here. At a few reads per site a beta is frequently
# exactly 0 or 1, so |diff| = 1 can be shared by far more sites than the cutoff
# admits, in which case the "top 1%" is whichever of them the sort happens to reach.
# The tie situation at the threshold is reported so it is visible rather than implied.

option_list <- list(
  make_option("--bulk-kz", type = "character"),
  make_option("--bulk-mcc", type = "character"),
  make_option("--top-pct", type = "double", default = 1),
  make_option("--min-coverage", type = "integer", default = 3L,
              help = "Reads required at the site in BOTH bulk profiles before ranking [default: %default]"),
  make_option("--output-keys", type = "character", help = "chrom,start of the selected sites")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["output-keys"]]), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, select = c(1, 2, 4, 5),
                        col.names = c("chrom", "start", "beta", "cov"))
a <- rd(opt[["bulk-kz"]]); b <- rd(opt[["bulk-mcc"]])
message("Bulk_KZ sites: ", nrow(a), "   Bulk_MCC sites: ", nrow(b))

d <- merge(a, b, by = c("chrom", "start"))
rm(a, b); invisible(gc())
message("Covered in both: ", nrow(d))

# Depth requirement applies to SELECTION only, and only to the bulk profiles. At one
# or two reads a beta is 0 or 1 outright, so |diff| hits its ceiling by sampling
# accident and the top percentile fills with the least reliable sites in the genome.
# The correlations themselves remain unfiltered, and the public data -- whose depths
# were not preserved -- plays no part here.
mc <- opt[["min-coverage"]]
if (mc > 1) {
  n0 <- nrow(d); d <- d[cov.x >= mc & cov.y >= mc]
  message("With >= ", mc, " reads in both bulk profiles: ", nrow(d),
          " (", round(100 * nrow(d) / n0, 1), "% of shared sites)")
}

d[, abs_diff := abs(beta.x - beta.y)]
setorder(d, -abs_diff)
n_keep <- ceiling(nrow(d) * opt[["top-pct"]] / 100)
cutoff <- d$abs_diff[n_keep]
n_at_cutoff <- sum(d$abs_diff == cutoff)
n_above <- sum(d$abs_diff > cutoff)

message("Top ", opt[["top-pct"]], "% = ", n_keep, " sites; |diff| cutoff = ", round(cutoff, 4))
message("  sites strictly above the cutoff: ", n_above)
message("  sites exactly AT the cutoff:     ", n_at_cutoff,
        "  -> ", n_keep - n_above, " of them taken, ",
        n_at_cutoff - (n_keep - n_above), " dropped by tie-break")
if (n_at_cutoff > 0.5 * n_keep)
  warning("over half the selection sits on a tie; the cut is largely arbitrary",
          call. = FALSE, immediate. = TRUE)

fwrite(d[seq_len(n_keep), .(chrom, start)], opt[["output-keys"]], sep = "\t", col.names = FALSE)
message("Wrote: ", opt[["output-keys"]])
