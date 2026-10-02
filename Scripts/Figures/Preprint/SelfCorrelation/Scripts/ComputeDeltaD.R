library(optparse)
library(data.table)

# Paired per-CpG difference in cell-to-cell methylation discordance, MCC minus
# keratinocytes:
#
#     D  = 2MU / (k(k-1))          discordant share of cell pairs at a site
#     dD = D_MCC - D_Kera          same site, same number of cells
#
# Pairing at the site is the point: it cancels site-to-site variation in methylation
# level, which dominates the marginal distributions of D and hid the comparison.
# A genome-wide distribution of dD shifted above zero means MCC cells disagree with
# each other more than keratinocytes do at the same CpGs.
#
# Both types are cut to k = min(k_MCC, k_Kera) so the estimator's granularity and
# variance are identical on the two sides of every difference.
#
# Downsampling is done ANALYTICALLY rather than by repeated random draws. Averaging R
# random subsamples converges to a closed form: for a type with N cells of which M are
# methylated, drawing k without replacement,
#
#     E[M'(k-M')] = k*p*(1-p) * [k - (N-k)/(N-1)]        p = M/N
#     E[D]        = 2*p*(1-p) * [k - (N-k)/(N-1)] / (k-1)
#
# so this is the R -> infinity limit: no Monte Carlo noise, one pass, and dD comes out
# continuous instead of confined to the coarse discrete grid a single draw allows.
# --validate-draws checks the algebra against real resampling on a subsample.

option_list <- list(
  make_option("--sites", type = "character",
              help = "PerSiteCellCounts.tsv.gz: chrom,start,k_kz,M_kz,k_mcc,M_mcc"),
  make_option("--min-cells", type = "integer", default = 5L,
              help = "Cells required in BOTH types before a site is used [default: %default]"),
  make_option("--bins", type = "integer", default = 2000L,
              help = "Histogram bins spanning dD in [-1,1] [default: %default]"),
  make_option("--validate-draws", type = "integer", default = 0L,
              help = "If >0, Monte-Carlo check of the closed form with this many draws"),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output-hist", type = "character", help = "(k, dD bin) counts"),
  make_option("--output-summary", type = "character", help = "Overall and per-k statistics")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
for (o in c("output-hist", "output-summary"))
  dir.create(dirname(opt[[o]]), recursive = TRUE, showWarnings = FALSE)

d <- fread(opt$sites)
message("Sites covered in both types (k >= 2): ", nrow(d))
mc <- opt[["min-cells"]]
d <- d[k_kz >= mc & k_mcc >= mc]
message("With >= ", mc, " cells in both: ", nrow(d))
if (!nrow(d)) stop("no sites left at --min-cells ", mc)

d[, k := pmin(k_kz, k_mcc)]

# Expected D after downsampling N cells to k, exactly.
ED <- function(M, N, k) {
  p <- M / N
  2 * p * (1 - p) * (k - (N - k) / (N - 1)) / (k - 1)
}
d[, D_kz := ED(M_kz, k_kz, k)]
d[, D_mcc := ED(M_mcc, k_mcc, k)]
d[, dD := D_mcc - D_kz]

if (opt[["validate-draws"]] > 0) {
  R <- opt[["validate-draws"]]
  v <- d[sample(.N, min(.N, 50000L))]
  disc <- function(M, k) 2 * M * (k - M) / (k * (k - 1))
  mc_kz <- Reduce(`+`, lapply(seq_len(R), function(i)
    disc(rhyper(nrow(v), v$M_kz, v$k_kz - v$M_kz, v$k), v$k))) / R
  mc_mcc <- Reduce(`+`, lapply(seq_len(R), function(i)
    disc(rhyper(nrow(v), v$M_mcc, v$k_mcc - v$M_mcc, v$k), v$k))) / R
  message("Closed form vs ", R, " random draws on ", nrow(v), " sites:")
  message(sprintf("  D_Kera  mean analytic %.6f  MC %.6f   max|diff| %.5f",
                  mean(v$D_kz), mean(mc_kz), max(abs(v$D_kz - mc_kz))))
  message(sprintf("  D_MCC   mean analytic %.6f  MC %.6f   max|diff| %.5f",
                  mean(v$D_mcc), mean(mc_mcc), max(abs(v$D_mcc - mc_mcc))))
  message(sprintf("  dD      mean analytic %.6f  MC %.6f",
                  mean(v$dD), mean(mc_mcc - mc_kz)))
}

# Histogram keyed on k so any coverage restriction can be applied at plot time
w <- 2 / opt$bins
d[, bin := round(floor(dD / w) * w + w / 2, 8)]
hist <- d[, .(n_sites = .N), by = .(k, bin)]
setorder(hist, k, bin)
fwrite(hist, opt[["output-hist"]])
message("Wrote: ", opt[["output-hist"]], " (", nrow(hist), " (k, bin) cells)")

# list(), not .(): the .() alias only works inside data.table's j expression.
# quantile() returns named values, which would become stray column names.
stats <- function(x) list(n_sites = length(x), mean_dD = mean(x), median_dD = median(x),
                          frac_pos = mean(x > 0), frac_neg = mean(x < 0),
                          q05 = unname(quantile(x, 0.05)),
                          q95 = unname(quantile(x, 0.95)))
overall <- d[, stats(dD)][, k := NA_integer_]
by_k <- d[, stats(dD), by = k][order(k)]
fwrite(rbind(overall, by_k, use.names = TRUE), opt[["output-summary"]])
message("Wrote: ", opt[["output-summary"]])
message(sprintf("\n  OVERALL  n=%s  mean dD=%+.5f  median=%+.5f  frac(dD>0)=%.4f",
                format(overall$n_sites, big.mark = ","), overall$mean_dD,
                overall$median_dD, overall$frac_pos))
