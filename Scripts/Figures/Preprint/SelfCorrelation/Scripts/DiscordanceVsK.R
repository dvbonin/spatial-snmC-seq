library(optparse)
library(data.table)

# Disagreement versus coverage on a FIXED set of sites.
#
# The per-coverage table (DisagreementByCoverage.R) compares a different set of sites at
# every k, and scores them with min(M,U)/k, which underestimates the minority at small k.
# Both effects push the curve up with k. Here they are removed:
#
#   fixed site set   only sites with k >= min-k in BOTH cell types, so every k below that
#                    is evaluated on the same sites
#   unbiased metric  pairwise discordance D = 2MU/(k(k-1)), the share of cell pairs at the
#                    site that disagree. E[D] = 2p(1-p) at every k, so it does not drift
#                    with coverage
#
# min(M,U)/k is reported alongside on the same sites and the same draws, as the contrast.
# Each k is an independent subsample of the fixed set, averaged over --reps draws.

option_list <- list(
  make_option("--kz", type = "character"), make_option("--mcc", type = "character"),
  make_option("--min-k", type = "integer", default = 10L,
              help = "cells required in both types; also the largest k plotted [default: %default]"),
  make_option("--reps", type = "integer", default = 10L),
  make_option("--label", type = "character", default = ""),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

K <- opt[["min-k"]]
rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
d <- merge(rd(opt$kz)[k >= K], rd(opt$mcc)[k >= K],
           by = c("chrom", "start"), suffixes = c("_kz", "_mcc"))
message(opt$label, ": ", nrow(d), " sites with k >= ", K, " in both types")

draw <- function(M, kk, k) rhyper(length(M), M, kk - M, k)
tab <- rbindlist(lapply(2:K, function(k) {
  r <- rbindlist(lapply(seq_len(opt$reps), function(i) {
    m_kz  <- draw(d$M_kz,  d$k_kz,  k)
    m_mcc <- draw(d$M_mcc, d$k_mcc, k)
    D   <- function(m) mean(2 * m * (k - m) / (k * (k - 1)))
    mn  <- function(m) mean(pmin(m, k - m) / k)
    data.table(D_KZ = D(m_kz), D_MCC = D(m_mcc), min_KZ = mn(m_kz), min_MCC = mn(m_mcc))
  }))
  cbind(data.table(k = k, n_sites = nrow(d)), r[, lapply(.SD, mean)])
}))

fwrite(tab, opt$output)
cat(sprintf("\n%s   -   %s sites, %d draws per k\n", opt$label,
            format(nrow(d), big.mark = ","), opt$reps))
cat(sprintf("%3s | %8s %8s | %8s %8s\n", "k", "D KZ", "D MCC", "min/k KZ", "min/k MCC"))
for (i in seq_len(nrow(tab))) cat(sprintf("%3d | %8.4f %8.4f | %8.4f %8.4f\n",
  tab$k[i], tab$D_KZ[i], tab$D_MCC[i], tab$min_KZ[i], tab$min_MCC[i]))
cat(sprintf("\nWrote: %s\n", opt$output))
