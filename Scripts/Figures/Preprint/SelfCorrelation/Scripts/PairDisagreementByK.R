library(optparse)
library(data.table)

# One realised random pairing of the cells at each site, pooled per coverage.
#
# Sites must carry >= 2 cells of both types, and both types are then cut to the shared
# count k = min(k_KZ, k_MCC) by a hypergeometric draw -- equivalent to picking k of the
# site's cells at random. A coverage level k therefore holds the same sites and the same
# number of cells in both types, so the two are directly comparable.
#
# Only even k survive, since an odd count cannot be split into pairs. The k cells are
# then split once, at random, into k/2 disjoint pairs; each pair agrees or disagrees, so
# a level k holds n_sites * k/2 binary outcomes and the reported rate is their mean.
#
# The pairing is not simulated cell by cell. The number of discordant pairs d in a random
# perfect matching depends only on (k, m): with u = k - m unmethylated cells,
#
#   P(d) = C(m,d) C(u,d) d! (m-d-1)!! (u-d-1)!! / (k-1)!!
#
# so d is drawn from that distribution per site, which is the same random variable as
# shuffling the cells. Its mean is 2mu/(k(k-1)) = D, reported alongside as the
# all-pairings expectation the single draw scatters around.

option_list <- list(
  make_option("--kz", type = "character"), make_option("--mcc", type = "character"),
  make_option("--blacklist", type = "character", default = NULL,
              help = "BED of regions to drop before anything else (chrom, start, end)"),
  make_option("--label", type = "character", default = ""),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
d <- merge(rd(opt$kz)[k >= 2], rd(opt$mcc)[k >= 2],
           by = c("chrom", "start"), suffixes = c("_kz", "_mcc"))
if (!is.null(opt$blacklist)) {
  bl <- fread(opt$blacklist, header = FALSE, col.names = c("chrom", "start", "end"))
  setkey(bl, chrom, start, end)
  d[, end := start + 1L]
  hit <- foverlaps(d, bl, by.x = c("chrom", "start", "end"), which = TRUE, nomatch = NULL)
  n0 <- nrow(d)
  if (length(hit$xid)) d <- d[-unique(hit$xid)]
  d[, end := NULL]
  message(sprintf("blacklist: dropped %s of %s sites (%.2f%%)",
                  format(n0 - nrow(d), big.mark = ","), format(n0, big.mark = ","),
                  100 * (n0 - nrow(d)) / n0))
}

d[, k := pmin(k_kz, k_mcc)]
d <- d[k %% 2 == 0]
d[, m_kz  := rhyper(.N, M_kz,  k_kz  - M_kz,  k)]
d[, m_mcc := rhyper(.N, M_mcc, k_mcc - M_mcc, k)]

# log of the odd double factorials, ldf(n) = log(n!!) for odd n, ldf(-1) = 0
tab <- cumsum(log(seq(1, max(d$k) + 1, by = 2)))
ldf <- function(n) ifelse(n < 1, 0, tab[(n + 1) / 2])

# one draw of the discordant-pair count for every site in a (k, m) group
draw_d <- function(k, m, n) {
  u <- k - m
  if (m == 0 || u == 0) return(rep(0L, n))
  dv <- seq(m %% 2, min(m, u), by = 2)
  lp <- lchoose(m, dv) + lchoose(u, dv) + lgamma(dv + 1) +
        ldf(m - dv - 1) + ldf(u - dv - 1) - ldf(k - 1)
  sample(dv, n, replace = TRUE, prob = exp(lp - max(lp)))
}

one_type <- function(mcol, name) {
  g <- d[, .N, by = c("k", mcol)]
  setnames(g, mcol, "m")
  g[, dis := mapply(function(k, m, n) sum(draw_d(k, m, n)), k, m, N)]
  g[, D_sum := N * 2 * m * (k - m) / (k * (k - 1))]
  out <- g[, .(cell_type = name, n_sites = sum(N), n_pairs = sum(N) * k[1] / 2,
               n_discordant = sum(dis), D_expected = sum(D_sum) / sum(N)), by = k][order(k)]
  out[, rate := n_discordant / n_pairs]
  out[]
}

res <- rbind(one_type("m_kz", "Keratinocytes"), one_type("m_mcc", "MCC"))
fwrite(res, opt$output)

cat(sprintf("\n%s   -   %s sites, matched to k = min(k_KZ, k_MCC), even k only\n",
            opt$label, format(nrow(d), big.mark = ",")))
cat(sprintf("%-14s %4s %12s %14s %14s %9s %9s\n",
            "cell type", "k", "sites", "pairs", "discordant", "rate", "E[D]"))
for (i in seq_len(nrow(res))) if (res$n_sites[i] >= 100)
  cat(sprintf("%-14s %4d %12s %14s %14s %9.5f %9.5f\n", res$cell_type[i], res$k[i],
      format(res$n_sites[i], big.mark = ","), format(res$n_pairs[i], big.mark = ","),
      format(res$n_discordant[i], big.mark = ","), res$rate[i], res$D_expected[i]))
o <- res[, .(pairs = sum(n_pairs), dis = sum(n_discordant)), by = cell_type]
cat("\n")
for (i in seq_len(nrow(o))) cat(sprintf("  %-14s pooled over all k: %s / %s pairs = %.5f\n",
  o$cell_type[i], format(o$dis[i], big.mark = ","), format(o$pairs[i], big.mark = ","),
  o$dis[i] / o$pairs[i]))
cat(sprintf("\nWrote: %s\n", opt$output))
