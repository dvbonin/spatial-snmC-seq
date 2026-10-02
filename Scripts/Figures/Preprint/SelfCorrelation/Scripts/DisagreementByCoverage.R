library(optparse)
library(data.table)

# Disagreeing cells per CpG, stratified by how many cells cover the site.
#
# At a site with M methylated and U unmethylated cells, min(M,U) cells sit in the
# minority -- they disagree with the rest. Summing that within each coverage level k
# keeps the comparison honest: the raw sum over all sites silently up-weights
# well-covered sites, because min(M,U) can reach k/2 while a 2-cell site can only ever
# contribute 1.
#
# Both cell types are cut to k = min(k_KZ, k_MCC) at each site, so a given row compares
# the same number of cells at the same sites in both types. Downsampling is one seeded
# hypergeometric draw per site -- equivalent to picking k cells at random and counting
# the methylated ones.
#
# Only cells whose reads at a site all agree are counted at all; a cell that disagrees
# with itself carries no direction and is excluded upstream.

option_list <- list(
  make_option("--kz", type = "character"), make_option("--mcc", type = "character"),
  make_option("--min-cells", type = "integer", default = 2L),
  make_option("--label", type = "character", default = ""),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output", type = "character", help = "CSV of the table")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
mc <- opt[["min-cells"]]
d <- merge(rd(opt$kz)[k >= mc], rd(opt$mcc)[k >= mc],
           by = c("chrom", "start"), suffixes = c("_kz", "_mcc"))
d[, k := pmin(k_kz, k_mcc)]

# subsample both types to the shared cell count, then count the minority
d[, m_kz  := rhyper(.N, M_kz,  k_kz  - M_kz,  k)]
d[, m_mcc := rhyper(.N, M_mcc, k_mcc - M_mcc, k)]
d[, dis_kz  := pmin(m_kz,  k - m_kz)]
d[, dis_mcc := pmin(m_mcc, k - m_mcc)]

tab <- d[, .(n_sites = .N,
             cells_compared = .N * k[1],
             disagreeing_KZ = sum(dis_kz),
             disagreeing_MCC = sum(dis_mcc)), by = k][order(k)]
tab[, `:=`(pct_KZ = 100 * disagreeing_KZ / cells_compared,
           pct_MCC = 100 * disagreeing_MCC / cells_compared)]
tab[, pct_diff := pct_MCC - pct_KZ]

fwrite(tab, opt$output)

cat("\n================================================================================\n")
cat(sprintf("  %s   -   disagreeing cells by matched coverage\n", opt$label))
cat("================================================================================\n")
cat(sprintf("%4s %13s %15s %13s %13s %9s %9s %8s\n",
            "k", "sites", "cells cmp", "disagree KZ", "disagree MCC", "% KZ", "% MCC", "diff"))
for (i in seq_len(nrow(tab))) cat(sprintf("%4d %13s %15s %13s %13s %8.3f%% %8.3f%% %+7.3f\n",
  tab$k[i], format(tab$n_sites[i], big.mark = ","), format(tab$cells_compared[i], big.mark = ","),
  format(tab$disagreeing_KZ[i], big.mark = ","), format(tab$disagreeing_MCC[i], big.mark = ","),
  tab$pct_KZ[i], tab$pct_MCC[i], tab$pct_diff[i]))
tot <- tab[, .(n_sites = sum(n_sites), cells = sum(cells_compared),
               kz = sum(disagreeing_KZ), mcc = sum(disagreeing_MCC))]
cat(sprintf("%4s %13s %15s %13s %13s %8.3f%% %8.3f%% %+7.3f\n", "all",
            format(tot$n_sites, big.mark = ","), format(tot$cells, big.mark = ","),
            format(tot$kz, big.mark = ","), format(tot$mcc, big.mark = ","),
            100 * tot$kz / tot$cells, 100 * tot$mcc / tot$cells,
            100 * (tot$mcc - tot$kz) / tot$cells))
cat(sprintf("\n  Wrote: %s\n", opt$output))
