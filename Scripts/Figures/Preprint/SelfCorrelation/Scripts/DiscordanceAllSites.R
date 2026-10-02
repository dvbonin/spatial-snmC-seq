library(optparse)
library(data.table)

# Pairwise discordance over every usable site, at each site's own coverage.
#
# D = 2MU/(k(k-1)) is the share of cell pairs at a site that disagree, and its
# expectation is 2p(1-p) whatever k is. Nothing therefore has to be subsampled and no
# coverage floor has to be imposed: each site is scored at its native k, in each cell
# type at that type's own k, and the site means are directly comparable.
#
# Sites need k >= 2 in both types -- k = 1 has no pair to compare, and requiring both
# keeps the two types on the same sites so the comparison is not a coverage contrast.
#
# The per-stratum table is a breakdown, not a correction: strata differ in which sites
# they contain, so a trend across them is site composition, not a coverage artefact.

option_list <- list(
  make_option("--kz", type = "character"), make_option("--mcc", type = "character"),
  make_option("--label", type = "character", default = ""),
  make_option("--output", type = "character", help = "CSV, one row per coverage stratum")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
d <- merge(rd(opt$kz)[k >= 2], rd(opt$mcc)[k >= 2],
           by = c("chrom", "start"), suffixes = c("_kz", "_mcc"))

d[, D_kz  := 2 * M_kz  * (k_kz  - M_kz)  / (k_kz  * (k_kz  - 1))]
d[, D_mcc := 2 * M_mcc * (k_mcc - M_mcc) / (k_mcc * (k_mcc - 1))]

tab <- d[, .(n_sites = .N, mean_k_KZ = mean(k_kz), mean_k_MCC = mean(k_mcc),
             D_KZ = mean(D_kz), D_MCC = mean(D_mcc)), by = .(k_min = pmin(k_kz, k_mcc))
         ][order(k_min)]
tab[, ratio := D_MCC / D_KZ]
fwrite(tab, opt$output)

o <- d[, .(D_KZ = mean(D_kz), D_MCC = mean(D_mcc),
           k_KZ = mean(k_kz), k_MCC = mean(k_mcc))]
cat(sprintf("\n%s   -   %s sites with k >= 2 in both types, scored at native coverage\n",
            opt$label, format(nrow(d), big.mark = ",")))
cat(sprintf("  mean coverage:  KZ %.2f cells   MCC %.2f cells\n", o$k_KZ, o$k_MCC))
cat(sprintf("  mean D:         KZ %.5f        MCC %.5f      MCC/KZ %.4f\n\n",
            o$D_KZ, o$D_MCC, o$D_MCC / o$D_KZ))
cat(sprintf("%5s %13s | %7s %7s | %8s %8s %7s\n",
            "k_min", "n_sites", "k KZ", "k MCC", "D KZ", "D MCC", "ratio"))
for (i in which(tab$n_sites >= 100)) cat(sprintf("%5d %13s | %7.2f %7.2f | %8.4f %8.4f %7.4f\n",
  tab$k_min[i], format(tab$n_sites[i], big.mark = ","), tab$mean_k_KZ[i], tab$mean_k_MCC[i],
  tab$D_KZ[i], tab$D_MCC[i], tab$ratio[i]))
cat(sprintf("\nWrote: %s\n", opt$output))
