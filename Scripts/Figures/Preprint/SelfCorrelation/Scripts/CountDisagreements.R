library(optparse)
library(data.table)

# Raw count of cell-to-cell disagreements, no normalisation.
#
# At a site with M methylated and U unmethylated cells, the minority count min(M,U) is
# the number of cells that disagree with the majority -- so summing it over sites gives
# the total weight of disagreement in each cell type.
#
# Both types are cut to k = min(k_MCC, k_Kera) at every site, so the two totals are
# built from the same number of cells at the same sites and are directly comparable.
#
# min() is not linear, so unlike the mean of D there is no closed form for its
# expectation under downsampling: this needs real draws. One seeded draw is the
# reported figure; --draws repeats it to show how much the total moves between draws.

option_list <- list(
  make_option("--kz", type = "character"), make_option("--mcc", type = "character"),
  make_option("--min-cells", type = "integer", default = 2L),
  make_option("--draws", type = "integer", default = 20L,
              help = "Extra draws used only to report sampling spread [default: %default]"),
  make_option("--seed", type = "integer", default = 42L)
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
kz <- rd(opt$kz); mcc <- rd(opt$mcc)
mc <- opt[["min-cells"]]
cat(sprintf("\nSites with >= 1 unmixed cell:      Kera %s   MCC %s\n",
            format(nrow(kz), big.mark = ","), format(nrow(mcc), big.mark = ",")))
cat(sprintf("Sites with >= %d unmixed cells:     Kera %s   MCC %s\n", mc,
            format(sum(kz$k >= mc), big.mark = ","), format(sum(mcc$k >= mc), big.mark = ",")))

d <- merge(kz[k >= mc], mcc[k >= mc], by = c("chrom", "start"),
           suffixes = c("_kz", "_mcc"))
rm(kz, mcc); invisible(gc())
d[, k := pmin(k_kz, k_mcc)]
cat(sprintf("Sites with >= %d in BOTH types:     %s\n", mc, format(nrow(d), big.mark = ",")))
cat(sprintf("Cells compared per type (sum of k): %s\n\n", format(sum(d$k), big.mark = ",")))

one_draw <- function() {
  m_kz  <- rhyper(nrow(d), d$M_kz,  d$k_kz  - d$M_kz,  d$k)
  m_mcc <- rhyper(nrow(d), d$M_mcc, d$k_mcc - d$M_mcc, d$k)
  c(kz = sum(pmin(m_kz, d$k - m_kz)), mcc = sum(pmin(m_mcc, d$k - m_mcc)))
}

first <- one_draw()
cat("=== TOTAL DISAGREEMENTS (sum of min(M,U) over sites, matched cell counts) ===\n")
cat(sprintf("  Keratinocytes        %15s\n", format(first[["kz"]], big.mark = ",")))
cat(sprintf("  Merkel Cell Carcinoma%15s\n", format(first[["mcc"]], big.mark = ",")))
cat(sprintf("  difference (MCC-Kera)%15s   (%+.2f%% relative to Kera)\n",
            format(first[["mcc"]] - first[["kz"]], big.mark = ","),
            100 * (first[["mcc"]] / first[["kz"]] - 1)))
cat(sprintf("  per compared cell:     Kera %.5f   MCC %.5f\n",
            first[["kz"]] / sum(d$k), first[["mcc"]] / sum(d$k)))

if (opt$draws > 1) {
  reps <- replicate(opt$draws, one_draw())
  cat(sprintf("\n  spread over %d draws:  Kera %s +/- %s   MCC %s +/- %s\n",
              opt$draws,
              format(round(mean(reps["kz", ])), big.mark = ","), format(round(sd(reps["kz", ])), big.mark = ","),
              format(round(mean(reps["mcc", ])), big.mark = ","), format(round(sd(reps["mcc", ])), big.mark = ",")))
  dif <- reps["mcc", ] - reps["kz", ]
  cat(sprintf("  difference over draws: %s +/- %s   (draws with MCC > Kera: %d/%d)\n",
              format(round(mean(dif)), big.mark = ","), format(round(sd(dif)), big.mark = ","),
              sum(dif > 0), opt$draws))
}
