library(optparse)
library(data.table)

# Pairwise methylation discordance per CpG, per cell type:
#
#     D = 2MU / (k(k-1))
#
# the fraction of cell pairs at that site that disagree. It is an unbiased estimate of
# 2p(1-p) -- the chance two independently drawn cells disagree -- at every k, because
# the k(k-1) denominator counts distinct pairs. What changes with k is granularity,
# not expectation: at k=2 only {0,1} is attainable, at k=50 a fine grid.
#
# That granularity is why the two cell types are matched site by site. The SHAPE of D's
# distribution at fixed k is fixed by k and p, so unequal k between the types would
# make their densities differ even with identical biology. Each site is therefore cut
# down to k_match = min(k_KZ, k_MCC) in both types.
#
# Subsampling is done by a hypergeometric draw rather than by picking cells: taking
# k_match of a type's k cells and counting the methylated ones IS a draw from
# Hypergeometric(k, M, k_match), so no per-cell table is needed. When a type already
# has k == k_match the draw is deterministic and returns M unchanged.
#
# One draw per site, never averaged over draws: averaging would pull each site's D
# toward the mean and narrow the very distribution being measured. Single draws are
# noisy per site and unbiased across millions of them.
#
# Unmatched D is also reported, since the MEAN needs no matching (unbiased at any k)
# and matching only throws away cells from the better-covered type.

option_list <- list(
  make_option("--kz", type = "character", help = "chrom,start,k,M for keratinocytes"),
  make_option("--mcc", type = "character", help = "chrom,start,k,M for Merkel carcinoma"),
  make_option("--min-cells", type = "integer", default = 2L,
              help = "Sites need at least this many cells in BOTH types [default: %default]"),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--output-sites", type = "character", help = "Per-site table, for re-analysis"),
  make_option("--output-hist", type = "character", help = "(cell_type, k_match, D) counts"),
  make_option("--output-summary", type = "character", help = "Mean D, overall and by k"),
  make_option("--output-cellspersite", type = "character",
              help = "Sites per cell count: unmatched per type, and matched")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
for (o in c("output-sites", "output-hist", "output-summary", "output-cellspersite"))
  dir.create(dirname(opt[[o]]), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, col.names = c("chrom", "start", "k", "M"))
kz <- rd(opt$kz); mcc <- rd(opt$mcc)
message("KZ sites: ", nrow(kz), "   MCC sites: ", nrow(mcc))

d <- merge(kz[, .(chrom, start, k_kz = k, M_kz = M)],
           mcc[, .(chrom, start, k_mcc = k, M_mcc = M)], by = c("chrom", "start"))
rm(kz, mcc); invisible(gc())
message("Covered in both: ", nrow(d))

mn <- opt[["min-cells"]]
d <- d[k_kz >= mn & k_mcc >= mn]
message("With >= ", mn, " cells in both: ", nrow(d))
fwrite(d, opt[["output-sites"]], sep = "\t", compress = "gzip")
message("Wrote: ", opt[["output-sites"]])

disc <- function(M, k) 2 * M * (k - M) / (k * (k - 1))

# unmatched: each type keeps all of its own cells
d[, D_kz_unmatched := disc(M_kz, k_kz)]
d[, D_mcc_unmatched := disc(M_mcc, k_mcc)]

# matched: both types cut to the smaller cell count at that site
d[, k_match := pmin(k_kz, k_mcc)]
d[, M_kz_s := rhyper(.N, M_kz, k_kz - M_kz, k_match)]
d[, M_mcc_s := rhyper(.N, M_mcc, k_mcc - M_mcc, k_match)]
d[, D_kz := disc(M_kz_s, k_match)]
d[, D_mcc := disc(M_mcc_s, k_match)]

# D is discrete given k_match, so the histogram is exact rather than binned
# cell_type is added after grouping: a length-1 constant cannot go in `by`, whose
# entries must all be as long as the table.
hist <- rbind(
  d[, .(n_sites = .N), by = .(k_match, D = D_kz)][, cell_type := "Keratinocytes"],
  d[, .(n_sites = .N), by = .(k_match, D = D_mcc)][, cell_type := "MerkelCarcinoma"])
setcolorder(hist, c("cell_type", "k_match", "D", "n_sites"))
setorder(hist, cell_type, k_match, D)
fwrite(hist, opt[["output-hist"]])
message("Wrote: ", opt[["output-hist"]], " (", nrow(hist), " distinct (k, D) combinations)")

by_k <- rbind(
  d[, .(cell_type = "Keratinocytes", n_sites = .N, mean_D = mean(D_kz)), by = k_match],
  d[, .(cell_type = "MerkelCarcinoma", n_sites = .N, mean_D = mean(D_mcc)), by = k_match])
overall <- data.table(
  cell_type = rep(c("Keratinocytes", "MerkelCarcinoma"), 2),
  matched = rep(c(TRUE, FALSE), each = 2),
  n_sites = nrow(d),
  mean_D = c(mean(d$D_kz), mean(d$D_mcc),
             mean(d$D_kz_unmatched), mean(d$D_mcc_unmatched)),
  mean_k = c(mean(d$k_match), mean(d$k_match), mean(d$k_kz), mean(d$k_mcc)))
# Cells per site. The matched counts are identical for both types by construction
# (k_match is one number per site), so the unmatched per-type counts are the ones that
# actually differ -- MCC cells cover sites more densely.
cps <- rbind(
  d[, .(n_sites = .N), by = .(k = k_kz)][, `:=`(cell_type = "Keratinocytes", matched = FALSE)],
  d[, .(n_sites = .N), by = .(k = k_mcc)][, `:=`(cell_type = "MerkelCarcinoma", matched = FALSE)],
  d[, .(n_sites = .N), by = .(k = k_match)][, `:=`(cell_type = "Matched", matched = TRUE)])
setorder(cps, cell_type, k)
fwrite(cps[, .(cell_type, matched, k, n_sites)], opt[["output-cellspersite"]])
message("Wrote: ", opt[["output-cellspersite"]])

fwrite(rbind(overall[, .(cell_type, matched, k_match = NA_integer_, n_sites, mean_D, mean_k)],
             by_k[, .(cell_type, matched = TRUE, k_match, n_sites, mean_D, mean_k = k_match)]),
       opt[["output-summary"]])
message("Wrote: ", opt[["output-summary"]])
print(overall)
