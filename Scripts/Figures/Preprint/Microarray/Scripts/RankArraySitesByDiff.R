library(optparse)
library(data.table)

# Ranks every array probe by how differently the two reference cell types are
# methylated there, and writes out the top percentile as a BED.
#
# diff is beta_kc - beta_mc and the ranking is on |diff|: sites separating the two
# references in either direction are equally useful for correlating a single cell
# against them. Sites where both references agree carry no cell-type signal at all,
# which is why the correlation is computed on this subset rather than genome-wide.
#
# The full ranked table is kept so any other cutoff can be taken later by filtering
# on `percentile`, without redoing the ranking.

option_list <- list(
  make_option("--input", type = "character",
              help = "hg38 BED from liftOver: chrom,start,end,beta_mc,beta_kc"),
  make_option("--output-ranked", type = "character"),
  make_option("--output-top-bed", type = "character"),
  make_option("--percentile", type = "integer", default = 1L),
  make_option("--primary-only", action = "store_true", default = TRUE,
              help = "Drop alt/random/chrUn contigs -- the pipeline routes those reads into a separate stream, so a single cell has no coverage there to correlate against")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["output-ranked"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt[["output-top-bed"]]), recursive = TRUE, showWarnings = FALSE)

d <- fread(opt$input, header = FALSE,
           col.names = c("chrom", "start", "end", "beta_mc", "beta_kc"))
message("Read ", nrow(d), " lifted probes")

d <- d[!is.na(beta_mc) & !is.na(beta_kc)]
if (isTRUE(opt[["primary-only"]])) {
  keep <- paste0("chr", c(1:22, "X", "Y"))
  n0 <- nrow(d); d <- d[chrom %in% keep]
  message("  dropped ", n0 - nrow(d), " probes on non-primary contigs")
}
# liftOver can map two hg19 probes onto one hg38 base; keep the first so `key`
# stays unique when it is joined against the per-well CpG beds.
n0 <- nrow(d); d <- unique(d, by = c("chrom", "start"))
if (n0 > nrow(d)) message("  dropped ", n0 - nrow(d), " duplicate hg38 positions")

d[, diff := beta_kc - beta_mc]
d[, key := paste0(chrom, ":", start)]
d <- d[order(-abs(diff))]
n <- nrow(d)
d[, rank := .I]
d[, percentile := pmin(100L, ceiling(100 * rank / n))]

fwrite(d[, .(chrom, start, end, beta_mc, beta_kc, diff, key, rank, percentile)],
       opt[["output-ranked"]], compress = "gzip")
message("Wrote: ", opt[["output-ranked"]], " (", n, " probes ranked by |beta_kc - beta_mc|)")

top <- d[percentile <= opt$percentile]
fwrite(top[, .(chrom, start, end, beta_mc, beta_kc, diff, key)],
       opt[["output-top-bed"]], sep = "\t", col.names = FALSE, compress = "gzip")
message("Wrote: ", opt[["output-top-bed"]], " (", nrow(top), " sites, |diff| >= ",
        round(min(abs(top$diff)), 4), ")")
message("  balance: KC-hyper ", sum(top$diff > 0), " / MCC-hyper ", sum(top$diff < 0))
