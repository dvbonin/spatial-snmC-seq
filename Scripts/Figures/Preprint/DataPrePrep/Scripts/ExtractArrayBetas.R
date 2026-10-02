library(minfi)

# EPIC array betas for the two reference cell types, straight from the raw IDATs.
#
# Emits hg19 coordinates -- the EPIC annotation package is hg19-based, so RunAll
# pipes this through liftOver to reach the hg38 the sequencing uses.
#
# preprocessNoob + getBeta is deliberate. getBeta's default offset of 100 damps
# betas at probes with very low total intensity, where M/(M+U) is a ratio of two
# small noisy numbers; without it every beta inflates, and by a different amount
# per chip (+0.067 vs +0.030 on these two), which flips which reference looks
# globally more methylated and skews the |diff| ranking downstream.
#
# force = TRUE is required: the two chips report slightly different probe counts
# (1,051,943 vs 1,051,815), and read.metharray refuses to bind them otherwise.
#
# CAVEAT on which sample is which. The identities are presumed from filename and
# folder, never confirmed against a metadata record: 207003720077 is the in-house
# array taken to be the Merkel carcinoma reference, GSM4586444 the public GEO
# sample taken to be keratinocyte. If that is the wrong way round, every
# downstream corr_to_mcc / corr_to_kz swaps and the UMAP panel's colour direction
# inverts. Assigned explicitly below rather than relying on column order.

# Argument parsing is by hand: this runs in the SNMC_Array env, which carries minfi
# and the EPIC annotations but not optparse.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
opt <- list(
  "idat-dir"     = arg("idat-dir", file.path(Sys.getenv("SNMC_DATA"), "SNMC_B1/Microarrays")),
  # 207003720077 = in-house, presumed Merkel carcinoma; GSM4586444 = public GEO,
  # presumed keratinocyte. See the CAVEAT above.
  "mcc-basename" = arg("mcc-basename", "207003720077_R08C01"),
  "kc-basename"  = arg("kc-basename",  "GSM4586444_201516320034_R01C01"),
  output         = arg("output")
)
stopifnot(!is.null(opt$output))

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

basenames <- file.path(opt[["idat-dir"]], c(opt[["mcc-basename"]], opt[["kc-basename"]]))
message("Reading IDATs:\n  ", paste(basenames, collapse = "\n  "))

rgset <- read.metharray(basenames = basenames, force = TRUE)
rgset@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")

mset  <- preprocessNoob(rgset)
betas <- getBeta(mset)
colnames(betas) <- c("beta_mc", "beta_kc")
message("Betas: ", nrow(betas), " probes x ", ncol(betas), " samples")

anno <- getAnnotation(mset)
stopifnot(identical(rownames(anno), rownames(betas)))

out <- data.frame(
  chrom   = anno$chr,
  start   = anno$pos - 1L,
  end     = anno$pos,
  beta_mc = betas[, "beta_mc"],
  beta_kc = betas[, "beta_kc"]
)
n_all <- nrow(out)
out <- out[!is.na(out$beta_mc) & !is.na(out$beta_kc), ]
out <- out[order(out$chrom, out$start), ]
message("Kept ", nrow(out), " / ", n_all, " probes with a beta in both samples")
message("  mean beta_mc = ", round(mean(out$beta_mc), 4),
        "   mean beta_kc = ", round(mean(out$beta_kc), 4))

write.table(out, opt$output, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
message("Wrote: ", opt$output)
