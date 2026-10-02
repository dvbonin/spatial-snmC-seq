library(optparse)
library(data.table)

# Pearson r between two per-CpG profiles, over the sites both cover.
#
# One pair per invocation so the matrix can be filled by a job array -- each pair
# holds two multi-million-row tables, which is why they are not all done at once.
# Only chrom/start/beta are read; coverage is deliberately unused, as nothing here
# filters on depth.

option_list <- list(
  make_option("--label-a", type = "character"), make_option("--file-a", type = "character"),
  make_option("--label-b", type = "character"), make_option("--file-b", type = "character"),
  make_option("--keys", type = "character", default = NULL,
              help = "Optional chrom,start list restricting which sites count"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

rd <- function(p) fread(p, header = FALSE, select = c(1, 2, 4),
                        col.names = c("chrom", "start", "beta"))

same <- identical(opt[["file-a"]], opt[["file-b"]])
a <- rd(opt[["file-a"]])
d <- if (same) a[, .(chrom, start, beta.x = beta, beta.y = beta)] else {
  b <- rd(opt[["file-b"]]); m <- merge(a, b, by = c("chrom", "start")); rm(b); m
}
rm(a); invisible(gc())

if (!is.null(opt$keys)) {
  k <- fread(opt$keys, header = FALSE, col.names = c("chrom", "start"))
  d <- merge(d, k, by = c("chrom", "start"))
  rm(k); invisible(gc())
}

r <- if (nrow(d) > 1) cor(d$beta.x, d$beta.y) else NA_real_
fwrite(data.table(source1 = opt[["label-a"]], source2 = opt[["label-b"]],
                  r = r, n = nrow(d)), opt$output)
message(opt[["label-a"]], " x ", opt[["label-b"]], ": r=", round(r, 4), " n=", nrow(d))
