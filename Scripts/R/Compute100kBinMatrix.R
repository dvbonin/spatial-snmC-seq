library(optparse)
library(tidyverse)
library(data.table)

option_list <- list(
  make_option("--beds-pattern", type = "character",
              help = "Glob pattern for per-well *_MethCPG.bed.gz files"),
  make_option("--prefix",       type = "character",
              help = "Cell ID prefix: <sample>_<plate>"),
  make_option("--output",       type = "character",
              help = "Output gzipped cell x bin matrix CSV"),
  make_option("--bin-size",     type = "integer", default = 100000L,
              help = "Genomic bin size in bp [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

beds_pattern <- opt[["beds-pattern"]]
prefix       <- opt$prefix
output_path  <- opt$output
bin_size     <- opt[["bin-size"]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# BED format (biscuit vcf2bed -e): chr, start, end, base, ctx2, ctx3, seq, beta, coverage
bed_files <- Sys.glob(beds_pattern)
if (length(bed_files) == 0) stop("No BED files matched: ", beds_pattern)
message("Binning ", length(bed_files), " cells into ", bin_size, "bp bins")

agg_list <- vector("list", length(bed_files))
for (i in seq_along(bed_files)) {
  path    <- bed_files[[i]]
  well    <- basename(path) %>% str_remove("_MethCPG\\.bed\\.gz$")
  cell_id <- paste0(prefix, "_", well)

  df <- tryCatch(fread(path, sep = "\t", header = FALSE, data.table = TRUE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  # Merged-CpG bed: beta in column 4, coverage in column 5, rows are dyads. No context
  # column to filter on -- vcf2bed emits CG only.
  if (nrow(df) == 0) next

  df[, bin := paste0(V1, ":", floor(V2 / bin_size) * bin_size)]
  agg <- df[, .(wmean = sum(V4 * V5, na.rm = TRUE) / sum(V5, na.rm = TRUE)), by = bin]
  agg[, cell_id := cell_id]
  agg_list[[i]] <- agg

  if (i %% 20 == 0) message("  ", i, "/", length(bed_files), " cells binned")
}

dt   <- rbindlist(Filter(Negate(is.null), agg_list))
wide <- dcast(dt, cell_id ~ bin, value.var = "wmean")

message("Wrote matrix: ", nrow(wide), " cells x ", ncol(wide) - 1, " bins")
fwrite(wide, output_path, compress = "gzip")
message("Wrote: ", output_path)
