library(optparse)
library(data.table)
library(ggplot2)
library(patchwork)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# Each 10kb bin as one point: single-cell pseudobulk methylation on x, the public
# bulk keratinocyte dataset on y, one panel per single-cell cell type.
#
# The comparison is the point of the figure: keratinocytes should track the public
# keratinocyte data along the identity line, and Merkel carcinoma should not. The
# gap between the two panels' r is the readout, not either value alone.
#
# n differs slightly between panels because the coverage threshold is applied to
# each pair separately -- a bin kept in one panel can fail in the other.

option_list <- list(
  make_option("--panels", type = "character",
              help = "Comma-separated Label=binned_bed.gz, one per x-axis panel"),
  make_option("--y-label", type = "character"),
  make_option("--y-file", type = "character"),
  make_option("--min-coverage", type = "integer", default = 50L,
              help = "Minimum total_coverage in BOTH samples per bin [default: %default]"),
  make_option("--title", type = "character", default = ""),
  make_option("--tiff-dpi", type = "integer", default = 300L,
              help = "Resolution of the publication TIFF [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
min_cov <- opt[["min-coverage"]]

read_bins <- function(path) {
  d <- fread(path, header = FALSE,
             col.names = c("chrom", "bin_start", "bin_end", "mean_beta", "total_coverage", "n_cpgs"))
  d[, .(key = paste0(chrom, ":", bin_start), mean_beta, total_coverage)]
}

y_dt <- read_bins(opt[["y-file"]])

make_panel <- function(x_label, x_path) {
  merged <- merge(read_bins(x_path), y_dt, by = "key", suffixes = c("_x", "_y"))
  merged <- merged[total_coverage_x >= min_cov & total_coverage_y >= min_cov]
  r_val <- cor(merged$mean_beta_x, merged$mean_beta_y)
  message(x_label, " vs ", opt[["y-label"]], ": r=", round(r_val, 3), ", n=", nrow(merged))

  ggplot(merged, aes(mean_beta_x, mean_beta_y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#e34948") +
    geom_point(alpha = 0.167, size = 0.53, color = "#1c5cab") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = x_label, y = opt[["y-label"]],
         subtitle = paste0("r = ", sprintf("%.2f", r_val),
                           "   n = ", format(nrow(merged), big.mark = ","))) +
    theme_bw(base_size = 11) +
    theme(plot.background = element_rect(fill = "#fcfcfb", color = NA),
          panel.background = element_rect(fill = "#fcfcfb", color = NA))
}

panels <- lapply(strsplit(opt$panels, ",")[[1]], function(spec) {
  parts <- strsplit(spec, "=")[[1]]
  make_panel(parts[1], paste(parts[-1], collapse = "="))
})

combined <- Reduce(`|`, panels) + plot_annotation(title = opt$title)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
save_figure( combined, opt$output, width = 5 * length(panels), height = 5.2, dpi = 200, bg = "#fcfcfb", svg = FALSE)
# Publication raster. This figure draws ~560k points, so the vector formats are
# unusable for it -- the SVG runs to 240 MB and will not open in Illustrator. A
# compressed TIFF at print resolution is the practical deliverable instead.
ggsave(sub("\\.png$", ".tiff", opt$output), combined,
       width = 5 * length(panels), height = 5.2, dpi = opt[["tiff-dpi"]],
       bg = "#fcfcfb", device = "tiff", compression = "lzw")
message("Wrote: ", opt$output)
