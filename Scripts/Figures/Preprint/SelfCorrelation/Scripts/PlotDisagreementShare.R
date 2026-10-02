library(optparse)
library(data.table)
library(ggplot2)

# Mean per-site share of disagreeing cells against matched coverage k, one line per
# cell type and plate, point size showing how many sites back each k.
#
# Within a fixed k the share at a site is dis_i/k, so the mean over sites is
# sum(dis_i)/(n_sites * k) = disagreeing / cells_compared -- exactly the pct columns
# the coverage table already carries.

option_list <- list(
  make_option("--tables", type = "character", help = "label=csv,label=csv"),
  make_option("--min-sites", type = "integer", default = 100L),
  make_option("--y", type = "character", default = "share",
              help = "'share': disagreeing cells / (sites * k); 'per-site': disagreeing cells / sites [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

spec <- strsplit(strsplit(opt$tables, ",")[[1]], "=")
d <- rbindlist(lapply(spec, function(s) fread(s[2])[, plate := s[1]]))
d <- d[n_sites >= opt[["min-sites"]]]

if (opt$y == "per-site") {
  d[, `:=`(y_KZ = disagreeing_KZ / n_sites, y_MCC = disagreeing_MCC / n_sites)]
  y_lab <- "Mean disagreeing cells per site"
  ttl <- "Disagreeing cells per site"
} else {
  d[, `:=`(y_KZ = pct_KZ, y_MCC = pct_MCC)]
  y_lab <- "Mean share of disagreeing cells per site (%)"
  ttl <- "Per-site disagreement rate"
}

m <- melt(d, id.vars = c("plate", "k", "n_sites"),
          measure.vars = c("y_KZ", "y_MCC"),
          variable.name = "cell_type", value.name = "share")
m[, cell_type := factor(cell_type, c("y_KZ", "y_MCC"), c("Keratinocytes", "MCC"))]
m[, grp := paste(cell_type, plate)]

p <- ggplot(m, aes(k, share, colour = cell_type, linetype = plate, group = grp)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_colour_manual(values = c(Keratinocytes = "#1B6CA8", MCC = "#D1495B")) +
  labs(x = "Matched coverage k (cells per site)",
       y = y_lab,
       colour = "Cell type", linetype = "Plate",
       title = sprintf("%s (k with >= %s sites)", ttl,
                       format(opt[["min-sites"]], big.mark = ","))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(opt$output, p, width = 7.5, height = 4.5, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 7.5, height = 4.5, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
