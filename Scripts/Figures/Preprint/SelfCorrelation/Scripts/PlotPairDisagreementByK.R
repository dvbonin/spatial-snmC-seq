library(optparse)
library(data.table)
library(ggplot2)

# Realised pair-disagreement rate per coverage, with the all-pairings expectation as a
# thin reference line. Points are single random pairings, so they scatter around it.

option_list <- list(
  make_option("--tables", type = "character", help = "label=csv,label=csv"),
  make_option("--min-sites", type = "integer", default = 100L),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

spec <- strsplit(strsplit(opt$tables, ",")[[1]], "=")
d <- rbindlist(lapply(spec, function(s) fread(s[2])[, plate := s[1]]))
d <- d[n_sites >= opt[["min-sites"]]]

p <- ggplot(d, aes(k, rate, colour = cell_type, linetype = plate,
                   group = paste(cell_type, plate))) +
  geom_line(aes(y = D_expected), linewidth = 0.3, alpha = 0.45) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_colour_manual(values = c(Keratinocytes = "#1B6CA8", MCC = "#D1495B")) +
  labs(x = "Cells per site k (even only)", y = "Share of disagreeing cell pairs",
       colour = "Cell type", linetype = "Plate",
       title = "One random pairing per site, pooled over sites within each k",
       subtitle = "Thin lines: mean over all possible pairings") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

ggsave(opt$output, p, width = 8.5, height = 5, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 8.5, height = 5, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
