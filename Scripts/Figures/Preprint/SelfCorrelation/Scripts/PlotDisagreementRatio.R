library(optparse)
library(data.table)
library(ggplot2)

# MCC/KZ ratio of disagreeing cells against matched coverage k, one line per plate.
# Ratio of the counts, not of the percentages -- identical, since both types are
# compared over the same cells_compared within a k.

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
d[, ratio := disagreeing_MCC / disagreeing_KZ]

p <- ggplot(d, aes(k, ratio, colour = plate)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_colour_manual(values = c(P3 = "#1B6CA8", P4 = "#D1495B")) +
  labs(x = "Matched coverage k (cells per site)",
       y = "Disagreeing cells, MCC / KZ",
       colour = "Plate",
       title = sprintf("Relative rate of disagreeing cells (k with >= %s sites)",
                       format(opt[["min-sites"]], big.mark = ","))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(opt$output, p, width = 7, height = 4.5, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 7, height = 4.5, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
print(d[, .(plate, k, n_sites, disagreeing_KZ, disagreeing_MCC, ratio = round(ratio, 4))])
