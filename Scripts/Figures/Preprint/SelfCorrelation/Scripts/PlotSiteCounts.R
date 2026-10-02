library(optparse)
library(data.table)
library(ggplot2)

# How many sites back each coverage level. The count is a property of the matched site
# set, so it is the same for both cell types; one line per input table.

option_list <- list(
  make_option("--tables", type = "character", help = "label=csv,label=csv"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

spec <- strsplit(strsplit(opt$tables, ",")[[1]], "=")
d <- rbindlist(lapply(spec, function(s)
  unique(fread(s[2])[, .(k, n_sites)])[, set := factor(s[1], sapply(spec, `[`, 1))]))

p <- ggplot(d, aes(k, n_sites, colour = set, group = set)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_y_log10(breaks = 10^(0:7), labels = c("1", "10", "100", "1k", "10k", "100k", "1M", "10M"),
                minor_breaks = rep(1:9, 8) * 10^rep(0:7, each = 9)) +
  scale_colour_manual(values = c("#1B6CA8", "#7FB2D9", "#D1495B", "#E9A0A8")) +
  labs(x = "Cells per site k (even only)", y = "Sites", colour = NULL,
       title = "Sites per coverage level") +
  theme_bw(base_size = 12)

ggsave(opt$output, p, width = 8, height = 4.8, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 8, height = 4.8, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
print(dcast(d, k ~ set, value.var = "n_sites"))
