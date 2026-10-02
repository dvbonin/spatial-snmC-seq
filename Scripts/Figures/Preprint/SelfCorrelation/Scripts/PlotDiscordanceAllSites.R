library(optparse)
library(data.table)
library(ggplot2)

# Pairwise discordance by coverage stratum, at native coverage, with the all-site mean
# as a horizontal reference. Strata contain different sites, so the trend is composition.

option_list <- list(
  make_option("--tables", type = "character", help = "label=csv,label=csv"),
  make_option("--min-sites", type = "integer", default = 100L),
  make_option("--group-label", type = "character", default = "Plate"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

spec <- strsplit(strsplit(opt$tables, ",")[[1]], "=")
d <- rbindlist(lapply(spec, function(s) fread(s[2])[, plate := s[1]]))
# site-count weighted mean over all strata = the all-site mean
ov <- melt(d[, .(D_KZ = sum(D_KZ * n_sites) / sum(n_sites),
                 D_MCC = sum(D_MCC * n_sites) / sum(n_sites)), by = plate],
           id.vars = "plate", variable.name = "variable", value.name = "value")
m <- melt(d[n_sites >= opt[["min-sites"]]], id.vars = c("plate", "k_min"),
          measure.vars = c("D_KZ", "D_MCC"))
for (x in list(m, ov)) x[, cell_type := factor(grepl("MCC", variable), c(FALSE, TRUE),
                                               c("Keratinocytes", "MCC"))]

p <- ggplot(m, aes(k_min, value, colour = cell_type, linetype = plate,
                   group = paste(cell_type, plate))) +
  geom_hline(data = ov, aes(yintercept = value, colour = cell_type, linetype = plate),
             linewidth = 0.4, alpha = 0.6) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(8)) +
  scale_colour_manual(values = c(Keratinocytes = "#1B6CA8", MCC = "#D1495B")) +
  guides(linetype = if (uniqueN(d$plate) > 1) "legend" else "none") +
  labs(x = "Coverage stratum  min(k_KZ, k_MCC)", y = "Pairwise discordance  D",
       colour = "Cell type", linetype = opt[["group-label"]],
       title = "Discordance at native coverage, all sites with >= 2 cells in both types",
       subtitle = "Horizontal lines: all-site mean. Strata hold different sites, so the trend is composition") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

ggsave(opt$output, p, width = 8.5, height = 5, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 8.5, height = 5, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
