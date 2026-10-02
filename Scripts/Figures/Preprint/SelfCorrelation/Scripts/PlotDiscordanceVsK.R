library(optparse)
library(data.table)
library(ggplot2)

# Two panels on the same fixed site sets and the same draws: pairwise discordance D,
# which is flat in k, next to the minority share min(M,U)/k, which is not.

option_list <- list(
  make_option("--tables", type = "character", help = "label=csv,label=csv"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

spec <- strsplit(strsplit(opt$tables, ",")[[1]], "=")
d <- rbindlist(lapply(spec, function(s) fread(s[2])[, plate := s[1]]))

m <- melt(d, id.vars = c("plate", "k"),
          measure.vars = c("D_KZ", "D_MCC", "min_KZ", "min_MCC"))
m[, cell_type := factor(grepl("MCC", variable), c(FALSE, TRUE), c("Keratinocytes", "MCC"))]
m[, metric := factor(grepl("^D_", variable), c(TRUE, FALSE),
                     c("Pairwise discordance  D = 2MU / (k(k-1))",
                       "Minority share  min(M,U) / k"))]

p <- ggplot(m, aes(k, value, colour = cell_type, linetype = plate,
                   group = paste(cell_type, plate))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_colour_manual(values = c(Keratinocytes = "#1B6CA8", MCC = "#D1495B")) +
  labs(x = "Cells subsampled per site (k)", y = "Mean over sites",
       colour = "Cell type", linetype = "Plate",
       title = "Same sites at every k: D is flat in coverage, min(M,U)/k is not") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

ggsave(opt$output, p, width = 9.5, height = 4.3, dpi = 300)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p, width = 9.5, height = 4.3, dpi = 300, device = grDevices::svg)
cat(sprintf("Wrote: %s\n", opt$output))
