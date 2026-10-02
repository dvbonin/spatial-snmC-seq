library(optparse)
library(tidyverse)
library(data.table)
library(patchwork)

# Three views of per-CpG pairwise methylation discordance, keratinocytes vs Merkel
# carcinoma. Reads only the small tabulated intermediates, never the per-site table,
# so restyling costs nothing.
#
#   1 density of D        kernel density over the site-weighted D values, bounded to
#                         [0,1] and renormalised. Note D is discrete given k -- only
#                         {0,1} at k=2 -- so the smooth curve interpolates across
#                         values D cannot actually take; it shows where the mass sits,
#                         not a continuous underlying variable.
#   2 cells per site      the k distribution, which sets how coarse D can be
#   3 mean D by k         one number per k per cell type. Point area is proportional
#                         to how many sites back it: the high-k tail rests on orders
#                         of magnitude fewer sites than the flat low-k region, and
#                         without that encoding the tail reads as a strong trend.

option_list <- list(
  make_option("--hist", type = "character", help = "(cell_type, k_match, D, n_sites)"),
  make_option("--summary", type = "character", help = "Mean D overall and by k"),
  make_option("--cells-per-site", type = "character",
              help = "Sites per cell count, from ComputeDiscordance.R"),
  make_option("--min-k", type = "integer", default = 2L,
              help = "Restrict panel 1 to sites with at least this many matched cells [default: %default]"),
  make_option("--min-sites", type = "integer", default = 1000L,
              help = "Drop k values backed by fewer sites than this from the mean-D panel [default: %default]"),
  make_option("--max-k-plot", type = "integer", default = 20L,
              help = "Upper k shown in panels 2-3 [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

COLS <- c(Keratinocytes = "#4393c3", MerkelCarcinoma = "#d95f02")
LABS <- c(Keratinocytes = "Keratinocytes", MerkelCarcinoma = "Merkel Cell Carcinoma")

h <- fread(opt$hist)
s <- fread(opt$summary)

# panel 1 -- sites at each attainable value of D, on a log count axis.
# D is discrete given k (only {0,1} at k=2), so these are exact counts at the values
# D can actually take -- no smoothing, and the log axis keeps the small values legible
# next to the ~18M sites sitting at D = 0.
hk <- h[k_match >= opt[["min-k"]]]
dens <- hk[, .(n_sites = sum(n_sites)), by = .(cell_type, D)]
fmt_n <- function(x) ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                     ifelse(x >= 1e3, paste0(x / 1e3, "k"), as.character(x)))
p1 <- ggplot(dens, aes(D, n_sites, colour = cell_type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  scale_y_log10(labels = fmt_n) +
  labs(x = "Pairwise discordance D", y = "Sites (log scale)",
       title = paste0("Discordance distribution (k >= ", opt[["min-k"]], ")")) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
# panels 2 and 3 -- how many cells cover a site. Shown per cell type BEFORE matching:
# after matching k is one number per site and so identical for both types, which is
# why the two lines used to sit exactly on top of each other.
cps <- fread(opt[["cells-per-site"]])[matched == FALSE & k <= opt[["max-k-plot"]]]
cps[, share := n_sites / sum(n_sites), by = cell_type]

p2 <- ggplot(cps, aes(k, share, colour = cell_type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  labs(x = "Cells covering the site", y = "Share of sites",
       title = "Cells per site (share)") +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Log y: the counts span four orders of magnitude, and the collapse at high k is the
# point -- it is what makes the high-k end of panel 4 rest on almost no data.
p3 <- ggplot(cps, aes(k, n_sites, colour = cell_type)) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  scale_y_log10(labels = function(x) ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                                     ifelse(x >= 1e3, paste0(x / 1e3, "k"), as.character(x)))) +
  labs(x = "Cells covering the site", y = "Sites (log scale)",
       title = "Cells per site (count)") +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# panel 3 -- mean D at each k; the comparison without granularity artefacts
# Only k values with enough sites behind them: the tail was dominated by points
# resting on a few hundred sites, which read as a trend rather than as noise.
bk <- s[!is.na(k_match) & k_match <= opt[["max-k-plot"]] & n_sites >= opt[["min-sites"]]]
p4 <- ggplot(bk, aes(k_match, mean_D, colour = cell_type)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  labs(x = "Cells covering the site (matched)", y = "Mean D",
       title = paste0("Mean discordance by cell count (k with \u2265 ",
                      format(opt[["min-sites"]], big.mark = ","), " sites)")) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ov <- s[is.na(k_match)]
for (i in seq_len(nrow(ov)))
  message(sprintf("  %-16s matched=%-5s mean D = %.4f   (mean k = %.2f, %s sites)",
                  ov$cell_type[i], ov$matched[i], ov$mean_D[i], ov$mean_k[i],
                  format(ov$n_sites[i], big.mark = ",")))

ggsave(opt$output, (p1 | p2) / (p3 | p4), width = 11, height = 9, dpi = 150)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), (p1 | p2) / (p3 | p4), width = 11, height = 9, dpi = 150, device = grDevices::svg)
message("Wrote: ", opt$output)
