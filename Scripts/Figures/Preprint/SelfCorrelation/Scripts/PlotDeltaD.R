library(optparse)
library(tidyverse)
library(data.table)
library(patchwork)

# Main result: the genome-wide distribution of dD = D_MCC - D_Kera, paired per CpG.
# Mass to the right of zero means MCC cells disagree with each other more than
# keratinocytes do at the same sites, with the same number of cells on both sides.
#
# The dD-versus-k panel is a sensitivity check, not a result: it shows whether the
# shift holds across coverage levels or is carried by sparsely covered sites. Only k
# values with enough sites behind them are drawn.

option_list <- list(
  make_option("--hist", type = "character", help = "(k, bin, n_sites) from ComputeDeltaD.R"),
  make_option("--summary", type = "character", help = "Overall and per-k statistics"),
  make_option("--min-sites", type = "integer", default = 1000L,
              help = "Drop k values backed by fewer sites than this [default: %default]"),
  make_option("--max-k-plot", type = "integer", default = 20L),
  make_option("--xlim", type = "double", default = 0.6,
              help = "Half-width of the dD axis [default: %default]"),
  make_option("--label", type = "character", default = "",
              help = "Suffix for panel titles, e.g. the min-cells threshold"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

h <- fread(opt$hist)
s <- fread(opt$summary)
ov <- s[is.na(k)]
lab <- if (nchar(opt$label)) paste0(" ", opt$label) else ""

# panel 1 -- the distribution itself
dd <- h[, .(n_sites = sum(n_sites)), by = bin]
dd[, share := n_sites / sum(n_sites)]
p1 <- ggplot(dd[abs(bin) <= opt$xlim], aes(bin, share)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_area(fill = "#4393c3", alpha = 0.35) +
  geom_line(colour = "#1c5cab", linewidth = 0.6) +
  geom_vline(xintercept = ov$mean_dD, colour = "#B2182B", linewidth = 0.6) +
  labs(x = expression(Delta*D~"="~D[MCC]-D[Kera]), y = "Share of sites",
       title = paste0("Paired discordance difference", lab),
       subtitle = sprintf("mean %+.4f (red), median %+.4f, %.1f%% of sites > 0, n = %s",
                          ov$mean_dD, ov$median_dD, 100 * ov$frac_pos,
                          format(ov$n_sites, big.mark = ","))) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(panel.grid.minor = element_blank())

# panel 2 -- same, log y, so the tails are visible rather than flattened
p2 <- ggplot(dd[abs(bin) <= opt$xlim & n_sites > 0], aes(bin, n_sites)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#1c5cab", linewidth = 0.6) +
  scale_y_log10(labels = function(x) ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                                     ifelse(x >= 1e3, paste0(x / 1e3, "k"), as.character(x)))) +
  labs(x = expression(Delta*D), y = "Sites (log scale)",
       title = paste0("Same distribution, log counts", lab)) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(panel.grid.minor = element_blank())

# panel 3 -- sensitivity: is the shift consistent across coverage?
bk <- s[!is.na(k) & k <= opt[["max-k-plot"]] & n_sites >= opt[["min-sites"]]]
p3 <- ggplot(bk, aes(k, mean_dD)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = ov$mean_dD, colour = "#B2182B", linewidth = 0.4, alpha = 0.6) +
  geom_line(colour = "#1c5cab", linewidth = 0.8) +
  labs(x = "Matched cells per site (k)", y = expression("mean"~Delta*D),
       title = paste0("Sensitivity to coverage", lab),
       subtitle = paste0("k with >= ", format(opt[["min-sites"]], big.mark = ","),
                         " sites; red = genome-wide mean")) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(panel.grid.minor = element_blank())

# panel 4 -- how many sites each k contributes, so panel 3 can be weighted by eye
p4 <- ggplot(s[!is.na(k) & k <= opt[["max-k-plot"]]], aes(k, n_sites)) +
  geom_line(colour = "#1c5cab", linewidth = 0.8) +
  scale_y_log10(labels = function(x) ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                                     ifelse(x >= 1e3, paste0(x / 1e3, "k"), as.character(x)))) +
  labs(x = "Matched cells per site (k)", y = "Sites (log scale)",
       title = paste0("Sites per coverage level", lab)) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(panel.grid.minor = element_blank())

message(sprintf("  mean dD = %+.5f   median = %+.5f   frac>0 = %.4f   frac<0 = %.4f   n = %s",
                ov$mean_dD, ov$median_dD, ov$frac_pos, ov$frac_neg,
                format(ov$n_sites, big.mark = ",")))
ggsave(opt$output, (p1 | p2) / (p3 | p4), width = 12, height = 9, dpi = 150)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), (p1 | p2) / (p3 | p4), width = 12, height = 9, dpi = 150, device = grDevices::svg)
message("Wrote: ", opt$output)
