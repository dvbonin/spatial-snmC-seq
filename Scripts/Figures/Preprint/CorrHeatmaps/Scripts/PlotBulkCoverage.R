library(optparse)
library(tidyverse)
library(data.table)
library(patchwork)

# Per-site read depth in the two bulk profiles -- the data the top-1% site selection
# is drawn from, and therefore what decides whether that selection is meaningful.
#
# Two panels:
#   left   how many sites sit at each depth
#   right  how many survive a given minimum, which is the number that matters: a
#          site needs enough reads for its beta to be an estimate rather than a coin
#          flip. At 1-2 reads a beta can only be 0, 0.5 or 1, so |KZ - MCC| reaches
#          its ceiling by sampling accident and a top-percentile cut fills up with
#          the least reliable sites in the genome.
#
# Reads a small pre-tabulated table rather than the 28M-row beds, so the figure can
# be restyled without touching the data.

option_list <- list(
  make_option("--input", type = "character", help = "Tabulated counts: sample,coverage,n_sites"),
  make_option("--mark", type = "character", default = "3,5",
              help = "Thresholds to mark on the cumulative panel [default: %default]"),
  make_option("--max-cov", type = "integer", default = 30L,
              help = "Depths above this are pooled into the last bin [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

COLS <- c(Bulk_KZ = "#4393c3", Bulk_MCC = "#d95f02")
LABS <- c(Bulk_KZ = "Keratinocytes (bulk)", Bulk_MCC = "Merkel Cell Carcinoma (bulk)")
marks <- as.integer(strsplit(opt$mark, ",")[[1]])
cap <- opt[["max-cov"]]

d <- fread(opt$input)
tot <- d[, .(total = sum(n_sites)), by = sample]

hist_d <- d[, .(n_sites = sum(n_sites)), by = .(sample, cov = pmin(coverage, cap))]
p_hist <- ggplot(hist_d, aes(cov, n_sites / 1e6, colour = sample)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  scale_x_continuous(breaks = c(1, seq(5, cap, 5)),
                     labels = c("1", as.character(seq(5, cap - 5, 5)), paste0(cap, "+"))) +
  labs(x = "Reads at site", y = "Sites (millions)",
       title = "Per-site read depth") +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Sites retained at each minimum depth, as a fraction of that profile's covered sites
cum_d <- d[order(coverage), .(coverage, kept = rev(cumsum(rev(n_sites)))), by = sample]
cum_d <- merge(cum_d, tot, by = "sample")[coverage <= cap]
p_cum <- ggplot(cum_d, aes(coverage, 100 * kept / total, colour = sample)) +
  geom_vline(xintercept = marks, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = COLS, labels = LABS, name = NULL) +
  scale_x_continuous(breaks = sort(unique(c(1, marks, seq(10, cap, 10))))) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Minimum reads required", y = "Sites retained (%)",
       title = paste0("Retained by depth cutoff (dashed: ",
                      paste(marks, collapse = ", "), ")")) +
  theme_bw(base_size = 12, base_family = "Nimbus Sans") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

for (s in tot$sample) {
  k <- cum_d[sample == s]
  msg <- vapply(marks, function(m) {
    r <- k[coverage == m]
    if (nrow(r)) sprintf(">=%d: %.1f%%", m, 100 * r$kept / r$total) else sprintf(">=%d: -", m)
  }, character(1))
  message("  ", s, "  ", format(tot[sample == s]$total, big.mark = ","),
          " sites   ", paste(msg, collapse = "   "))
}

ggsave(opt$output, p_hist | p_cum, width = 11, height = 4.6, dpi = 150)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt$output), p_hist | p_cum, width = 11, height = 4.6, dpi = 150, device = grDevices::svg)
message("Wrote: ", opt$output)
