library(optparse)
library(tidyverse)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# Empty-vs-single-cell read-count densities, one figure per plate group with the
# group's plates pooled. Only wells the layout says hold 0 or 1 cells are shown;
# no QC filtering, since the point is how well read count alone separates the two.
#
# The x axis is deliberately NOT shared between plate groups: the preps differ in
# median yield by roughly 10x, so a common axis would squash the deeper one. The
# groups are therefore comparable within a panel, not across panels.
#
# geom_density normalises each curve to unit area, so the two groups look equally
# populated however different their sizes are. The n is put in the legend for
# that reason.
#
# ReadCounts_Main counts primary-chromosome reads only: chrM and the unplaced
# scaffolds are routed to a separate Alt stream at the well split. That stream
# would add well under 0.1% here, so the axis is unaffected in practice.

option_list <- list(
  make_option("--well-stats-dir", type = "character",
              default = file.path(Sys.getenv("SNMC_PIPE"), "Analysis_Combined/Results/WellStats")),
  make_option("--plates", type = "character",
              help = "Comma-separated <sample>_<plate> prefixes to pool, e.g. FirstPrep_P3,FirstPrep_P4"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

WELL_TYPE_COLORS <- c(Empty = "#e34948", "Single Cell" = "#3aa657")

fmt_count <- function(x) {
  ifelse(x >= 1e9, paste0(signif(x / 1e9, 3), "B"),
  ifelse(x >= 1e6, paste0(signif(x / 1e6, 3), "M"),
  ifelse(x >= 1e3, paste0(signif(x / 1e3, 3), "k"),
                   as.character(signif(x, 3)))))
}

prefixes <- strsplit(opt$plates, ",")[[1]]
well_stats <- map_dfr(prefixes, function(prefix) {
  read_csv(file.path(opt[["well-stats-dir"]], paste0(prefix, "_well_stats.csv")),
           show_col_types = FALSE) %>% mutate(prefix = prefix)
}) %>%
  filter(CellCount %in% c(0, 1)) %>%
  mutate(WellType = if_else(CellCount == 0, "Empty", "Single Cell"))

counts <- well_stats %>% count(WellType)
print(well_stats %>% group_by(WellType) %>%
        summarise(n = n(), min = min(ReadCounts_Main), median = median(ReadCounts_Main),
                  max = max(ReadCounts_Main), .groups = "drop"))

# n folded into the legend key, since the densities themselves carry no size
legend_labels <- counts %>%
  transmute(WellType, label = paste0(WellType, " (n = ", n, ")")) %>%
  deframe()

p <- ggplot(well_stats, aes(x = ReadCounts_Main, fill = WellType, color = WellType)) +
  geom_density(alpha = 0.4, linewidth = 0.8, na.rm = TRUE) +
  scale_x_log10(labels = fmt_count) +
  scale_fill_manual(values = WELL_TYPE_COLORS, labels = legend_labels) +
  scale_color_manual(values = WELL_TYPE_COLORS, labels = legend_labels) +
  theme_bw(base_size = 13, base_family = FIGURE_FONT) +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)) +
  labs(x = "Reads (log scale)", y = "Density", fill = NULL, color = NULL,
       title = "Read counts - Empty vs. Single cell wells")

save_figure( p, opt$output, width = 7, height = 5, dpi = 150)
message("Wrote: ", opt$output, " (", nrow(well_stats), " wells from ",
        paste(prefixes, collapse = " + "), ")")
