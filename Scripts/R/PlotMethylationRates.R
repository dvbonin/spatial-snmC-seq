library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--well-stats", type = "character",
              help = "Per-well stats CSV (from BuildWellStats.R) -- CellCount, MeanGlobMeth_CpG, MeanGlobMeth_CHG, MeanGlobMeth_CHH"),
  make_option("--prefix", type = "character",
              help = "Cell ID prefix: <sample>_<plate>"),
  make_option("--output", type = "character",
              help = "Output path for the methylation-rates boxplot PNG (one box per context, all wells pooled)")
)
opt <- parse_args(OptionParser(option_list = option_list))

prefix      <- opt$prefix
output_path <- opt$output
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

annotation <- read_csv(opt[["well-stats"]], show_col_types = FALSE)

pd <- annotation %>%
  select(WellPosition, CellCount, CpG = MeanGlobMeth_CpG, CHG = MeanGlobMeth_CHG, CHH = MeanGlobMeth_CHH) %>%
  pivot_longer(c(CpG, CHG, CHH), names_to = "context", values_to = "mean_meth") %>%
  mutate(context = factor(context, levels = c("CpG", "CHG", "CHH")))

# Same discrete CellCount palette as the pipeline's clustering panel: 0 = red, non-zero counts
# get a lightblue -> darkgreen gradient.
cell_count_levels <- sort(unique(pd$CellCount))
non_zero_levels   <- cell_count_levels[cell_count_levels != 0]
cell_count_colors <- c("0" = "red",
                        setNames(colorRampPalette(c("lightblue", "darkgreen"))(length(non_zero_levels)),
                                 as.character(non_zero_levels)))

medians <- pd %>%
  group_by(context) %>%
  summarise(median_meth = median(mean_meth, na.rm = TRUE), .groups = "drop")

p <- ggplot(pd, aes(x = context, y = mean_meth)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(aes(colour = factor(CellCount, levels = cell_count_levels)),
              width = 0.15, alpha = 0.6, size = 0.9, na.rm = TRUE) +
  geom_text(data = medians, aes(x = context, y = 1, label = sprintf("%.3f", median_meth)),
            vjust = 0, size = 4, inherit.aes = FALSE) +
  scale_colour_manual(values = cell_count_colors, drop = FALSE, name = "Cells in well") +
  theme_bw(base_size = 13) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(x = NULL, y = "Mean methylation (site-averaged)",
       title = paste0("Methylation rates by context (all wells) — ", prefix))

ggsave(output_path, p, width = 7, height = 5.5, dpi = 150)
message("Wrote: ", output_path)
