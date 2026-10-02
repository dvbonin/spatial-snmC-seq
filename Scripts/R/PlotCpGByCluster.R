library(optparse)
library(tidyverse)
library(data.table)

# Per-cluster CpG-coverage boxplot, coloured by well type. Multi-cell wells are
# now part of the clustering, so they get their own category instead of falling
# into an NA well_class. The plate QC threshold line is gone -- QC is no longer
# decided in this pipeline.

option_list <- list(
  make_option("--clusters",   type = "character", help = "Cluster table from LeidenCluster.R"),
  make_option("--well-stats", type = "character", help = "Per-well stats CSV from BuildWellStats.R"),
  make_option("--prefix",     type = "character", help = "<sample>_<plate>"),
  make_option("--output",     type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

ann <- fread(opt$clusters, data.table = FALSE) %>%
  mutate(cluster = factor(cluster)) %>%
  left_join(fread(opt[["well-stats"]], data.table = FALSE) %>%
              select(WellPosition, CellCount, CpGsCovered),
            by = c("well" = "WellPosition")) %>%
  mutate(WellType = case_when(CellCount == 0 ~ "Empty",
                              CellCount == 1 ~ "Singlet",
                              CellCount >  1 ~ "Multi",
                              TRUE           ~ NA_character_))

plot_width <- max(6, 1.2 * nlevels(ann$cluster))
p <- ggplot(ann, aes(x = cluster, y = CpGsCovered)) +
  geom_boxplot(fill = "grey90", outlier.shape = NA) +
  geom_jitter(aes(colour = WellType), width = 0.15, alpha = 0.6, size = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = c(Empty = "#E8998D", Singlet = "#4393c3", Multi = "#2F6B5E"),
                      na.value = "grey70", name = NULL) +
  scale_y_log10() +
  theme_bw(base_size = 13) +
  labs(x = "Leiden cluster", y = "CpGs covered (log scale)",
       title = paste0("CpG coverage by Leiden cluster — ", opt$prefix))
ggsave(opt$output, p, width = plot_width, height = 5.5, dpi = 150)
message("Wrote: ", opt$output)
