library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--well-stats", type = "character",
              help = "Per-well stats CSV (from BuildWellStats.R) -- CellType, CellCount, CpGsCovered"),
  make_option("--prefix", type = "character",
              help = "Cell ID prefix: <sample>_<plate>"),
  make_option("--output", type = "character",
              help = "Output path for CpG coverage boxplot PNG")
)
opt <- parse_args(OptionParser(option_list = option_list))

prefix      <- opt$prefix
output_path <- opt$output
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

annotation <- read_csv(opt[["well-stats"]], show_col_types = FALSE)

type_order <- sort(unique(annotation$CellType))
pd <- annotation %>%
  mutate(composition = paste0(CellType, "\nn=", CellCount))
level_order <- pd %>%
  distinct(CellType, CellCount, composition) %>%
  mutate(CellType  = factor(CellType, levels = type_order),
         CellCount = as.integer(CellCount)) %>%
  arrange(CellType, CellCount) %>%
  pull(composition)
pd <- pd %>%
  mutate(composition = factor(composition, levels = level_order),
         CellType    = factor(CellType,    levels = type_order))

n_groups   <- n_distinct(paste(annotation$CellType, annotation$CellCount))
plot_width <- max(7, 1.4 * n_groups)

fmt_count <- function(x) {
  ifelse(x >= 1e9, paste0(signif(x / 1e9, 3), "B"),
  ifelse(x >= 1e6, paste0(signif(x / 1e6, 3), "M"),
  ifelse(x >= 1e3, paste0(signif(x / 1e3, 3), "k"),
                   as.character(signif(x, 3)))))
}

p <- ggplot(pd, aes(x = composition, y = CpGsCovered, fill = CellType)) +
  geom_boxplot(outlier.size = 0.6, na.rm = TRUE) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.6, na.rm = TRUE) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  scale_x_discrete(labels = function(x) sub("^.*\n", "", x)) +
  scale_y_log10(labels = fmt_count) +
  labs(x = "Cells in Well", y = "Genome-wide covered CpGs (log scale)",
       title = paste0("CpG coverage by composition — ", prefix))

ggsave(output_path, p, width = plot_width, height = 5, dpi = 150)
message("Wrote: ", output_path)
