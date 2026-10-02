library(optparse)
library(tidyverse)
library(data.table)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# Renders a pairwise table as a heatmap, faceted by cell type so that within-type
# comparisons sit in the diagonal blocks and cross-type ones off it.
#
# --value r  the correlation. Scale fixed to -1..1 rather than scaled to the data, so
#            every figure in this directory stays comparable: a colour means the same
#            correlation in all of them.
# --value n  how many sites the correlation rested on. Scaled to its own range, since
#            there is no meaningful absolute reference. Diagonal cells give each
#            profile its own site count; off-diagonal cells give the intersection,
#            which is what the correlation actually used.
#
# Cells with no value (Public exists only for keratinocytes) render grey.

option_list <- list(
  make_option("--input", type = "character", help = "Pairs CSV from BuildCorrelationMatrix.R"),
  make_option("--title", type = "character", default = "Pairwise Methylation Correlation"),
  make_option("--value", type = "character", default = "r",
              help = "'r' (correlation) or 'n' (sites used) [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

PREFERRED_ORDER <- c("Bulk_KZ", "Pseudobulk_KZ", "Public_KZ", "Bulk_MCC", "Pseudobulk_MCC")
CELL_LABELS <- c(KZ = "Keratinocytes", MCC = "Merkel Cell Carcinoma")
SOURCE_LABELS <- c(Bulk = "Bulk", Pseudobulk = "Pseudobulk", Public = "Public")

pairs <- fread(opt$input, data.table = FALSE)
present <- unique(c(pairs$source1, pairs$source2))
SOURCE_ORDER <- PREFERRED_ORDER[PREFERRED_ORDER %in% present]
SHORT_ORDER  <- unique(sub("_(KZ|MCC)$", "", SOURCE_ORDER))

# The table holds each pair once; mirroring fills the other triangle. distinct() keeps
# the diagonal from being duplicated by its own mirror image.
full <- bind_rows(pairs, pairs %>% rename(source1 = source2, source2 = source1)) %>%
  distinct(source1, source2, .keep_all = TRUE) %>%
  mutate(source1_cell = factor(sub(".*_", "", source1), levels = c("KZ", "MCC")),
         source2_cell = factor(sub(".*_", "", source2), levels = c("KZ", "MCC")),
         source1 = factor(sub("_(KZ|MCC)$", "", source1), levels = SHORT_ORDER),
         source2 = factor(sub("_(KZ|MCC)$", "", source2), levels = rev(SHORT_ORDER)))

fmt_count <- function(x) {
  ifelse(x >= 1e6, paste0(signif(x / 1e6, 3), "M"),
  ifelse(x >= 1e3, paste0(signif(x / 1e3, 3), "k"), as.character(signif(x, 3))))
}

p <- if (opt$value == "n") {
  ggplot(full, aes(source1, source2, fill = n)) +
    geom_tile() +
    geom_text(aes(label = fmt_count(n)), size = 4.6) +
    scale_fill_gradient(low = "#EFF6FB", high = "#2171B5", na.value = "grey90",
                        labels = fmt_count, name = "Sites")
} else {
  ggplot(full, aes(source1, source2, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)), size = 5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                        limits = c(-1, 1), na.value = "grey90", name = "r")
}
p <- p +
  scale_x_discrete(labels = SOURCE_LABELS) +
  scale_y_discrete(labels = SOURCE_LABELS) +
  facet_grid(source2_cell ~ source1_cell,
             labeller = labeller(source1_cell = CELL_LABELS, source2_cell = CELL_LABELS),
             scales = "free", space = "free") +
  coord_equal() +
  theme_minimal(base_size = 13, base_family = FIGURE_FONT) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.15, "lines")) +
  labs(x = NULL, y = NULL, title = str_wrap(opt$title, 45))

n <- length(SOURCE_ORDER)
save_figure( p, opt$output, width = 1.3 * n + 1.5, height = 1.3 * n + 0.8, dpi = 150)
message("Wrote: ", opt$output)
