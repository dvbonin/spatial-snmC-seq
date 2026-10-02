library(optparse)
library(tidyverse)
library(data.table)
library(umap)
library(patchwork)

# Four-panel UMAP for one (sample, plate): cell type, cell count, mean global CpG
# methylation, CpG coverage. Built on the cached PCA from ComputePCA.R -- it
# used to recompute its own PCA for one well-filter variant and reuse the cached
# one for the other, so the two panels were not strictly comparable.
#
# The old fourth panel showed pass/fail against the plate's fixed CpG-coverage
# threshold. That QC no longer exists in this pipeline, and mean methylation is
# the more informative thing to see next to coverage.

option_list <- list(
  make_option("--pca",        type = "character", help = "PCA CSV.gz from ComputePCA.R"),
  make_option("--well-stats", type = "character", help = "Per-well stats CSV from BuildWellStats.R"),
  make_option("--sample",     type = "character"),
  make_option("--plate",      type = "character"),
  make_option("--k-neighbors", type = "integer", default = 15L),
  make_option("--output",     type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

pca_df  <- fread(opt$pca, data.table = FALSE)
pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
pcs     <- as.matrix(pca_df[, pc_cols])

well_stats <- fread(opt[["well-stats"]], data.table = FALSE) %>%
  select(WellPosition, CellType, CellCount, CpGsCovered, MeanGlobMeth_CpG)

set.seed(42)
nb <- min(opt[["k-neighbors"]], nrow(pcs) - 1)
um <- umap(pcs, n_neighbors = nb, min_dist = 0.3, metric = "euclidean")

ann <- pca_df %>%
  select(cell_id, well) %>%
  mutate(UMAP1 = um$layout[, 1], UMAP2 = um$layout[, 2]) %>%
  left_join(well_stats, by = c("well" = "WellPosition"))

cell_count_levels <- sort(union(c(0, 1, 2, 5, 10), unique(ann$CellCount)))
non_zero_levels   <- cell_count_levels[cell_count_levels != 0]
cell_count_colors <- c("0" = "#E8998D",
                       setNames(colorRampPalette(c("#89C2D9", "#2F6B5E"))(length(non_zero_levels)),
                                as.character(non_zero_levels)))

cell_type_palette <- c("#4393c3", "#d95f02", "#7570b3", "#66a61e", "#e6ab02", "#a6761d")
cell_types <- sort(setdiff(unique(ann$CellType), NA))
cell_type_colors <- setNames(cell_type_palette[seq_along(cell_types)], cell_types)

title_prefix <- paste0(opt$sample, " - ", opt$plate, " - ")
base_theme <- theme_minimal(base_size = 13)

p_ct <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = cell_type_colors, na.value = "grey70", name = "Cell type") +
  base_theme + labs(title = paste0(title_prefix, "Cell type"))

p_cc <- ggplot(ann, aes(UMAP1, UMAP2, colour = factor(CellCount, levels = cell_count_levels))) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = cell_count_colors, drop = FALSE, name = "Cell count") +
  guides(colour = guide_legend(reverse = TRUE)) +
  base_theme + labs(title = paste0(title_prefix, "Cell count"))

p_cov <- ggplot(ann, aes(UMAP1, UMAP2, colour = log10(CpGsCovered))) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_gradientn(colours = unname(cell_count_colors), na.value = "grey70",
                         name = expression(log[10]~"CpGs covered")) +
  base_theme + labs(title = paste0(title_prefix, "CpG coverage"))

p_meth <- ggplot(ann, aes(UMAP1, UMAP2, colour = MeanGlobMeth_CpG)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_gradientn(colours = unname(cell_count_colors), na.value = "grey70",
                         name = "Mean mCpG") +
  base_theme + labs(title = paste0(title_prefix, "Mean global CpG methylation"))

ggsave(opt$output, (p_ct | p_cc) / (p_meth | p_cov), width = 14, height = 12, dpi = 150)
message("Done. Plot written to ", opt$output)
