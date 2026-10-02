library(optparse)
library(tidyverse)
library(data.table)
library(patchwork)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# UMAP panel for FirstPrep P3+P4 pooled into ONE MethSCAn run (joint VMR calling,
# see ../../DataPrePrep/Scripts/BuildMergedFirstPrepMethScan.sbatch). Plotted on
# the representative-combo embedding exported by SweepLeidenClusterUMAP.R
# (one real clustering from the sweep, not an independently re-computed UMAP)
# so this panel and the sweep grid agree on what the data actually looks like.
# 6 panels, 3x2: top Cell Type / Plate, middle Cell Count / Coverage, bottom
# Mean Methylation / QC (passes = single-cell & mean_empty_frac < threshold).

option_list <- list(
  make_option("--embedding", type = "character",
              help = "Representative-combo embedding CSV from SweepLeidenClusterUMAP.R (cell_id, well, Plate, cluster, UMAP1, UMAP2)"),
  make_option("--well-empty-frac-csv", type = "character", help = "FirstPrep_Merged_well_empty_frac.csv from SweepLeidenClusterUMAP.R"),
  make_option("--well-stats-p3", type = "character"),
  make_option("--well-stats-p4", type = "character"),
  make_option("--threshold", type = "double", default = 0.1,
              help = "Wells with mean_empty_frac below this pass QC [default: %default]"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

well_stats <- bind_rows(
  fread(opt[["well-stats-p3"]], data.table = FALSE) %>% mutate(Plate = "P3"),
  fread(opt[["well-stats-p4"]], data.table = FALSE) %>% mutate(Plate = "P4")
) %>% select(Plate, WellPosition, CellType, CellCount, CpGsCovered, MeanGlobMeth_CpG)

empty_frac <- fread(opt[["well-empty-frac-csv"]], data.table = FALSE) %>% select(cell_id, mean_empty_frac)

message("Loading representative embedding: ", opt$embedding)
ann <- fread(opt$embedding, data.table = FALSE) %>%
  left_join(well_stats, by = c("Plate", "well" = "WellPosition")) %>%
  left_join(empty_frac, by = "cell_id")
message(nrow(ann), " cells")

# ---- Plotting ----
cell_count_levels <- sort(union(c(0, 1, 2, 5, 10), unique(ann$CellCount)))
non_zero_levels   <- cell_count_levels[cell_count_levels != 0]
cell_count_colors <- c("0" = "#E8998D",
                        setNames(colorRampPalette(c("#89C2D9", "#2F6B5E"))(length(non_zero_levels)),
                                 as.character(non_zero_levels)))

cell_type_palette <- c("#4393c3", "#d95f02", "#7570b3", "#66a61e", "#e6ab02", "#a6761d")
cell_types <- sort(setdiff(unique(ann$CellType), NA))
cell_type_colors <- setNames(cell_type_palette[seq_along(cell_types)], cell_types)
CELLTYPE_DISPLAY <- c(HealthyKeratinocytes = "Keratinocytes", MerkelCarcinoma = "Merkel Cell Carcinoma")
cell_type_labels <- ifelse(cell_types %in% names(CELLTYPE_DISPLAY), CELLTYPE_DISPLAY[cell_types], cell_types)

plate_colors <- c(P3 = "#8B7CC0", P4 = "#D9A441")

base_theme <- theme_minimal(base_size = 13, base_family = FIGURE_FONT) +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_blank(), axis.ticks = element_blank())

p_ct <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_manual(values = cell_type_colors, labels = cell_type_labels, na.value = "grey70", name = "Cell Type") +
  base_theme + labs(title = "Cell Type")

p_plate <- ggplot(ann, aes(UMAP1, UMAP2, colour = Plate)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_manual(values = plate_colors, labels = c(P3 = "P1", P4 = "P2"), name = "Plate") +
  base_theme + labs(title = "Plate")

p_cc <- ggplot(ann, aes(UMAP1, UMAP2, colour = factor(CellCount, levels = cell_count_levels))) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_manual(values = cell_count_colors, drop = FALSE, name = "Cell Count") +
  guides(colour = guide_legend(reverse = TRUE)) +
  base_theme + labs(title = "Cell Count")

p_cov <- ggplot(ann, aes(UMAP1, UMAP2, colour = log10(CpGsCovered))) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_gradientn(colours = unname(cell_count_colors), na.value = "grey70",
                         name = expression(log[10]~"CpGs")) +
  base_theme + labs(title = "Coverage")

p_meanmeth <- ggplot(ann, aes(UMAP1, UMAP2, colour = MeanGlobMeth_CpG)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_gradientn(colours = unname(cell_count_colors), na.value = "grey70", name = "Mean Methylation") +
  base_theme + labs(title = "Mean Methylation")

ann <- ann %>%
  mutate(qc_group = case_when(
    CellCount != 1                              ~ "Empty/Multi",
    mean_empty_frac < opt$threshold              ~ "Passed",
    TRUE                                          ~ "Removed"
  ))
qc_colors <- c(Passed = "#2ca02c", Removed = "#d62728", "Empty/Multi" = "grey70")

p_qc <- ggplot(ann, aes(UMAP1, UMAP2, colour = qc_group)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_manual(values = qc_colors, name = "QC") +
  base_theme + labs(title = paste0("QC Status: Single Cells (mean empty-frac < ",
                                   opt$threshold, ")"))

panel6 <- (p_ct | p_plate) / (p_cc | p_cov) / (p_meanmeth | p_qc)
save_figure( panel6, opt$output, width = 14, height = 18, dpi = 150)
message("Wrote: ", opt$output)
