library(optparse)
library(tidyverse)
library(umap)
library(data.table)
library(patchwork)

# UMAP parameter sweep (n_neighbors x min_dist) for the actual QC-passed
# clustering, on two populations, each with its own fresh PCA:
#   1. QC-passed single cells only (CellCount == 1 & mean_empty_frac < threshold).
#   2. All single cells (CellCount == 1) -- coloured/shaped by whether they
#      pass, so a fixed PCA doesn't have to be recomputed to see where the
#      removed ones actually sit relative to the kept ones.

option_list <- list(
  make_option("--matrix", type = "character", help = "Merged FirstPrep mean_shrunken_residuals.csv.gz"),
  make_option("--well-empty-frac-csv", type = "character", help = "FirstPrep_Merged_well_empty_frac.csv from SweepLeidenClusterUMAP.R"),
  make_option("--well-stats-p3", type = "character"),
  make_option("--well-stats-p4", type = "character"),
  make_option("--threshold", type = "double", default = 0.1,
              help = "Wells with mean_empty_frac below this pass QC [default: %default]"),
  make_option("--n-neighbors", type = "character", default = "20,25,30"),
  make_option("--min-dist", type = "character", default = "0.25,0.5,0.75"),
  make_option("--output-qc-passed-grid", type = "character"),
  make_option("--output-all-single-grid", type = "character"),
  make_option("--output-qc-passed-pca", type = "character", default = NULL,
              help = "Persist the QC-passed single-cell PCA so downstream figures reuse it")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["output-qc-passed-grid"]]),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt[["output-all-single-grid"]]), recursive = TRUE, showWarnings = FALSE)

K_VALUES  <- as.integer(strsplit(opt[["n-neighbors"]], ",")[[1]])
MD_VALUES <- as.numeric(strsplit(opt[["min-dist"]], ",")[[1]])

well_stats <- bind_rows(
  fread(opt[["well-stats-p3"]], data.table = FALSE) %>% mutate(Plate = "P3"),
  fread(opt[["well-stats-p4"]], data.table = FALSE) %>% mutate(Plate = "P4")
) %>% select(Plate, WellPosition, CellType, CellCount)

empty_frac <- fread(opt[["well-empty-frac-csv"]], data.table = FALSE) %>% select(cell_id, mean_empty_frac)

mat_df <- fread(opt$matrix, data.table = FALSE)
rownames(mat_df) <- mat_df[[1]]
mat_df[[1]] <- NULL
mat_df <- as.matrix(mat_df)

annotation_all <- tibble(cell_id = rownames(mat_df)) %>%
  mutate(Plate = str_extract(cell_id, "(?<=FirstPrep_)P[34](?=_)"),
         well  = str_extract(cell_id, "[A-P][0-9]{2}$")) %>%
  left_join(well_stats, by = c("Plate", "well" = "WellPosition")) %>%
  left_join(empty_frac, by = "cell_id") %>%
  filter(CellCount == 1) %>%
  mutate(passes = mean_empty_frac < opt$threshold)
message("All single cells: ", nrow(annotation_all))

annotation_pass <- annotation_all %>% filter(passes)
message("QC-passed single cells: ", nrow(annotation_pass))

run_pca <- function(annotation) {
  mat_filt <- mat_df[annotation$cell_id, , drop = FALSE]
  keep_cols <- colSums(!is.na(mat_filt)) > 0
  mat_filt  <- mat_filt[, keep_cols, drop = FALSE]
  for (j in seq_len(ncol(mat_filt))) {
    v <- mat_filt[, j]
    if (anyNA(v)) v[is.na(v)] <- median(v, na.rm = TRUE)
    mat_filt[, j] <- v
  }
  keep_cols2 <- apply(mat_filt, 2, function(v) var(v, na.rm = TRUE) > 0)
  mat_filt   <- mat_filt[, keep_cols2, drop = FALSE]

  X   <- scale(mat_filt)
  X   <- X[, colSums(is.nan(X)) == 0, drop = FALSE]
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  n_pcs <- min(30, ncol(pca$x))
  pca$x[, 1:n_pcs, drop = FALSE]
}

message("PCA (QC-passed)...")
pcs_pass <- run_pca(annotation_pass)

# Persisted so figures built on the QC-passed cells (see ../../UMAP) reuse this
# embedding instead of recomputing their own from the 230 x ~76k matrix, which
# would silently drift from the QC decision made here.
if (!is.null(opt[["output-qc-passed-pca"]])) {
  dir.create(dirname(opt[["output-qc-passed-pca"]]), recursive = TRUE, showWarnings = FALSE)
  as_tibble(pcs_pass, rownames = "cell_id") %>%
    mutate(Plate = str_extract(cell_id, "(?<=FirstPrep_)P[34](?=_)"),
           well  = str_extract(cell_id, "[A-P][0-9]{2}$"), .after = cell_id) %>%
    write_csv(opt[["output-qc-passed-pca"]])
  message("Wrote: ", opt[["output-qc-passed-pca"]], " (", nrow(pcs_pass), " cells x ",
          ncol(pcs_pass), " PCs)")
}
message("PCA (all single cells)...")
pcs_all <- run_pca(annotation_all)

CELLTYPE_DISPLAY <- c(HealthyKeratinocytes = "Keratinocytes", MerkelCarcinoma = "Merkel Cell Carcinoma")
cell_type_palette <- c("#4393c3", "#d95f02")

make_grid <- function(annotation, pcs, shape_by_pass, title, output) {
  cell_types <- sort(setdiff(unique(annotation$CellType), NA))
  cell_type_colors <- setNames(cell_type_palette[seq_along(cell_types)], cell_types)
  cell_type_labels <- ifelse(cell_types %in% names(CELLTYPE_DISPLAY), CELLTYPE_DISPLAY[cell_types], cell_types)

  plots <- list()
  for (k in K_VALUES) {
    for (md in MD_VALUES) {
      set.seed(42)
      nb <- min(k, nrow(pcs) - 1)
      um <- umap(pcs, n_neighbors = nb, min_dist = md, metric = "euclidean")
      ann <- annotation %>% mutate(UMAP1 = um$layout[, 1], UMAP2 = um$layout[, 2])

      if (shape_by_pass) {
        p <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType, shape = passes)) +
          geom_point(size = 1.5, alpha = 0.6) +
          scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 4),
                             labels = c(`TRUE` = "Passed", `FALSE` = "Failed"),
                             name = paste0("QC (mean empty-frac < ", opt$threshold, ")"))
      } else {
        p <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType)) +
          geom_point(shape = 16, size = 1.5, alpha = 0.6)
      }
      p <- p +
        scale_colour_manual(values = cell_type_colors, labels = cell_type_labels,
                            name = "Cell type") +
        guides(colour = guide_legend(override.aes = list(size = 3.2)),
               shape  = guide_legend(override.aes = list(size = 3.2))) +
        theme_minimal(base_size = 11, base_family = "Nimbus Sans") +
        theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
              axis.text = element_blank(), axis.ticks = element_blank()) +
        labs(title = paste0("n_neighbors=", k, ", min_dist=", md))
      plots[[paste0(k, "_", md)]] <- p
    }
  }
  ordered <- unlist(lapply(K_VALUES, function(k) lapply(MD_VALUES, function(md) paste0(k, "_", md))))
  grid <- wrap_plots(plots[ordered], ncol = length(MD_VALUES)) +
    plot_annotation(title = title,
                    subtitle = "UMAP n_neighbors (rows) x min_dist (cols)") +
    plot_layout(guides = "collect")
  grid <- grid & theme(legend.position = "bottom",
                       legend.text     = element_text(size = 13),
                       legend.title    = element_text(size = 13.5),
                       legend.key.size = unit(1.1, "lines"))
  ggsave(output, grid, width = 4.5 * length(MD_VALUES),
         height = 4.2 * length(K_VALUES) + 0.7, dpi = 150)
  # SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
  # is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
  if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", output), grid, width = 4.5 * length(MD_VALUES),
         height = 4.2 * length(K_VALUES) + 0.7, dpi = 150, device = grDevices::svg)
  message("Wrote: ", output)
}

make_grid(annotation_pass, pcs_pass, shape_by_pass = FALSE,
          title = "UMAP Sweep: QC-Passed Single Cells",
          output = opt[["output-qc-passed-grid"]])
make_grid(annotation_all,  pcs_all,  shape_by_pass = TRUE,
          title = "UMAP Sweep: All Single Cells",
          output = opt[["output-all-single-grid"]])
