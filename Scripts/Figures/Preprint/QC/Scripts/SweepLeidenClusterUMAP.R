library(optparse)
library(tidyverse)
library(data.table)
library(igraph)
library(umap)
library(patchwork)

# Leiden clustering sweep over n_neighbors x resolution, on the all-wells PCA
# from ComputeMethScanPCA.R, to see how stable the cluster structure is
# and to build a QC signal that doesn't depend on any one clustering's exact
# boundaries. Three outputs:
#   1. A UMAP grid (rows = n_neighbors, cols = resolution), one clustering per
#      panel, open circles.
#   2. Per-well mean fraction-empty-in-own-cluster, averaged across every
#      (n_neighbors, resolution) combination in the sweep -- a well that
#      consistently lands in a mostly-empty cluster is a QC-worthy well
#      regardless of any single clustering's exact boundaries. Written as a
#      CSV; plotted separately by PlotEmptyFractionBeeswarm.py.
#   3. The (--representative-k, --representative-resolution) combination's own
#      UMAP embedding + cluster, exported for PlotMergedFirstPrepPanel6.R to
#      plot on -- so panel6 shows one real clustering from the sweep rather
#      than yet another independently-parameterized embedding.
# The pairwise distance matrix (independent of n_neighbors) is computed once
# and reused across the whole sweep.

option_list <- list(
  make_option("--pca", type = "character", help = "All-wells PCA CSV.gz from ComputeMethScanPCA.R"),
  make_option("--well-stats-p3", type = "character"),
  make_option("--well-stats-p4", type = "character"),
  make_option("--n-neighbors", type = "character", default = "15,20,25"),
  make_option("--resolution", type = "character", default = "0.25,0.5,0.75"),
  make_option("--objective", type = "character", default = "modularity"),
  make_option("--representative-k", type = "integer", default = 25L),
  make_option("--representative-resolution", type = "double", default = 0.5),
  make_option("--output-umap-grid", type = "character"),
  make_option("--output-well-empty-frac-csv", type = "character"),
  make_option("--output-representative-embedding-csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["output-umap-grid"]]),                    recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt[["output-well-empty-frac-csv"]]),          recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt[["output-representative-embedding-csv"]]), recursive = TRUE, showWarnings = FALSE)

K_VALUES   <- as.integer(strsplit(opt[["n-neighbors"]], ",")[[1]])
RES_VALUES <- as.numeric(strsplit(opt$resolution, ",")[[1]])

well_stats <- bind_rows(
  fread(opt[["well-stats-p3"]], data.table = FALSE) %>% mutate(Plate = "P3"),
  fread(opt[["well-stats-p4"]], data.table = FALSE) %>% mutate(Plate = "P4")
) %>% select(Plate, WellPosition, CellType, CellCount)

message("Loading cached PCA: ", opt$pca)
pca_df  <- fread(opt$pca, data.table = FALSE)
pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)

annotation <- pca_df %>%
  select(cell_id, well, all_of(pc_cols)) %>%
  mutate(Plate = str_extract(cell_id, "(?<=FirstPrep_)P[34](?=_)")) %>%
  left_join(well_stats, by = c("Plate", "well" = "WellPosition"))
pcs <- as.matrix(annotation[, pc_cols])
rownames(pcs) <- annotation$cell_id
message(nrow(pcs), " cells x ", ncol(pcs), " PCs")

# ---- SNN graph builder (Jaccard-weighted); the distance matrix doesn't ----
# ---- depend on k, so compute it once and reuse                        ----
d <- as.matrix(dist(pcs))
diag(d) <- Inf
n <- nrow(d)

build_snn_graph <- function(k) {
  knn_idx <- t(apply(d, 1, function(row) order(row)[1:k]))
  neighbor_sets <- lapply(seq_len(n), function(i) knn_idx[i, ])

  el <- combn(n, 2)
  is_edge <- apply(el, 2, function(pair) {
    (pair[2] %in% neighbor_sets[[pair[1]]]) || (pair[1] %in% neighbor_sets[[pair[2]]])
  })
  el <- el[, is_edge, drop = FALSE]
  jac <- apply(el, 2, function(pair) {
    length(intersect(neighbor_sets[[pair[1]]], neighbor_sets[[pair[2]]])) /
      length(union(neighbor_sets[[pair[1]]], neighbor_sets[[pair[2]]]))
  })
  el  <- el[, jac > 0, drop = FALSE]
  jac <- jac[jac > 0]

  g <- graph_from_edgelist(t(el), directed = FALSE)
  E(g)$weight <- jac
  g
}

representative_ann <- NULL
empty_frac_long <- list()
plots <- list()

for (k in K_VALUES) {
  message("n_neighbors = ", k)
  g  <- build_snn_graph(k)
  set.seed(42)
  um <- umap(pcs, n_neighbors = k, min_dist = 0.3, metric = "euclidean")

  for (res in RES_VALUES) {
    set.seed(42)
    cl <- cluster_leiden(g, objective_function = opt$objective, resolution = res,
                          weights = E(g)$weight, n_iterations = 2)
    ann <- annotation %>%
      mutate(cluster = factor(membership(cl)), UMAP1 = um$layout[, 1], UMAP2 = um$layout[, 2])

    cluster_empty_frac <- ann %>%
      group_by(cluster) %>%
      summarise(empty_frac = mean(CellCount == 0, na.rm = TRUE), .groups = "drop")
    ann <- ann %>% left_join(cluster_empty_frac, by = "cluster")

    empty_frac_long[[length(empty_frac_long) + 1]] <-
      ann %>% select(cell_id, Plate, well, CellCount, empty_frac) %>%
      mutate(k_neighbors = k, resolution = res)

    if (k == opt[["representative-k"]] && res == opt[["representative-resolution"]]) {
      representative_ann <- ann
    }

    p <- ggplot(ann, aes(UMAP1, UMAP2, colour = cluster)) +
      geom_point(shape = 1, size = 1.8, alpha = 0.6) +
      theme_minimal(base_size = 11, base_family = "Nimbus Sans") +
      theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
            axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none") +
      labs(title = paste0("n_neighbors=", k, ", resolution=", res))
    plots[[paste0(k, "_", res)]] <- p
  }
}

if (is.null(representative_ann)) {
  stop("representative-k/-resolution (", opt[["representative-k"]], ", ", opt[["representative-resolution"]],
       ") is not part of the swept grid -- add it to --n-neighbors/--resolution or pick a value from there.")
}
write_csv(representative_ann %>% select(cell_id, well, Plate, cluster, UMAP1, UMAP2),
          opt[["output-representative-embedding-csv"]])
message("Wrote: ", opt[["output-representative-embedding-csv"]])

ordered <- unlist(lapply(K_VALUES, function(k) lapply(RES_VALUES, function(res) paste0(k, "_", res))))
grid <- wrap_plots(plots[ordered], ncol = length(RES_VALUES)) +
  plot_annotation(title = "Leiden Sweep: All Wells",
                  subtitle = "n_neighbors (rows) x resolution (cols); colour = Leiden cluster")
n_row <- length(K_VALUES); n_col <- length(RES_VALUES)
ggsave(opt[["output-umap-grid"]], grid, width = 4.5 * n_col, height = 4.2 * n_row, dpi = 150)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if (nzchar(Sys.getenv("SNMC_SVG"))) ggsave( sub("\\.png$", ".svg", opt[["output-umap-grid"]]), grid, width = 4.5 * n_col, height = 4.2 * n_row, dpi = 150, device = grDevices::svg)
message("Wrote: ", opt[["output-umap-grid"]])

# ---- Per-well mean fraction-empty-in-own-cluster across the whole sweep ----
empty_frac_df <- bind_rows(empty_frac_long)
well_rank <- empty_frac_df %>%
  group_by(cell_id, Plate, well, CellCount) %>%
  summarise(mean_empty_frac = mean(empty_frac), .groups = "drop")

write_csv(well_rank, opt[["output-well-empty-frac-csv"]])
message("Wrote: ", opt[["output-well-empty-frac-csv"]])
