library(optparse)
library(tidyverse)
library(data.table)
library(igraph)

# Leiden clustering on a cached PCA (see ComputePCA.R). Emits cluster labels
# only -- the UMAP scatter and the per-cluster coverage boxplot are separate
# rules, and the per-cluster pass/fail QC this rule used to apply is gone: QC is
# decided downstream at the figure level, not here.

option_list <- list(
  make_option("--pca",          type = "character", help = "PCA CSV.gz from ComputePCA.R"),
  make_option("--prefix",       type = "character", help = "Cell ID prefix: <sample>_<plate>"),
  make_option("--clusters-out", type = "character", help = "Output per-well cluster table CSV"),
  make_option("--k-neighbors",  type = "integer", default = 15L,
              help = "Neighbours used to build the SNN graph [default: %default]"),
  make_option("--resolution",   type = "double", default = 1.0,
              help = "Leiden resolution: higher = more, smaller clusters [default: %default]"),
  make_option("--objective",    type = "character", default = "modularity",
              help = "Leiden objective function: modularity or CPM [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["clusters-out"]]), recursive = TRUE, showWarnings = FALSE)

pca_df  <- fread(opt$pca, data.table = FALSE)
pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
pcs     <- as.matrix(pca_df[, pc_cols])
rownames(pcs) <- pca_df$cell_id
message(nrow(pcs), " cells x ", length(pc_cols), " PCs")

# Shared-nearest-neighbour graph (Jaccard-weighted), the construction Seurat and
# scanpy use ahead of Leiden.
build_snn_graph <- function(pcs, k) {
  d <- as.matrix(dist(pcs))
  diag(d) <- Inf
  n <- nrow(d)
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

  # Build on a fixed vertex set: deriving vertices from the edge list drops any
  # cell whose edges all had Jaccard 0, silently misaligning membership().
  g <- make_empty_graph(n = n, directed = FALSE) %>% add_edges(as.vector(el))
  E(g)$weight <- jac
  g
}

g  <- build_snn_graph(pcs, k = opt[["k-neighbors"]])
cl <- cluster_leiden(g, objective_function = opt$objective, resolution = opt$resolution,
                     weights = E(g)$weight, n_iterations = 2)

clusters <- pca_df %>%
  select(cell_id, well) %>%
  mutate(cluster = factor(membership(cl)), sample_plate = opt$prefix)
message("Leiden (k=", opt[["k-neighbors"]], ", resolution=", opt$resolution, ", ",
        opt$objective, "): ", nlevels(clusters$cluster), " clusters")

write_csv(clusters, opt[["clusters-out"]])
message("Wrote: ", opt[["clusters-out"]])
