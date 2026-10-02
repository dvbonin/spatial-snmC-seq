library(optparse)
library(tidyverse)
library(data.table)
library(umap)
library(ggtext)

# UMAP scatter coloured by Leiden cluster. The embedding is visualisation only --
# the clustering in LeidenCluster.R never uses these coordinates. Legend shows
# each cluster's #Empty wells out of its total; the empty-heaviest cluster is bold.

option_list <- list(
  make_option("--pca",        type = "character", help = "PCA CSV.gz from ComputePCA.R"),
  make_option("--clusters",   type = "character", help = "Cluster table from LeidenCluster.R"),
  make_option("--well-stats", type = "character", help = "Per-well stats CSV from BuildWellStats.R"),
  make_option("--prefix",     type = "character", help = "<sample>_<plate>"),
  make_option("--k-neighbors", type = "integer", default = 15L),
  make_option("--output",     type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

pca_df   <- fread(opt$pca, data.table = FALSE)
pc_cols  <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
pcs      <- as.matrix(pca_df[, pc_cols])
clusters <- fread(opt$clusters, data.table = FALSE)
well_stats <- fread(opt[["well-stats"]], data.table = FALSE) %>%
  select(WellPosition, CellCount, well_class)

set.seed(42)
nb <- min(opt[["k-neighbors"]], nrow(pcs) - 1)
um <- umap(pcs, n_neighbors = nb, min_dist = 0.3, metric = "euclidean")

ann <- clusters %>%
  mutate(cluster = factor(cluster), UMAP1 = um$layout[, 1], UMAP2 = um$layout[, 2]) %>%
  left_join(well_stats, by = c("well" = "WellPosition"))

cluster_stats <- ann %>%
  group_by(cluster) %>%
  summarise(empty_n = sum(CellCount == 0, na.rm = TRUE), n_wells = n(), .groups = "drop")
empty_cluster_id <- cluster_stats %>% slice_max(empty_n, n = 1, with_ties = FALSE) %>% pull(cluster)

legend_labels <- cluster_stats %>%
  mutate(cluster = as.character(cluster),
         label = if_else(cluster == as.character(empty_cluster_id),
                         paste0(cluster, " (**", empty_n, "/", n_wells, "**)"),
                         paste0(cluster, " (", empty_n, "/", n_wells, ")"))) %>%
  select(cluster, label) %>% deframe()

p <- ggplot(ann, aes(UMAP1, UMAP2, colour = cluster)) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_discrete(labels = legend_labels, name = "Leiden cluster\n(Empty / total wells)") +
  theme_minimal(base_size = 13) +
  theme(legend.text = element_markdown()) +
  labs(title = paste0("UMAP coloured by Leiden cluster — ", opt$prefix))
ggsave(opt$output, p, width = 9, height = 6, dpi = 150)
message("Wrote: ", opt$output)
