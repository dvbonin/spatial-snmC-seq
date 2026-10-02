library(optparse)
library(tidyverse)
library(data.table)

# PCA on a cells x features matrix (MethSCAn residuals or the 100kb-bin matrix).
# Split out of LeidenCluster.R so the clustering, the UMAP plot and the
# clustering panel all read one persisted PCA instead of recomputing their own.
#
# All wells are kept: Empty, single-cell and multi-cell alike. Restricting to
# wells with a known well_class silently dropped every multi-cell well, which is
# a population choice rather than a QC one, and the downstream figure work looks
# at the full plate.

option_list <- list(
  make_option("--matrix",    type = "character",
              help = "mean_shrunken_residuals.csv.gz or mean_meth_100k.csv.gz"),
  make_option("--min-features", type = "integer", default = 1L,
              help = "Minimum covered features per cell to include [default: %default]"),
  make_option("--n-pcs",     type = "integer", default = 30L,
              help = "Number of principal components to keep [default: %default]"),
  make_option("--pca-out",   type = "character",
              help = "Output per-well PCA coordinates CSV.gz (cell_id, well, PC1..PCn)")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["pca-out"]]), recursive = TRUE, showWarnings = FALSE)

message("Loading matrix: ", opt$matrix)
mat_df <- fread(opt$matrix, data.table = FALSE)
rownames(mat_df) <- mat_df[[1]]
mat_df[[1]] <- NULL
mat_df <- as.matrix(mat_df)

annotation <- tibble(cell_id = rownames(mat_df)) %>%
  mutate(well       = str_extract(cell_id, "[A-P][0-9]{2}$"),
         n_features = rowSums(!is.na(mat_df))) %>%
  filter(n_features >= opt[["min-features"]])
message("Keeping ", nrow(annotation), " wells (all classes: Empty / single / multi-cell)")

mat_filt  <- mat_df[annotation$cell_id, , drop = FALSE]
keep_cols <- colSums(!is.na(mat_filt)) > 0
mat_filt  <- mat_filt[, keep_cols, drop = FALSE]
for (j in seq_len(ncol(mat_filt))) {
  v <- mat_filt[, j]
  if (anyNA(v)) v[is.na(v)] <- median(v, na.rm = TRUE)
  mat_filt[, j] <- v
}
keep_cols2 <- apply(mat_filt, 2, function(v) var(v, na.rm = TRUE) > 0)
mat_filt   <- mat_filt[, keep_cols2, drop = FALSE]

X     <- scale(mat_filt)
X     <- X[, colSums(is.nan(X)) == 0, drop = FALSE]
pca   <- prcomp(X, center = FALSE, scale. = FALSE)
n_pcs <- min(opt[["n-pcs"]], ncol(pca$x))
pcs   <- pca$x[, 1:n_pcs, drop = FALSE]
message("PCA: ", nrow(pcs), " cells x ", n_pcs, " PCs")

as_tibble(pcs, rownames = "cell_id") %>%
  mutate(well = str_extract(cell_id, "[A-P][0-9]{2}$"), .after = cell_id) %>%
  write_csv(opt[["pca-out"]])
message("Wrote: ", opt[["pca-out"]])
