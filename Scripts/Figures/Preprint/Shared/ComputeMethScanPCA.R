library(optparse)
library(tidyverse)
library(data.table)

# PCA on a MethSCAn mean-shrunken-residuals matrix, over ALL wells in it -- Empty,
# single-cell and multi-cell alike, any cell type. QC filtering, where it happens at
# all, happens downstream of this.
#
# Shared by two figures, on two different matrices:
#   ../QC     the merged FirstPrep P3+P4 matrix from ../DataPrePrep
#   (a per-plate matrix from the Analysis pipeline works the same way)
#
# cell_id is whatever the matrix calls the row; well is the trailing A01-style position,
# which is the same in both layouts.

option_list <- list(
  make_option("--matrix", type = "character", help = "mean_shrunken_residuals.csv.gz"),
  make_option("--min-vmrs", type = "integer", default = 1L),
  make_option("--pca-out", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt[["pca-out"]]), recursive = TRUE, showWarnings = FALSE)

message("Loading residuals matrix: ", opt$matrix)
mat_df <- fread(opt$matrix, data.table = FALSE)
rownames(mat_df) <- mat_df[[1]]
mat_df[[1]] <- NULL
mat_df <- as.matrix(mat_df)

annotation <- tibble(cell_id = rownames(mat_df)) %>%
  mutate(n_vmrs = rowSums(!is.na(mat_df))) %>%
  filter(n_vmrs >= opt[["min-vmrs"]])
message("Keeping ", nrow(annotation), " wells (Empty / single-cell / multi-cell, any cell type)")

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

X     <- scale(mat_filt)
X     <- X[, colSums(is.nan(X)) == 0, drop = FALSE]
pca   <- prcomp(X, center = FALSE, scale. = FALSE)
n_pcs <- min(30, ncol(pca$x))
pcs   <- pca$x[, 1:n_pcs, drop = FALSE]
message("PCA: ", nrow(pcs), " cells x ", n_pcs, " PCs")

pca_df <- as_tibble(pcs, rownames = "cell_id") %>%
  mutate(well = str_extract(cell_id, "[A-P][0-9]{2}$"), .after = cell_id)
write_csv(pca_df, opt[["pca-out"]])
message("Wrote: ", opt[["pca-out"]])
