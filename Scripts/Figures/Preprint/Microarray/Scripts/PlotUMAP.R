library(optparse)
library(tidyverse)
library(umap)
library(data.table)
library(patchwork)
library(grid)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# UMAP of the QC-passed FirstPrep single cells (P3+P4, merged MethSCAn run), in
# two panels sharing one embedding:
#
#   Cell Type            the dispensed annotation
#   Corr. to Microarray  each cell's Pearson r to the Merkel-carcinoma EPIC array
#                        minus its r to the keratinocyte array, over the top-diff
#                        sites. Positive = looks like the carcinoma reference.
#
# The second panel is the independent check on the first: cell identity read off
# the methylation itself, against nothing but bulk array references, so the two
# panels agreeing is evidence the clusters are real cell types and not a batch or
# coverage artefact.
#
# The PCA is not recomputed here -- it is the one the QC directory used to make the
# pass/fail call (../../QC/Data/FirstPrep_Merged_qcpassed_pca.csv.gz), so the cells
# plotted and the cells judged are the same set by construction.

option_list <- list(
  make_option("--pca", type = "character",
              help = "FirstPrep_Merged_qcpassed_pca.csv.gz from ../../QC"),
  make_option("--well-stats-dir", type = "character",
              default = file.path(Sys.getenv("SNMC_PIPE"), "Analysis_Combined/Results/WellStats")),
  make_option("--correlations-p3", type = "character"),
  make_option("--correlations-p4", type = "character"),
  make_option("--shared", type = "character",
              default = file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/FirstPrepCoverage.R")),
  make_option("--n-neighbors", type = "integer", default = 25),
  make_option("--min-dist", type = "double", default = 0.7),
  make_option("--topdiff-label", type = "character", default = "Top 1% Diff."),
  make_option("--output", type = "character"),
  make_option("--output-fused", type = "character",
              help = "Single-panel variant with the cell-type UMAP inset bottom-right"),
  make_option("--fused-title", type = "character",
              default = "Methylation Similarity to Merkel vs. Keratinocyte Arrays")
)
opt <- parse_args(OptionParser(option_list = option_list))
source(opt$shared)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

pcs <- fread(opt$pca, data.table = FALSE)
pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
message("PCA: ", nrow(pcs), " cells x ", length(pc_cols), " PCs")

well_stats <- read_firstprep_well_stats(opt[["well-stats-dir"]]) %>%
  select(Plate = plate, WellPosition, CellType, CellCount, passes_cpg_qc)

correlations <- bind_rows(
  fread(opt[["correlations-p3"]], data.table = FALSE) %>% mutate(Plate = "P3"),
  fread(opt[["correlations-p4"]], data.table = FALSE) %>% mutate(Plate = "P4")
)

ann <- pcs %>%
  select(cell_id, Plate, well) %>%
  left_join(well_stats, by = c("Plate", "well" = "WellPosition")) %>%
  left_join(correlations, by = c("Plate", "well"))

stopifnot(all(ann$CellCount == 1 & ann$passes_cpg_qc))
message("Cells: ", nrow(ann), " (all QC-passed single cells)")
message("With a microarray correlation: ", sum(!is.na(ann$corr_to_mcc)))

set.seed(42)
um <- umap(as.matrix(pcs[, pc_cols]),
           n_neighbors = min(opt[["n-neighbors"]], nrow(pcs) - 1),
           min_dist = opt[["min-dist"]], metric = "euclidean")
ann$UMAP1 <- um$layout[, 1]
ann$UMAP2 <- um$layout[, 2]

umap_theme <- function() {
  theme_minimal(base_size = 13, base_family = FIGURE_FONT) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text = element_blank(), axis.ticks = element_blank())
}

# Colours are assigned over the alphabetically sorted types, so each type keeps the
# same colour no matter what order the legend lists them in.
cell_types <- sort(setdiff(unique(ann$CellType), NA))
CELLTYPE_COLORS <- setNames(
  c("#4393c3", "#d95f02", "#7570b3", "#66a61e", "#e6ab02", "#a6761d")[seq_along(cell_types)],
  cell_types)
# Named, so ggplot matches label to type rather than to position -- which is what
# lets `breaks` below reorder the legend safely.
ct_labels <- setNames(ifelse(cell_types %in% names(CELLTYPE_DISPLAY),
                             CELLTYPE_DISPLAY[cell_types], cell_types),
                      cell_types)
# Legend order only: Merkel carcinoma first. Any type not named here follows.
CT_ORDER  <- intersect(c("MerkelCarcinoma", "HealthyKeratinocytes"), cell_types)
CT_ORDER  <- c(CT_ORDER, setdiff(cell_types, CT_ORDER))

p_ct <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_manual(values = CELLTYPE_COLORS, labels = ct_labels,
                      breaks = CT_ORDER, na.value = "grey70", name = "Cell Type") +
  umap_theme() +
  labs(title = "Cell Type")

# Difference of the two correlations rather than either alone: it cancels the
# per-cell coverage effect that pushes both r's down together in sparse cells.
ann$corr_diff <- ann$corr_to_mcc - ann$corr_to_kz

p_diff <- ggplot(ann, aes(UMAP1, UMAP2, colour = corr_diff)) +
  geom_point(size = 1.8, alpha = 0.5) +
  scale_colour_gradient2(low = "#0868AC", mid = "grey90", high = "#D94801",
                         midpoint = 0, na.value = "grey85", name = "Δr") +
  umap_theme() +
  labs(title = paste0("Corr. to Microarray (", opt[["topdiff-label"]], ")"))

save_figure( p_ct | p_diff, opt$output, width = 14, height = 6, dpi = 150)
message("Wrote: ", opt$output)

# ---- Fused variant: one correlation panel, cell types inset bottom-right ----
#
# The same embedding as above, drawn once. The axis range is widened so the
# bottom-right corner is empty enough to hold the inset without covering cells.
#
# The inset is inset from the panel edge by the same PHYSICAL distance on both
# axes, which is why the offset is given in cm rather than as a fraction: the
# panel is wider than it is tall, so equal fractions would read as unequal gaps.
#
# The colour ramp runs blue -> light -> red across the observed range, so red
# reads as "most correlated with the carcinoma reference". Deliberately not
# centred on zero: this variant is meant to be read as a sequential scale.
#
# Legends are extracted as grobs and placed in their own column with explicit
# spacer heights, rather than collected by patchwork. That is the only way to put
# each one at the mid-height of the plot it describes -- Delta-r against the main
# panel's centre, Cell Type against the inset's.
if (!is.null(opt[["output-fused"]])) {
  # ---- tunable geometry -------------------------------------------------
  # PAD_* are fractions of the data range added outside it, so smaller values
  # zoom in on the clusters. The right/top values are larger than the
  # left/bottom ones only to open up the corner the inset sits in.
  PAD_X      <- c(left = 0.05, right = 0.18)
  PAD_Y      <- c(bottom = 0.05, top = 0.18)
  SHIFT_X    <- -0.10                  # pan the window, as a fraction of its width;
                                       # negative pans left, so the cells sit further right
  CROP_X_LEFT <- 0.10                  # drop this fraction off the left of the window
  CROP_Y_TOP <- 0.10                   # and this off the top
  INSET_FRAC <- 0.32                   # inset size, as a fraction of the panel
  GAP_CM     <- 0.32                   # gap from the panel edge, same both axes
  INSET_BORDER <- "grey55"             # inset frame colour
  FIG_W      <- 10.0; FIG_H <- 7       # inches; wider than the panel needs so
                                     # the legend column can hold the offsets below
  LEGEND_W   <- 0.27                   # legend column width, relative to panel
  # Shared horizontal centre for both legends, measured into the legend column.
  # Because the two blocks differ in width (Delta-r ~2 cm, Cell Type ~4.5 cm),
  # aligning their centres necessarily insets the narrow one: the floor is half
  # the wider block, so anything smaller would push Cell Type off the column.
  LEG_CENTER_CM <- 2.35
  legend_col_cm  <- FIG_W * 2.54 * LEGEND_W / (1 + LEGEND_W)
  PT_MAIN    <- 2.1;  A_MAIN  <- 0.55  # point size / alpha, main panel
  PT_INSET   <- 1.0;  A_INSET <- 0.6   # point size / alpha, inset
  INSET_LIGHTEN <- 0.42                # blend the cell-type colours this far toward
                                       # white, so they do not outweigh the pale
                                       # Delta-r points in the main panel
  RAMP       <- c("#2166AC", "#F7F7F7", "#B2182B")   # low -> mid -> high
  # -----------------------------------------------------------------------

  # Lightened only for the inset; the two-panel figure keeps the full-strength
  # palette. The legend is extracted from this plot, so its keys lighten to match.
  lighten <- function(cols, f) {
    m <- col2rgb(cols)
    setNames(rgb(t(m + (255 - m) * f), maxColorValue = 255), names(cols))
  }

  pad <- function(v, lo, hi) {
    r <- range(v, na.rm = TRUE); s <- diff(r)
    c(r[1] - s * lo, r[2] + s * hi)
  }
  xlim <- pad(ann$UMAP1, PAD_X[["left"]],   PAD_X[["right"]])
  xlim <- xlim + SHIFT_X * diff(xlim)
  xlim[1] <- xlim[1] + CROP_X_LEFT * diff(xlim)
  ylim <- pad(ann$UMAP2, PAD_Y[["bottom"]], PAD_Y[["top"]])
  ylim[2] <- ylim[2] - CROP_Y_TOP * diff(ylim)
  GAP  <- unit(GAP_CM, "cm")

  message("fused geometry: inset ", INSET_FRAC, " of panel, gap ", GAP_CM, " cm, ",
          "pad x ", PAD_X[["left"]], "/", PAD_X[["right"]],
          ", pad y ", PAD_Y[["bottom"]], "/", PAD_Y[["top"]],
          ", shift x ", SHIFT_X, ", crop left ", CROP_X_LEFT, "/top ", CROP_Y_TOP, ", legend w ", LEGEND_W,
          ", inset lighten ", INSET_LIGHTEN,
          "\n  legend col ", round(legend_col_cm, 2), " cm, shared centre ",
          LEG_CENTER_CM, " cm")

  # No title on the inset: the Cell Type legend outside already names it.
  p_inset <- ggplot(ann, aes(UMAP1, UMAP2, colour = CellType)) +
    geom_point(size = PT_INSET, alpha = A_INSET) +
    scale_colour_manual(values = lighten(CELLTYPE_COLORS, INSET_LIGHTEN),
                        labels = ct_labels, breaks = CT_ORDER,
                        na.value = "grey80", name = "Cell Type") +
    # The legend keys need their own size and full opacity: at the inset's point
    # size these lightened dots would be almost invisible in the legend.
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    theme_minimal(base_size = 11, base_family = FIGURE_FONT) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = INSET_BORDER, fill = NA, linewidth = 0.7),
          plot.background = element_rect(fill = "white", colour = NA),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0))

  p_fused <- ggplot(ann, aes(UMAP1, UMAP2, colour = corr_diff)) +
    geom_point(size = PT_MAIN, alpha = A_MAIN) +
    scale_colour_gradientn(colours = RAMP, name = "\u0394r") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    umap_theme() +
    labs(title = opt[["fused-title"]])

  main <- p_fused + theme(legend.position = "none") +
    inset_element(p_inset + theme(legend.position = "none"),
                  left   = unit(1 - INSET_FRAC, "npc") - GAP,
                  bottom = unit(1 - INSET_FRAC, "npc") - GAP,
                  right  = unit(1, "npc") - GAP,
                  top    = unit(1, "npc") - GAP,
                  align_to = "panel")

  # suppressWarnings is narrow and deliberate: converting to a grob probes font
  # metrics through a PostScript device that has no Nimbus Sans, emitting one
  # warning per text element (~85 per run). The PNG device resolves the font
  # normally, so the output is unaffected.
  legend_of <- function(p) {
    g <- suppressWarnings(ggplotGrob(p + theme(legend.position = "right")))
    g$grobs[[which(g$layout$name == "guide-box-right")]]
  }

  # wrap_elements centres a grob, which ties horizontal position to the grob's own
  # width. Padding it out to the full column width instead makes the centring a
  # no-op, so the left pad alone fixes where the content sits.
  # Pads on the left so the block's own centre lands on center_cm, then fills to
  # the column width so wrap_elements' centring becomes a no-op.
  place_h <- function(g, center_cm) {
    w0 <- convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
    left_cm <- center_cm - w0 / 2
    if (left_cm < 0) {
      warning("legend is ", round(w0, 2), " cm wide; centring it at ", center_cm,
              " cm would start it off-column -- raise LEG_CENTER_CM to at least ",
              round(w0 / 2, 2), call. = FALSE, immediate. = TRUE)
      left_cm <- 0
    }
    g <- gtable::gtable_add_cols(g, unit(left_cm, "cm"), pos = 0)
    w <- convertWidth(sum(g$widths), "cm", valueOnly = TRUE)
    if (w > legend_col_cm)
      warning("legend (", round(w0, 2), " cm + ", left_cm, " cm pad = ", round(w, 2),
              " cm) exceeds the ", round(legend_col_cm, 2),
              " cm column; it cannot shift that far -- raise FIG_W or LEGEND_W",
              call. = FALSE, immediate. = TRUE)
    message("    legend ", round(w0, 2), " cm wide -> left pad ", round(left_cm, 2),
            " cm, centre ", round(left_cm + w0 / 2, 2),
            " cm  (column ", round(legend_col_cm, 2), " cm)")
    if (legend_col_cm > w) g <- gtable::gtable_add_cols(g, unit(legend_col_cm - w, "cm"), pos = -1)
    g
  }

  # Vertical targets as a fraction from the top of the figure. The panel spans
  # almost the full height (a title above, an axis label below, roughly equal),
  # so its centre is ~0.50; the inset's centre sits INSET_FRAC/2 above the panel
  # floor. Both are approximations of the rendered geometry -- close enough that
  # each legend reads as belonging to its plot.
  panel_lo <- 0.069; panel_hi <- 0.945         # panel top/bottom, fig fractions
  panel_h  <- panel_hi - panel_lo
  gap_frac <- GAP_CM / (FIG_H * 2.54)
  y_corr   <- panel_lo + panel_h * 0.50
  # The inset's centre is half its height below the panel ceiling, itself dropped
  # by GAP -- dropping that term puts this legend a visible notch too low.
  y_type   <- panel_lo + gap_frac + panel_h * (INSET_FRAC / 2)
  # A guide box carries padding below its content, so wrap_elements centres the
  # grob a little above where the visible block reads as centred. Measured off the
  # rendered figure; revisit if the legend contents change substantially.
  NUDGE    <- 0.014
  NUDGE_CT <- 0.011   # the discrete legend needs a little more than the colourbar

  h_dr <- 0.18; h_ct <- 0.11
  leg_dr <- place_h(legend_of(p_fused), LEG_CENTER_CM)
  leg_ct <- place_h(legend_of(p_inset), LEG_CENTER_CM)

  # Ordered by target rather than hardcoded: the Cell Type block belongs beside the
  # inset, so it leads when the inset is at the top and trails when it is at the
  # bottom. Targets are fractions from the top of the figure, so smaller is higher.
  blocks <- list(list(y = y_corr + NUDGE,             h = h_dr, g = leg_dr),
                 list(y = y_type + NUDGE + NUDGE_CT,  h = h_ct, g = leg_ct))
  blocks <- blocks[order(vapply(blocks, function(b) b$y, numeric(1)))]
  s1 <- blocks[[1]]$y - blocks[[1]]$h / 2
  s2 <- blocks[[2]]$y - blocks[[2]]$h / 2 - (s1 + blocks[[1]]$h)
  s3 <- 1 - (s1 + blocks[[1]]$h + s2 + blocks[[2]]$h)
  stopifnot(s1 > 0, s2 > 0, s3 > 0)

  legend_col <- plot_spacer() / wrap_elements(full = blocks[[1]]$g) /
    plot_spacer() / wrap_elements(full = blocks[[2]]$g) / plot_spacer() +
    plot_layout(heights = c(s1, blocks[[1]]$h, s2, blocks[[2]]$h, s3))

  # Closest approach of any cell to the inset region, in panel fractions -- if
  # either margin goes negative the inset is covering data.
  fx <- (ann$UMAP1 - xlim[1]) / diff(xlim)
  fy <- (ann$UMAP2 - ylim[1]) / diff(ylim)
  inside <- fx > (1 - INSET_FRAC) & fy > (1 - INSET_FRAC)
  message("  cells inside the inset region: ", sum(inside),
          "  (nearest approach: x ", round(1 - INSET_FRAC - max(fx[fy > (1 - INSET_FRAC)]), 3),
          ", y ", round(1 - INSET_FRAC - max(fy[fx > (1 - INSET_FRAC)]), 3), ")")

  save_figure( (main | legend_col) + plot_layout(widths = c(1, LEGEND_W)), opt[["output-fused"]],
         width = FIG_W, height = FIG_H, dpi = 150)
  message("Wrote: ", opt[["output-fused"]])
}
