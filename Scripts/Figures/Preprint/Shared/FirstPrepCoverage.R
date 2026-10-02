# Shared by the CovBoxplots and CellTypeCounts figures, which are meant to be
# placed side by side and therefore must use an identical y axis.
#
# The axis is derived from the data rather than hardcoded: both figures read the
# same well_stats tables, so each computing the range itself gives the same
# answer and they cannot drift apart. An earlier version had the range pasted as
# literal numbers into one figure after being read off the other's rendered PNG;
# it went stale the moment the QC definition changed.

FIRSTPREP_PLATES <- tibble::tribble(~sample, ~plate, "FirstPrep", "P3", "FirstPrep", "P4")

read_firstprep_well_stats <- function(well_stats_dir) {
  purrr::pmap_dfr(FIRSTPREP_PLATES, function(sample, plate) {
    readr::read_csv(file.path(well_stats_dir, paste0(sample, "_", plate, "_well_stats.csv")),
                    show_col_types = FALSE) %>%
      dplyr::mutate(plate = plate)
  })
}

# Range spanned by every QC-passed well that contains at least one cell, padded
# outwards. Deliberately not restricted to single cells: CovBoxplots shows the
# 2/5/10-cell groups too, and the shared axis has to hold all of them.
#
# The padding is a fraction of the LOG span, not a multiplier on the value. A
# plain `max * 1.2` is only log10(1.2) = 0.08 of a ~2.4-decade axis, i.e. about
# 3% of the panel height, which leaves the tallest points looking pinned to the
# top edge. The ceiling is additionally floored at `max * 1.8` so an n label at
# `group max * 1.5` always has room above the highest point.
qc_passed_y_range <- function(well_stats, pad_lower = 0.04, pad_upper = 0.14,
                              n_label_mult = 1.5) {
  cov <- well_stats %>%
    dplyr::filter(CellCount > 0, passes_cpg_qc == TRUE) %>%
    dplyr::pull(CpGsCovered)
  mn <- min(cov, na.rm = TRUE); mx <- max(cov, na.rm = TRUE)
  span <- log10(mx) - log10(mn)
  c(10^(log10(mn) - pad_lower * span),
    max(10^(log10(mx) + pad_upper * span), mx * n_label_mult * 1.2))
}

# n label for a group sits a fixed factor above that group's highest point, so
# every panel places them the same way. qc_passed_y_range() reserves enough
# headroom above the global maximum for this to stay inside the view.
n_label_height <- function(group_max, mult = 1.5) group_max * mult

fmt_count <- function(x) {
  ifelse(x >= 1e9, paste0(signif(x / 1e9, 3), "B"),
  ifelse(x >= 1e6, paste0(signif(x / 1e6, 3), "M"),
  ifelse(x >= 1e3, paste0(signif(x / 1e3, 3), "k"),
                   as.character(signif(x, 3)))))
}

PLATE_COLORS         <- c(P3 = "#8B7CC0", P4 = "#D9A441")
PLATE_DISPLAY_LABELS <- c(P3 = "P1", P4 = "P2")   # display-only; data stays P3/P4
CELLTYPE_DISPLAY     <- c(HealthyKeratinocytes = "Keratinocytes",
                          MerkelCarcinoma = "Merkel Cell Carcinoma")
