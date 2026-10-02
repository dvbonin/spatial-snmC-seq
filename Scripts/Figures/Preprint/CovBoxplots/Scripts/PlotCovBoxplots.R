library(optparse)
library(tidyverse)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# CpG coverage against the number of cells dispensed into the well. FirstPrep
# P3+P4. Two panels:
#
#   QCPassed  -- QC-passed wells containing at least one cell. Shares its y axis
#                with the CellTypeCounts figure so the two can be placed side by
#                side; the range comes from Shared/FirstPrepCoverage.R, computed
#                from this same table.
#   AllWells  -- every well, unfiltered, empties included. Supplementary, so its
#                axis is left free to fill the panel.

option_list <- list(
  make_option("--well-stats-dir", type = "character",
              default = file.path(Sys.getenv("SNMC_PIPE"), "Analysis_Combined/Results/WellStats")),
  make_option("--shared", type = "character",
              default = file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/FirstPrepCoverage.R")),
  make_option("--output-qc-passed", type = "character"),
  make_option("--output-all-wells", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
source(opt$shared)
dir.create(dirname(opt[["output-qc-passed"]]), recursive = TRUE, showWarnings = FALSE)

well_stats <- read_firstprep_well_stats(opt[["well-stats-dir"]])
Y_RANGE <- qc_passed_y_range(well_stats)
message("shared y range: ", format(Y_RANGE[1], big.mark = ","), " .. ",
        format(Y_RANGE[2], big.mark = ","))

make_plot <- function(df, title, y_range = NULL) {
  n_labels <- df %>%
    group_by(CellCount) %>%
    summarise(n = n(), y_lab = n_label_height(max(CpGsCovered, na.rm = TRUE)), .groups = "drop")

  set.seed(42)  # geom_jitter is random; fix it so the figure is reproducible
  p <- ggplot(df, aes(x = factor(CellCount), y = CpGsCovered)) +
    geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    geom_jitter(aes(color = plate), width = 0.15, alpha = 0.6, size = 1.3, na.rm = TRUE) +
    geom_text(data = n_labels, aes(x = factor(CellCount), y = y_lab, label = n),
              inherit.aes = FALSE, size = 3.5) +
    # y: zero expansion when limits are supplied, so coord_cartesian's window is
    # exactly the shared range; the free panel keeps a little room for its labels.
    scale_y_log10(labels = fmt_count,
                  expand = if (is.null(y_range)) expansion(mult = c(0.03, 0.08))
                           else expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    scale_color_manual(values = PLATE_COLORS, labels = PLATE_DISPLAY_LABELS, name = "Plate") +
    theme_bw(base_size = 13, base_family = FIGURE_FONT) +
    labs(x = "Cells in Well", y = "CpGs (Log Scale)", title = title)

  if (!is.null(y_range)) {
    # expand stays TRUE here: the y scale above already contributes no expansion,
    # and turning it off would also strip the x padding that keeps the outer
    # boxes clear of the panel frame.
    p <- p + coord_cartesian(ylim = y_range)
    outside <- sum(df$CpGsCovered < y_range[1] | df$CpGsCovered > y_range[2], na.rm = TRUE)
    message("  wells outside the view: ", outside)
  }
  p
}

qc_passed <- well_stats %>% filter(CellCount > 0, passes_cpg_qc == TRUE)
message("QC-passed wells (CellCount > 0): ", nrow(qc_passed))
print(qc_passed %>% count(CellCount, plate))
save_figure(
       make_plot(qc_passed, "CpGs per Cell Count (QC-Passed Wells)", Y_RANGE), opt[["output-qc-passed"]],
       width = 7, height = 5, dpi = 150)
message("Wrote: ", opt[["output-qc-passed"]])

message("All wells: ", nrow(well_stats))
save_figure(
       make_plot(well_stats, "CpGs per Cell Count (All Wells)"), opt[["output-all-wells"]],
       width = 7, height = 5, dpi = 150)
message("Wrote: ", opt[["output-all-wells"]])
