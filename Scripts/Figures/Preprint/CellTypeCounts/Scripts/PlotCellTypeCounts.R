library(optparse)
library(tidyverse)


source(file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/figure_style.R"))
# CpG coverage of QC-passed single cells, split by cell type. FirstPrep P3+P4.
# Points coloured by plate of origin (displayed as P1/P2), n annotated per box.
#
# Sits next to the CovBoxplots QC-passed panel, which shows the same measurement
# split by cell count instead. Both derive their y axis from the same rule in
# Shared/FirstPrepCoverage.R over the same table, so the two axes are identical
# by construction rather than by anyone keeping numbers in sync.
#
# The limits are applied with coord_cartesian, which restricts the view only:
# wells outside the range still contribute to the boxplot statistics rather than
# being dropped from them.

option_list <- list(
  make_option("--well-stats-dir", type = "character",
              default = file.path(Sys.getenv("SNMC_PIPE"), "Analysis_Combined/Results/WellStats")),
  make_option("--shared", type = "character",
              default = file.path(Sys.getenv("SNMC_ROOT"), "Scripts/Figures/Preprint/Shared/FirstPrepCoverage.R")),
  make_option("--output", type = "character",
              help = "Output PNG")
)
opt <- parse_args(OptionParser(option_list = option_list))
source(opt$shared)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

all_wells <- read_firstprep_well_stats(opt[["well-stats-dir"]])
# Range over every QC-passed well with a cell in it -- including the 2/5/10-cell
# groups CovBoxplots shows -- so both panels end up on the same axis.
Y_RANGE <- qc_passed_y_range(all_wells)
message("shared y range: ", format(Y_RANGE[1], big.mark = ","), " .. ",
        format(Y_RANGE[2], big.mark = ","))

well_stats <- all_wells %>%
  filter(CellCount == 1, passes_cpg_qc == TRUE, CellType %in% names(CELLTYPE_DISPLAY)) %>%
  mutate(CellType = unname(CELLTYPE_DISPLAY[CellType]))

message("QC-passed single cells: ", nrow(well_stats))
print(well_stats %>% count(CellType, plate))
message("CpGsCovered range: ", format(min(well_stats$CpGsCovered), big.mark = ","),
        " .. ", format(max(well_stats$CpGsCovered), big.mark = ","))

n_labels <- well_stats %>%
  group_by(CellType) %>%
  summarise(n = n(), y_lab = n_label_height(max(CpGsCovered, na.rm = TRUE)), .groups = "drop")

set.seed(42)  # geom_jitter is random; fix it so the figure is reproducible
p <- ggplot(well_stats, aes(x = CellType, y = CpGsCovered)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = plate), width = 0.15, alpha = 0.6, size = 1.3) +
  geom_text(data = n_labels, aes(x = CellType, y = y_lab, label = n),
            inherit.aes = FALSE, size = 3.5) +
  # y expansion zero so coord_cartesian below gives exactly the shared range;
  # x padding keeps the outer boxes off the panel frame.
  scale_y_log10(labels = fmt_count, expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_color_manual(values = PLATE_COLORS, labels = PLATE_DISPLAY_LABELS, name = "Plate") +
  theme_bw(base_size = 13, base_family = FIGURE_FONT) +
  labs(x = "Cell Type", y = "CpGs (Log Scale)",
       title = "CpGs per Cell Type (QC-Passed)")

p <- p + coord_cartesian(ylim = Y_RANGE)
message("wells outside the view: ",
        sum(well_stats$CpGsCovered < Y_RANGE[1] | well_stats$CpGsCovered > Y_RANGE[2]))

save_figure( p, opt$output, width = 4.5, height = 5, dpi = 150)
message("Wrote: ", opt$output)
