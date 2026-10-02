library(optparse)
library(tidyverse)
library(data.table)

# Fills passes_cpg_qc in the FirstPrep well_stats tables from this directory's QC
# metric: a well passes if its mean fraction-empty-in-own-cluster, averaged over
# the Leiden sweep, is below --threshold.
#
# The verdict is on the WELL's data quality alone and deliberately ignores cell
# count, so the column composes rather than pre-empting: downstream code asks for
# `CellCount == 1 & passes_cpg_qc` to get QC-passed single cells, which is how
# every consumer already writes it. Marking the column purely on the metric means
# empty and multi-cell wells also get a verdict, which is what makes the two
# conditions separable.
#
# The same table is also written out on its own (--qc-pass-out). well_stats.csv is
# an Analysis.smk output, so a rerun of build_well_stats will reset this column to
# NA; that standalone copy is the durable record of the decision, and re-running
# this script restores the column.

option_list <- list(
  make_option("--well-empty-frac-csv", type = "character",
              help = "FirstPrep_Merged_well_empty_frac.csv from SweepLeidenClusterUMAP.R"),
  make_option("--well-stats-p3", type = "character"),
  make_option("--well-stats-p4", type = "character"),
  make_option("--threshold", type = "double", default = 0.1,
              help = "Wells with mean_empty_frac below this pass QC [default: %default]"),
  make_option("--qc-pass-out", type = "character",
              help = "Standalone copy of the QC verdict per well")
)
opt <- parse_args(OptionParser(option_list = option_list))

ef <- fread(opt[["well-empty-frac-csv"]], data.table = FALSE) %>%
  mutate(passes_cpg_qc = mean_empty_frac < opt$threshold)

apply_one <- function(path, plate) {
  stats <- fread(path, data.table = FALSE)
  verdict <- ef %>% filter(Plate == plate) %>% select(well, passes_cpg_qc)
  stats$passes_cpg_qc <- NULL
  stats <- stats %>% left_join(verdict, by = c("WellPosition" = "well"))
  write_csv(stats, path)
  message(plate, ": ", sum(stats$passes_cpg_qc, na.rm = TRUE), " / ", nrow(stats),
          " wells pass (mean_empty_frac < ", opt$threshold, "); ",
          sum(is.na(stats$passes_cpg_qc)), " unscored")
  stats %>% transmute(Plate = plate, WellPosition, CellType, CellCount, passes_cpg_qc)
}

out <- bind_rows(apply_one(opt[["well-stats-p3"]], "P3"),
                 apply_one(opt[["well-stats-p4"]], "P4")) %>%
  left_join(ef %>% select(Plate, well, mean_empty_frac),
            by = c("Plate", "WellPosition" = "well")) %>%
  relocate(mean_empty_frac, .before = passes_cpg_qc)

dir.create(dirname(opt[["qc-pass-out"]]), recursive = TRUE, showWarnings = FALSE)
write_csv(out, opt[["qc-pass-out"]])
message("Wrote: ", opt[["qc-pass-out"]])

single <- out %>% filter(CellCount == 1)
message("QC-passed single cells (CellCount == 1 & passes_cpg_qc): ",
        sum(single$passes_cpg_qc, na.rm = TRUE), " / ", nrow(single))
