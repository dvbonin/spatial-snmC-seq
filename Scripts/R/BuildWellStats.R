library(optparse)
library(tidyverse)
library(data.table)

# Consolidated per-well stats table. Methylation numbers are NOT recomputed here:
# they come straight from the PP pipeline's MethStats_Summary.csv, which
# WellMethStats.py already derived from the same per-well VCFs. The Analysis
# pipeline used to rescan every VCF a second time to produce three columns that
# were byte-identical to that table and a fourth that differed only by estimator.
#
# MeanGlobMeth_* is the site-weighted mean: each covered cytosine's own M/(M+U)
# averaged over sites, so every site counts once regardless of coverage. All
# three contexts use the same estimator (the old table mixed a site mean for CpG
# with read-pooled ratios for CHG/CHH).
#
# passes_cpg_qc is deliberately empty. QC is no longer decided in this pipeline;
# it is a downstream, figure-level decision, and the column is kept only so
# consumers do not have to be restructured before that decision is wired in.

option_list <- list(
  make_option("--mappings",      type = "character",
              help = "Comma-separated list of Mapping_*.csv files"),
  make_option("--main-counts",   type = "character",
              help = "Main_read_counts.csv (WellPosition;ReadCounts)"),
  make_option("--alt-counts",    type = "character", default = NULL,
              help = "Alt_read_counts.csv (non-primary contigs: chrM/EBV/decoy/scaffolds)"),
  make_option("--lambda-counts", type = "character", default = NULL,
              help = "Lambda_read_counts.csv [optional, spike-ins only]"),
  make_option("--puc19-counts",  type = "character", default = NULL,
              help = "pUC19_read_counts.csv [optional, spike-ins only]"),
  make_option("--meth-summary",  type = "character",
              help = "PP pipeline MethStats_Summary.csv (well, {CG,CHG,CHH}_{n_sites,n_calls,n_methylated,mean_beta})"),
  make_option("--output",        type = "character",
              help = "Output consolidated per-well stats CSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

mapping_paths <- strsplit(opt$mappings, ",")[[1]]
output_path   <- opt$output
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

well_meta <- map_dfr(mapping_paths, ~ fread(.x, sep = ";", data.table = FALSE)) %>%
  select(WellPosition, CellType, CellCount) %>%
  distinct() %>%
  mutate(well_class = case_when(
    CellCount == 0 ~ "Empty",
    CellCount == 1 ~ CellType,
    TRUE           ~ NA_character_
  ))

read_counts <- function(path, col_name) {
  df <- fread(path, sep = ";", data.table = FALSE) %>% select(WellPosition, ReadCounts)
  names(df)[names(df) == "ReadCounts"] <- col_name
  df
}

stats <- well_meta %>%
  left_join(read_counts(opt[["main-counts"]], "ReadCounts_Main"), by = "WellPosition")

for (spec in list(c("alt-counts", "ReadCounts_Alt"),
                  c("lambda-counts", "ReadCounts_Lambda"),
                  c("puc19-counts", "ReadCounts_pUC19"))) {
  if (!is.null(opt[[spec[1]]])) {
    stats <- stats %>% left_join(read_counts(opt[[spec[1]]], spec[2]), by = "WellPosition")
  }
}
# A stream only lists wells it actually observed, so a well missing from Alt /
# Lambda / pUC19 genuinely had zero such reads rather than being unmeasured.
count_cols <- intersect(c("ReadCounts_Alt", "ReadCounts_Lambda", "ReadCounts_pUC19"), names(stats))
stats <- stats %>% mutate(across(all_of(count_cols), ~ replace_na(.x, 0L)))

meth <- fread(opt[["meth-summary"]], data.table = FALSE) %>%
  transmute(WellPosition     = well,
            CpGsCovered      = CG_n_sites,
            MeanGlobMeth_CpG = CG_mean_beta,
            MeanGlobMeth_CHG = CHG_mean_beta,
            MeanGlobMeth_CHH = CHH_mean_beta)

stats <- stats %>%
  left_join(meth, by = "WellPosition") %>%
  mutate(passes_cpg_qc = NA)

write_csv(stats, output_path)
message("Wrote: ", output_path, " (", nrow(stats), " wells, ", ncol(stats), " columns)")
