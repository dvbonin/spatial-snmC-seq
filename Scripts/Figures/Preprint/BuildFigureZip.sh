#!/bin/bash
# Collects the figure panels for the preprint into one flat zip for download and
# layout work: editable SVG for every panel, plus the PubSeq TIFF (that panel's
# vector form runs to ~240 MB, so it ships as a print-resolution raster instead).
#
# The SVGs only exist if the pipelines were run with SNMC_SVG set -- see
# Shared/figure_style.{R,py}. A missing panel aborts the script rather than
# producing a quietly incomplete archive.
#
# Usage: bash BuildFigureZip.sh [svg|pdf] [output.zip]
#
# PDF is the more reliable import into Illustrator (embedded TrueType Arial
# subsets); SVG keeps text as markup. Both are editable text, neither outlines.

set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

B=${SNMC_FIGS}/Preprint
FMT=${1:-svg}
OUT=${2:-$B/PreprintFigures_${FMT}.zip}
TMP=${SNMC_TMP}

PANELS=(
  CellTypeCounts/Output/CellTypeCounts_FirstPrep
  CorrHeatmaps/Output/CorrHeatmap_Sites_AllSites
  CorrHeatmaps/Output/CorrHeatmap_Sites_Top1pct
  CovBoxplots/Output/CovBoxplots_FirstPrep_AllWells
  CovBoxplots/Output/CovBoxplots_FirstPrep_QCPassed
  SpikeIn/Output/FirstPrep_Lambda_coverage_box
  SpikeIn/Output/FirstPrep_pUC19_coverage_box
  SpikeIn/Output/FirstPrep_spikein_scatter
  QC/Output/FirstPrep_Merged_coverage_vs_empty_frac
  QC/Output/FirstPrep_Merged_panel6
  Microarray/Output/FirstPrep_Merged_umap_fused
  ExpressionData/Output/genebody_metagene_by_quartile
  ExpressionData/Output/promoter_metagene_by_quartile_5kb
  ExpressionData/Output/expression_percentiles
  MarkerGenes/Output/marker_gene_violins_promoter
  MarkerGenes/Output/marker_gene_violins_body
  ReadCountHistogram/Output/ReadCountHistogram_FirstPrep
  UsableReads/Output/usable_reads_firstprep
)

STAGE=$(mktemp -d "$TMP/figzip.XXXXXX")
trap 'rm -rf "$STAGE"' EXIT

for p in "${PANELS[@]}"; do cp "$B/$p.$FMT" "$STAGE/"; done
# PubSeq ships as a raster in either case: ~560k points make its vector form unusable.
cp "$B/PubSeq/Output/SingleCell_vs_Public_10kbBins.tiff" "$STAGE/"

rm -f "$OUT"
(cd "$STAGE" && zip -q -9 "$OUT" *)

echo "Wrote $OUT"
unzip -l "$OUT" | tail -3
