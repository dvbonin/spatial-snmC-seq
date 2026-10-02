#!/bin/bash
# Redraws both heatmaps from the saved correlation tables. No computation -- use this
# for title/colour/layout changes instead of rerunning the array.
# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"
S=${SNMC_ROOT}/Scripts/Figures/Preprint/CorrHeatmaps/Scripts
W=${SNMC_FIGS}/Preprint/CorrHeatmaps
conda activate ${SNMC_ENV_R}
Rscript "$S/PlotCorrHeatmap.R" --input "$W/Data/SiteCorrelationPairs_AllSites.csv" \
  --title "Pairwise Methylation Correlation (All Shared CpG Sites)" \
  --output "$W/Output/CorrHeatmap_Sites_AllSites.png"
Rscript "$S/PlotCorrHeatmap.R" --input "$W/Data/SiteCorrelationPairs_Top1pct.csv" \
  --title "Pairwise Methylation Correlation (Top 1% Differential CpG Sites)" \
  --output "$W/Output/CorrHeatmap_Sites_Top1pct.png"

# Site counts behind each correlation
Rscript "$S/PlotCorrHeatmap.R" --input "$W/Data/SiteCorrelationPairs_AllSites.csv" --value n \
  --title "CpG Sites Shared by Each Pair (All Shared CpG Sites)" \
  --output "$W/Output/SiteCount_Sites_AllSites.png"
Rscript "$S/PlotCorrHeatmap.R" --input "$W/Data/SiteCorrelationPairs_Top1pct.csv" --value n \
  --title "CpG Sites Shared by Each Pair (Top 1% Differential CpG Sites)" \
  --output "$W/Output/SiteCount_Sites_Top1pct.png"
