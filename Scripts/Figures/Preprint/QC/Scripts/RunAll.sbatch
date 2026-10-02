#!/bin/bash
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o Logs/RunAll_%j.log
#SBATCH -e Logs/RunAll_%j.log
set -e
: "${SNMC_ROOT:?source the repository config.sh first}"

# R and python both put scratch files under TMPDIR; without this they fall back to
# /tmp on the compute node, which fills up and fails the job with errors like
# "cannot create 'R_TempDir'".
export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

# 1) PCA on all FirstPrep P3+P4 wells, from DataPrePrep's merged MethSCAn matrix.
# 2) Leiden clustering sweep (n_neighbors x resolution grid) -- the per-well
#    mean-empty-fraction CSV (plotted separately, in Python -- see below), and
#    the representative (25, 0.5) combo's own embedding, exported for panel6.
# 3) panel6 UMAP (Cell Type / Plate / Cell Count / Coverage / Mean Methylation
#    / QC), on that representative embedding.
# 4) The actual QC-passed clustering: an n_neighbors x min_dist UMAP sweep,
#    fresh PCA on just the threshold-passed single cells, plus the same sweep
#    on all single cells (passed or not, for comparison).
# THRESHOLD is a default, not a human decision each run -- override by editing
# it here or rerunning PlotMergedFirstPrepPanel6.R / SweepQCThresholdClustering.R
# by hand with --threshold.

THRESHOLD=0.1

ROOT=${SNMC_FIGS}/Preprint/QC
SCRIPTS=${SNMC_ROOT}/Scripts/Figures/Preprint/QC/Scripts
SHARED=${SNMC_ROOT}/Scripts/Figures/Preprint/Shared
DATA="$ROOT/Data"
OUTPUT="$ROOT/Output"
MERGED_MATRIX=${SNMC_FIGS}/Preprint/DataPrePrep/Data/MergedFirstPrep/Matrix/mean_shrunken_residuals.csv.gz
WELLSTATS_DIR=${SNMC_PIPE}/Analysis_Combined/Results/WellStats
PQ_PYTHON=${SNMC_ENV_PQ}/bin/python
mkdir -p "$DATA" "$OUTPUT"
source ~/.bashrc

conda activate ${SNMC_ENV_R}

echo "=== FirstPrep Merged (P3+P4): PCA on all wells ==="
Rscript "$SHARED/ComputeMethScanPCA.R" \
  --matrix "$MERGED_MATRIX" \
  --pca-out "$DATA/FirstPrep_Merged_methscan_pca.csv.gz"

echo "=== FirstPrep Merged (P3+P4): Leiden clustering sweep ==="
Rscript "$SCRIPTS/SweepLeidenClusterUMAP.R" \
  --pca "$DATA/FirstPrep_Merged_methscan_pca.csv.gz" \
  --well-stats-p3 "$WELLSTATS_DIR/FirstPrep_P3_well_stats.csv" \
  --well-stats-p4 "$WELLSTATS_DIR/FirstPrep_P4_well_stats.csv" \
  --output-umap-grid "$OUTPUT/FirstPrep_Merged_leiden_sweep_umap_grid.png" \
  --output-well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --output-representative-embedding-csv "$DATA/FirstPrep_Merged_representative_embedding.csv"

echo "=== Write the QC verdict into the well_stats tables ==="
Rscript "$SCRIPTS/ApplyQCToWellStats.R" \
  --well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --well-stats-p3 "$WELLSTATS_DIR/FirstPrep_P3_well_stats.csv" \
  --well-stats-p4 "$WELLSTATS_DIR/FirstPrep_P4_well_stats.csv" \
  --threshold "$THRESHOLD" \
  --qc-pass-out "$DATA/FirstPrep_Merged_qc_pass.csv"

echo "=== FirstPrep Merged (P3+P4): per-well empty-fraction beeswarm (Python -- see project memory on the broken R Makeconf) ==="
"$PQ_PYTHON" "$SCRIPTS/PlotEmptyFractionBeeswarm.py" \
  --well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --n-sweeps 9 \
  --output "$OUTPUT/FirstPrep_Merged_leiden_sweep_empty_frac_rank.png"

echo "=== FirstPrep Merged (P3+P4): coverage vs own-cluster empty fraction (Python) ==="
"$PQ_PYTHON" "$SCRIPTS/PlotCoverageVsEmptyFraction.py" \
  --well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --well-stats-p3 "$WELLSTATS_DIR/FirstPrep_P3_well_stats.csv" \
  --well-stats-p4 "$WELLSTATS_DIR/FirstPrep_P4_well_stats.csv" \
  --threshold "$THRESHOLD" \
  --n-sweeps 9 \
  --output "$OUTPUT/FirstPrep_Merged_coverage_vs_empty_frac.png"

echo "=== FirstPrep Merged (P3+P4): panel6 UMAP ==="
Rscript "$SCRIPTS/PlotMergedFirstPrepPanel6.R" \
  --embedding "$DATA/FirstPrep_Merged_representative_embedding.csv" \
  --well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --well-stats-p3 "$WELLSTATS_DIR/FirstPrep_P3_well_stats.csv" \
  --well-stats-p4 "$WELLSTATS_DIR/FirstPrep_P4_well_stats.csv" \
  --threshold "$THRESHOLD" \
  --output "$OUTPUT/FirstPrep_Merged_panel6.png"

echo "=== QC-passed clustering: n_neighbors x min_dist sweep (QC-passed only, and all single cells) ==="
Rscript "$SCRIPTS/SweepQCThresholdClustering.R" \
  --matrix "$MERGED_MATRIX" \
  --well-empty-frac-csv "$DATA/FirstPrep_Merged_well_empty_frac.csv" \
  --well-stats-p3 "$WELLSTATS_DIR/FirstPrep_P3_well_stats.csv" \
  --well-stats-p4 "$WELLSTATS_DIR/FirstPrep_P4_well_stats.csv" \
  --threshold "$THRESHOLD" \
  --output-qc-passed-grid "$OUTPUT/FirstPrep_Merged_qc_threshold_clustering_sweep_passed.png" \
  --output-all-single-grid "$OUTPUT/FirstPrep_Merged_qc_threshold_clustering_sweep_all_single.png" \
  --output-qc-passed-pca "$DATA/FirstPrep_Merged_qcpassed_pca.csv.gz"

echo "Done."
