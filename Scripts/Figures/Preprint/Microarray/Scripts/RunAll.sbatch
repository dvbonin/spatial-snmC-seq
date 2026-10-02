#!/bin/bash
#SBATCH -J microarray
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o Logs/RunAll_%j.log
#SBATCH -e Logs/RunAll_%j.log

# hg38 array betas -> the UMAP panel coloured by each cell's methylation similarity to
# the two reference arrays. Three external inputs: the array betas prepared by
# ../DataPrePrep (IDATs -> minfi -> liftOver), the per-well CpG VCFs (PP pipeline), and
# the QC-passed single-cell PCA (../QC) -- the PCA is reused rather than recomputed so
# the cells plotted are exactly the cells the QC pass/fail call was made on.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"
conda activate ${SNMC_ENV_R}

S=${SNMC_ROOT}/Scripts/Figures/Preprint/Microarray/Scripts
FIG=${SNMC_FIGS}/Preprint
W=$FIG/Microarray
ARRAY=$FIG/DataPrePrep/Data/ArrayProbes_hg38.bed
VCF=${SNMC_PIPE}/PP_USPrep_All/Workup/VCF/Main/FirstPrep
PCA=$FIG/QC/Data/FirstPrep_Merged_qcpassed_pca.csv.gz

mkdir -p "$W/Data" "$W/Output"
export TMPDIR=${SNMC_TMP}

echo "== 1) rank by |beta_kc - beta_mc|, take top 1% =="
Rscript "$S/RankArraySitesByDiff.R" \
  --input "$ARRAY" \
  --output-ranked "$W/Data/RankedArraySitesByDiff.csv.gz" \
  --output-top-bed "$W/Data/TopDiffSites_Top1pct/TopDiffSites_array.bed.gz" \
  --percentile 1

echo "== 2) per-cell correlation to each reference, over those sites =="
# Wells come from the PCA, so correlations cover exactly the cells that get plotted.
for PLATE in P3 P4; do
  WELLS=$(zcat "$PCA" | awk -F, -v p="$PLATE" 'NR>1 && $2==p {printf "%s%s", sep, $3; sep=","}')
  echo "   $PLATE: $(echo "$WELLS" | tr ',' '\n' | wc -l) wells"
  python "$S/ComputeMicroarrayCorrelation.py" \
    --vcf-dir "$VCF/$PLATE" \
    --wells "$WELLS" \
    --topdiff-bed "$W/Data/TopDiffSites_Top1pct/TopDiffSites_array.bed.gz" \
    --threads 16 \
    --output "$W/Data/FirstPrep_${PLATE}_microarray_correlations.csv"
done

echo "== 3) UMAP panel =="
Rscript "$S/PlotUMAP.R" \
  --pca "$PCA" \
  --correlations-p3 "$W/Data/FirstPrep_P3_microarray_correlations.csv" \
  --correlations-p4 "$W/Data/FirstPrep_P4_microarray_correlations.csv" \
  --topdiff-label "Top 1% Diff." \
  --output "$W/Output/FirstPrep_Merged_umap.png" \
  --output-fused "$W/Output/FirstPrep_Merged_umap_fused.png"
