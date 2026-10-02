#!/bin/bash
#SBATCH -J bulkcov
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o Logs/bulkcov_%j.log
#SBATCH -e Logs/bulkcov_%j.log

# Tabulates per-site depth with awk rather than loading 28M rows into R; the small
# table it writes is the intermediate the plot reads, so restyling costs nothing.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"
conda activate ${SNMC_ENV_R}
S=${SNMC_ROOT}/Scripts/Figures/Preprint/CorrHeatmaps/Scripts
W=${SNMC_FIGS}/Preprint/CorrHeatmaps
OUT=$W/Data/BulkCoverageDistribution.csv

if [ ! -s "$OUT" ]; then
  echo "sample,coverage,n_sites" > "$OUT"
  for L in Bulk_KZ Bulk_MCC; do
    zcat "$W/Work/$L.bed.gz" | awk -F'\t' -v s="$L" '{c[$5]++} END{for(k in c) print s","k","c[k]}' >> "$OUT"
  done
fi
echo "tabulated: $(( $(wc -l < "$OUT") - 1 )) rows"

Rscript "$S/PlotBulkCoverage.R" --input "$OUT" --mark 3,5 --max-cov 30 \
  --output "$W/Output/BulkCoverageDistribution.png"
