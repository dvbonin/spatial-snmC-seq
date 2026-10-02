#!/bin/bash
#SBATCH -J sitecorr
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o Logs/sitecorr_%A_%a.log
#SBATCH -e Logs/sitecorr_%A_%a.log

# One array task = one cell of one matrix. The manifest line gives the variant, the
# two profiles and the optional site restriction.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"
conda activate ${SNMC_ENV_R}

S=${SNMC_ROOT}/Scripts/Figures/Preprint/CorrHeatmaps/Scripts
W=${SNMC_FIGS}/Preprint/CorrHeatmaps
export TMPDIR=${SNMC_TMP}

read -r VARIANT LA FA LB FB KEYS < <(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i' "$W/Work/pairs_manifest.tsv")
OUT="$W/Work/SiteCorr/${VARIANT}__${LA}__${LB}.csv"
mkdir -p "$(dirname "$OUT")"

# A finished cell is exactly a header plus one row. Checking the row count rather
# than mere existence means a cell interrupted mid-write is redone, not trusted.
if [ -f "$OUT" ] && [ "$(wc -l < "$OUT")" -eq 2 ]; then
  echo "cached: $(basename "$OUT")"; exit 0
fi

ARGS=(--label-a "$LA" --file-a "$FA" --label-b "$LB" --file-b "$FB" --output "$OUT")
[ "$KEYS" != "-" ] && ARGS+=(--keys "$KEYS")
Rscript "$S/CorrelateSitePair.R" "${ARGS[@]}"
