#!/bin/bash
#SBATCH -J MergedFPMethScan
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o Logs/MergedFPMethScan_%j.log
#SBATCH -e Logs/MergedFPMethScan_%j.log

# Instead of trying to reconcile two independently-called per-plate VMR sets
# (FirstPrep P3 and P4 share only ~4.5% of VMR coordinates exactly), pool both
# plates' single-cell beds and run MethSCAn's prepare/smooth/scan/matrix once
# on the combined ~768 cells, so VMRs are called jointly and every cell gets a
# real (non-imputed) value at every retained VMR. All FirstPrep wells are
# used here (not just QC-passed) -- QC filtering happens downstream, after
# the Leiden-clustering QC decision (see ../../QC).

ROOT=${SNMC_FIGS}/Preprint/DataPrePrep/Data/MergedFirstPrep
STAGING="$ROOT/Staged"
PREPARED="$ROOT/Prepared"
BEDS_FP=${SNMC_PIPE}/PP_USPrep_All/Workup/Beds/Main/FirstPrep

source ~/.bashrc
conda activate ${SNMC_ENV_METHSCAN}
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

rm -rf "$STAGING" "$PREPARED"
mkdir -p "$STAGING"
for PLATE in P3 P4; do
  for bed in "$BEDS_FP/$PLATE"/*_MethCPG.bed.gz; do
    well=$(basename "$bed" _MethCPG.bed.gz)
    ln -s "$bed" "$STAGING/FirstPrep_${PLATE}_${well}.bed.gz"
  done
done
echo "Staged $(ls "$STAGING" | wc -l) wells (P3 + P4)"

methscan prepare --input-format biscuit_short "$STAGING"/*.bed.gz "$PREPARED"
echo "Prepared."

methscan smooth "$PREPARED"
echo "Smoothed."

methscan scan --threads 32 "$PREPARED" "$ROOT/FirstPrep_Merged_VMRs.bed"
echo "Scanned VMRs."

mkdir -p "$ROOT/Matrix"
methscan matrix --threads 32 "$ROOT/FirstPrep_Merged_VMRs.bed" "$PREPARED" "$ROOT/Matrix"
echo "Wrote matrix outputs to $ROOT/Matrix"
