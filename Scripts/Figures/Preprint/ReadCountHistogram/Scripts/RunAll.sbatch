#!/bin/bash
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH -o Logs/RunAll_%j.log
#SBATCH -e Logs/RunAll_%j.log
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

# R and python both put scratch files under TMPDIR; without this they fall back to
# /tmp on the compute node, which fills up and fails the job with errors like
# "cannot create 'R_TempDir'".
export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

# NOTE: SLURM opens the -o/-e log file before this script runs, so
# ReadCountHistogram/Logs must already exist -- it cannot be created from here.
#
# No QC dependency: this figure shows how well read count alone separates empty
# from single-cell wells, so it deliberately uses no QC filter.

ROOT=${SNMC_FIGS}/Preprint/ReadCountHistogram
SCRIPTS=${SNMC_ROOT}/Scripts/Figures/Preprint/ReadCountHistogram/Scripts
R_ENV=${SNMC_ENV_R}

# Own x axis rather than one shared across plate groups: the preps differ ~10x in
# median yield, so a shared axis would squash the deeper one.
echo "=== FirstPrep P3+P4 ==="
"$R_ENV/bin/Rscript" "$SCRIPTS/PlotReadCountHistogram.R" \
  --plates "FirstPrep_P3,FirstPrep_P4" \
  --output "$ROOT/Output/ReadCountHistogram_FirstPrep.png"

echo "Done."
