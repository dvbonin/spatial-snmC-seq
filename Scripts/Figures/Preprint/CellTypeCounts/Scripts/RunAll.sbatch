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
# CellTypeCounts/Logs must already exist -- it cannot be created from in here.

ROOT=${SNMC_FIGS}/Preprint/CellTypeCounts
SCRIPTS=${SNMC_ROOT}/Scripts/Figures/Preprint/CellTypeCounts/Scripts
R_ENV=${SNMC_ENV_R}

# Depends on the QC directory having written passes_cpg_qc into the well_stats
# tables (QC/Scripts/ApplyQCToWellStats.R). Without it the column is NA and the
# QC-passed filter matches nothing.

# The y range is derived in Shared/FirstPrepCoverage.R and matches the
# CovBoxplots QC-passed panel by construction.
"$R_ENV/bin/Rscript" "$SCRIPTS/PlotCellTypeCounts.R" \
  --output "$ROOT/Output/CellTypeCounts_FirstPrep.png"

echo "Done."
