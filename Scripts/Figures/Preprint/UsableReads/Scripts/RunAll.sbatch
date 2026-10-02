#!/bin/bash
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH -o Logs/RunAll_%j.log
#SBATCH -e Logs/RunAll_%j.log
# ~/.bashrc must be sourced before `set -u`: /etc/bashrc reads an unset
# variable, which -u turns into a fatal error.
source ~/.bashrc
conda activate ${SNMC_ENV_PQ}
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

# R and python both put scratch files under TMPDIR; without this they fall back to
# /tmp on the compute node, which fills up and fails the job with errors like
# "cannot create 'R_TempDir'".
export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

# NOTE: SLURM opens the -o/-e log file before this script runs, so
# UsableReads/Logs must already exist -- it cannot be created from in here.

ROOT=${SNMC_FIGS}/Preprint/UsableReads
SCRIPTS=${SNMC_ROOT}/Scripts/Figures/Preprint/UsableReads/Scripts
OUTPUT="$ROOT/Output"

RC_USPREP=${SNMC_PIPE}/PP_USPrep_All/Results/ReadCounts

# Two variants: the three-stage headline funnel, and every
# stage the pipeline records (see PlotUsableReads.py for why the latter is
# grouped rather than overlaid).
for STAGES in main all; do
  case "$STAGES" in
    main) SUFFIX="" ;;
    all)  SUFFIX="_allstages" ;;
  esac

  echo "=== FirstPrep P3+P4 (${STAGES}) ==="
  python3 "$SCRIPTS/PlotUsableReads.py" \
    --read-counts "P3=$RC_USPREP/FirstPrep/P3/Main_read_counts.csv" \
                  "P4=$RC_USPREP/FirstPrep/P4/Main_read_counts.csv" \
    --stages "$STAGES" --title "Read Yield" \
    --output "$OUTPUT/usable_reads_firstprep${SUFFIX}.png"
done

echo "Done."
