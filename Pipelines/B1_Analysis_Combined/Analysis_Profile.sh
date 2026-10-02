#!/bin/bash
#SBATCH -J SNMC_B1_Analysis
#SBATCH -p cpu_p
#SBATCH --qos=cpu_preemptible
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3-00:00:00
#SBATCH -o Logs/controller_analysis.out
#SBATCH -e Logs/controller_analysis.err

source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

conda activate snakemake

snakemake \
  --directory ${SNMC_PIPE}/Analysis_Combined \
  --snakefile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Analysis.smk \
  --configfile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Configs/Config_Analysis.yaml \
  --workflow-profile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined \
  --slurm-logdir ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Logs/ \
  --cores 512 \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --latency-wait 60 \
  --scheduler greedy \
  --conda-frontend conda
