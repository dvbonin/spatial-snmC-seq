#!/bin/bash
#SBATCH -J SNMC_B1_FullPlate_ThenAnalysis
#SBATCH -p cpu_p
#SBATCH --qos=cpu_preemptible
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3-00:00:00
#SBATCH -o Logs/controller_fullplate_then_analysis.out
#SBATCH -e Logs/controller_fullplate_then_analysis.err

source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

conda activate snakemake

# 1. Finish FullPlate (first half) -- must succeed before Analysis starts, since Analysis
# reads FirstPrep/P3 and FirstPrep/P4's merged wells from its output.
snakemake \
  --directory ${SNMC_PIPE}/PP_USPrep_All \
  --snakefile ${SNMC_ROOT}/Pipelines/B1_PP_USPrep/PP_FullPlate.smk \
  --configfile ${SNMC_ROOT}/Pipelines/B1_PP_USPrep/Configs/Config_FullPlate.yaml \
  --workflow-profile ${SNMC_ROOT}/Pipelines/B1_PP_USPrep \
  --slurm-logdir ${SNMC_ROOT}/Pipelines/B1_PP_USPrep/Logs/ \
  --cores 512 \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --retries 10 \
  --keep-going \
  --latency-wait 60 \
  --scheduler greedy \
  --conda-frontend conda

# 2. Run Analysis (second half) on all configured plates, now including FirstPrep/P3 & P4.
snakemake \
  --directory ${SNMC_PIPE}/Analysis_Combined \
  --snakefile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Analysis.smk \
  --configfile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Configs/Config_Analysis.yaml \
  --workflow-profile ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined \
  --slurm-logdir ${SNMC_ROOT}/Pipelines/B1_Analysis_Combined/Logs/ \
  --cores 512 \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --retries 10 \
  --keep-going \
  --latency-wait 60 \
  --scheduler greedy \
  --conda-frontend conda
