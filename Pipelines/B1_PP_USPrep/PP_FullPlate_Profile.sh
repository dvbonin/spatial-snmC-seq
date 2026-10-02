#!/bin/bash
#SBATCH -J SNMC_B1_PP_FullPlate
#SBATCH -p cpu_p
#SBATCH --qos=cpu_preemptible
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3-00:00:00
#SBATCH -o Logs/controller_fullplate.out
#SBATCH -e Logs/controller_fullplate.err

source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

export TMPDIR="${SNMC_TMP}"
mkdir -p "$TMPDIR"

conda activate snakemake

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
