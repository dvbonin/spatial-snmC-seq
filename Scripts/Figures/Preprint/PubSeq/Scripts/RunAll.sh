#!/bin/bash
# Single-cell pseudobulk vs a public bulk keratinocyte WGBS dataset, in 10kb bins.
#
# The single-cell pseudobulks are prepared by ../../DataPrePrep, which must have been
# run before this script. The only input read directly here is the public dataset in the
# datasets tree.
#
# This is a submitter -- pure sbatch calls, no compute.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/PubSeq/Scripts
DPP=${SNMC_FIGS}/Preprint/DataPrePrep/Data
W=${SNMC_FIGS}/Preprint/PubSeq
PUBLIC=${SNMC_DATA}/SNMC_B1/PubSeq/Keratinocyte_public.bed.gz
ENV=${SNMC_ENV_R}
TMP=${SNMC_TMP}
BIN=10000

mkdir -p "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <command>
  local name=$1 mem=$2 tl=$3 dep=$4; shift 4
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $ENV; export TMPDIR=$TMP; $*"
}

echo "== stage 1: bin everything to ${BIN}bp =="
BIN_JOBS=()
BIN_JOBS+=($(sub bin_kz 32G 01:00:00 - \
  "Rscript $S/BinMethylation.R --input $DPP/Pseudobulk_KZ.bed.gz \
   --bin-size $BIN --output $W/Data/Keratinocytes_${BIN}.bed.gz"))
BIN_JOBS+=($(sub bin_mcc 32G 01:00:00 - \
  "Rscript $S/BinMethylation.R --input $DPP/Pseudobulk_MCC.bed.gz \
   --bin-size $BIN --output $W/Data/MerkelCarcinoma_${BIN}.bed.gz"))
BIN_JOBS+=($(sub bin_public 32G 01:00:00 - \
  "Rscript $S/BinMethylation.R --input $PUBLIC \
   --bin-size $BIN --output $W/Data/Keratinocyte_public_${BIN}.bed.gz"))
echo "   ${BIN_JOBS[*]}"

echo "== stage 2: scatter =="
DEP=$(IFS=:; echo "${BIN_JOBS[*]}")
J_PLOT=$(sub scatter 16G 00:30:00 "$DEP" \
  "Rscript $S/PlotPublicScatter.R \
   --panels 'Keratinocytes (Pseudobulk)=$W/Data/Keratinocytes_${BIN}.bed.gz,Merkel Cell Carcinoma (Pseudobulk)=$W/Data/MerkelCarcinoma_${BIN}.bed.gz' \
   --y-label 'Keratinocytes (Public)' \
   --y-file $W/Data/Keratinocyte_public_${BIN}.bed.gz \
   --min-coverage 50 \
   --title 'Pseudobulk vs. Public Keratinocyte (10kb Bins)' \
   --output $W/Output/SingleCell_vs_Public_10kbBins.png")
echo "   $J_PLOT"
echo
echo "watch: squeue -u \$USER"
