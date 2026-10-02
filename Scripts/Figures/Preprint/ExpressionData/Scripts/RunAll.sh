#!/bin/bash
# Methylation metagene profiles by bulk RNA-seq expression quartile, for QC-passed
# FirstPrep single cells of both cell types. Two figures:
#
#   promoter_metagene_by_quartile_5kb.png   TSS +/- 5kb, 40 bp bins
#   genebody_metagene_by_quartile.png       gene body scaled to 100 bins
#   expression_percentiles.png              the expression distribution those quartiles
#                                           are cut from, per cell type
#
# Inputs are the datasets tree (two public RNA-seq series) plus the gene regions and the
# well-mean single-cell pseudobulks prepared by ../../DataPrePrep, which must have been
# run before this script. Well-mean rather than read-weighted: a metagene
# by expression quartile is a statement about cells, so each cell contributes its call
# once regardless of how deeply it happened to be sequenced at that CpG.
#
# This is a submitter -- pure sbatch calls, no compute. Stages are chained by job
# dependency; the two pseudobulks and later the four profiles run concurrently.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/ExpressionData/Scripts
W=${SNMC_FIGS}/Preprint/ExpressionData
DPP=${SNMC_FIGS}/Preprint/DataPrePrep/Data
RNA=${SNMC_DATA}/SNMC_B1/RNAseq
PQ_ENV=${SNMC_ENV_PQ}
TMP=${SNMC_TMP}
FLANK=5000

mkdir -p "$W/Work" "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <env> <mem> <time> <deps|-> <command>
  local name=$1 env=$2 mem=$3 tl=$4 dep=$5; shift 5
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $env; export TMPDIR=$TMP; $*"
}

echo "== stage 1: expression table =="
J_EXP=$(sub expression "$PQ_ENV" 8G 00:15:00 - \
  "python3 $S/BuildExpression.py --mcc-counts $RNA/GSE223275_Raw_count_matrix.txt.gz \
   --kz-dir $RNA/GSE107871_extracted --output $W/Data/expression.csv")
echo "   $J_EXP"

echo "== stage 2: assign CpGs to gene bodies and TSS windows =="
PROF=()
for CT in KZ MCC; do
  PROF+=($(sub "genebody_$CT" "$PQ_ENV" 96G 02:00:00 - \
    "python3 $S/ProfileMetagene.py --pseudobulk $DPP/PseudobulkWellMean_$CT.bed.gz \
     --regions $DPP/gene_bodies.bed --mode genebody --n-bins 100 \
     --output $W/Data/genebody_$CT.parquet"))
  PROF+=($(sub "promoter_$CT" "$PQ_ENV" 96G 02:00:00 - \
    "python3 $S/ProfileMetagene.py --pseudobulk $DPP/PseudobulkWellMean_$CT.bed.gz \
     --regions $DPP/gene_tss.bed --mode promoter --flank $FLANK --bin-width 40 \
     --output $W/Data/promoter_$CT.parquet"))
done
echo "   ${PROF[*]}"

echo "== stage 3: figures =="
DEP=$(IFS=:; echo "${PROF[*]}"):$J_EXP
J_P1=$(sub plot_genebody "$PQ_ENV" 32G 00:30:00 "$DEP" \
  "python3 $S/PlotMetagene.py --mode genebody \
   --kz-profile $W/Data/genebody_KZ.parquet --mcc-profile $W/Data/genebody_MCC.parquet \
   --expression $W/Data/expression.csv \
   --output $W/Output/genebody_metagene_by_quartile.png")
J_P2=$(sub plot_promoter "$PQ_ENV" 32G 00:30:00 "$DEP" \
  "python3 $S/PlotMetagene.py --mode promoter --flank $FLANK \
   --kz-profile $W/Data/promoter_KZ.parquet --mcc-profile $W/Data/promoter_MCC.parquet \
   --expression $W/Data/expression.csv \
   --output $W/Output/promoter_metagene_by_quartile_5kb.png")
# Marked on the promoter gene set, the more inclusive of the two; passing the genebody
# profiles instead redraws it on that panel's slightly smaller set.
J_P3=$(sub plot_expression "$PQ_ENV" 16G 00:20:00 "$DEP" \
  "python3 $S/PlotExpressionPercentiles.py \
   --kz-profile $W/Data/promoter_KZ.parquet --mcc-profile $W/Data/promoter_MCC.parquet \
   --expression $W/Data/expression.csv \
   --output $W/Output/expression_percentiles.png")
echo "   $J_P1 $J_P2 $J_P3"
echo
echo "watch: squeue -u \$USER"
