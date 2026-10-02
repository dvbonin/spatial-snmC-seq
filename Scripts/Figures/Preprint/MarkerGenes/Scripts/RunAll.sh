#!/bin/bash
# Per-cell mean methylation at five marker gene loci, two region definitions:
#
#   marker_gene_violins_body.png        per-cell mean over the whole gene body
#   marker_gene_violins_promoter.png    per-cell mean over TSS-2kb to TSS+500bp
#
# Both region definitions were recovered by testing candidates against the surviving
# per-cell tables of the original figure -- see PlotMarkerGeneViolins.py.
#
# Inputs are the per-well beds (PP pipeline), the QC verdict in the Analysis well_stats,
# and the gene bodies prepared by ../../DataPrePrep -- so DataPrePrep pre-qc and QC must
# both have run first. Per-cell data cannot come from a pseudobulk, so the beds are read
# directly, by region.
#
# This is a submitter -- pure sbatch calls, no compute.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/MarkerGenes/Scripts
W=${SNMC_FIGS}/Preprint/MarkerGenes
DPP=${SNMC_FIGS}/Preprint/DataPrePrep/Data
PP=${SNMC_PIPE}/PP_USPrep_All/Workup/Beds/Main/FirstPrep
WS=${SNMC_PIPE}/Analysis_Combined/Results/WellStats
PQ_ENV=${SNMC_ENV_PQ}
TMP=${SNMC_TMP}
SRC="$PP/P3=$WS/FirstPrep_P3_well_stats.csv,$PP/P4=$WS/FirstPrep_P4_well_stats.csv"

# KRT5/KRT14 are the keratinocyte markers, KRT20/SYP/CHGA the Merkel ones. The
# violins put the Merkel markers left of the divider and the keratinocyte ones right.
MCC_GENES=KRT20,SYP,CHGA
KZ_GENES=KRT5,KRT14
UPSTREAM=2000
PROMOTER_UP=2000
PROMOTER_DOWN=500

mkdir -p "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <command>
  local name=$1 mem=$2 tl=$3 dep=$4; shift 4
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $PQ_ENV; export TMPDIR=$TMP; $*"
}

echo "== stage 1: pull every QC-passed cell's CpGs over the five loci =="
J_OBS=$(sub observations 16G 01:00:00 - \
  "python3 $S/ProfileMarkerGenes.py --well-stats '$SRC' --bodies $DPP/gene_bodies.bed \
   --genes $MCC_GENES,$KZ_GENES --upstream $UPSTREAM \
   --output $W/Data/marker_gene_observations.parquet")
echo "   $J_OBS"

echo "== stage 2: violin figures =="
J_VIO=()
# region, promoter flank (ignored for body), output suffix
for SPEC in body:body promoter:promoter; do
  IFS=: read -r REGION SUFFIX <<< "$SPEC"
  J_VIO+=($(sub "violins_$SUFFIX" 16G 00:20:00 "$J_OBS" \
    "python3 $S/PlotMarkerGeneViolins.py \
     --observations $W/Data/marker_gene_observations.parquet --bodies $DPP/gene_bodies.bed \
     --region $REGION --promoter-upstream $PROMOTER_UP --promoter-downstream $PROMOTER_DOWN \
     --left-genes $MCC_GENES --right-genes $KZ_GENES \
     --output $W/Output/marker_gene_violins_${SUFFIX}.png"))
done
echo "   ${J_VIO[*]}"
echo
echo "watch: squeue -u \$USER"
