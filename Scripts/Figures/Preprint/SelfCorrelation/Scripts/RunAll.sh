#!/bin/bash
# Per-CpG pairwise methylation discordance between single cells, keratinocytes vs
# Merkel carcinoma -- following up why MCC agrees less well with its own bulk.
#
#   stage 1  per site, per cell type: how many cells cover it and how many call it
#            methylated (the two cell types run in parallel)
#   stage 2  match the two types site by site to the smaller cell count, then compute
#            D = 2MU/(k(k-1)) for each
#   stage 3  three panels
#
# Inputs are only the per-well beds and well_stats from the PP/Analysis pipelines.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/SelfCorrelation/Scripts
W=${SNMC_FIGS}/Preprint/SelfCorrelation
PP=${SNMC_PIPE}/PP_USPrep_All/Workup/Beds/Main/FirstPrep
WS=${SNMC_PIPE}/Analysis_Combined/Results/WellStats
ENV=${SNMC_ENV_R}
TMP=${SNMC_TMP}
SRC="$PP/P3=$WS/FirstPrep_P3_well_stats.csv,$PP/P4=$WS/FirstPrep_P4_well_stats.csv"

mkdir -p "$W/Work" "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <command>
  local name=$1 mem=$2 tl=$3 dep=$4; shift 4
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $ENV; export TMPDIR=$TMP; $*"
}

echo "== stage 1: per-site cell counts and methylated-cell counts =="
J_KZ=$(sub agg_kz 96G 04:00:00 - \
  "Rscript $S/AggregateCellStates.R --well-stats '$SRC' --cell-type HealthyKeratinocytes \
   --output $W/Work/CellStates_KZ.tsv.gz")
J_MC=$(sub agg_mcc 96G 04:00:00 - \
  "Rscript $S/AggregateCellStates.R --well-stats '$SRC' --cell-type MerkelCarcinoma \
   --output $W/Work/CellStates_MCC.tsv.gz")
echo "   $J_KZ $J_MC"

echo "== stage 2: match cell counts and compute discordance =="
J_D=$(sub discordance 96G 03:00:00 "$J_KZ:$J_MC" \
  "Rscript $S/ComputeDiscordance.R --kz $W/Work/CellStates_KZ.tsv.gz \
   --mcc $W/Work/CellStates_MCC.tsv.gz --min-cells 2 \
   --output-sites $W/Work/PerSiteCellCounts.tsv.gz \
   --output-hist $W/Data/DiscordanceHistogram.csv \
   --output-summary $W/Data/DiscordanceSummary.csv \
   --output-cellspersite $W/Data/CellsPerSite.csv")
echo "   $J_D"

echo "== stage 3: plots =="
J_P=$(sub plot 16G 00:30:00 "$J_D" \
  "Rscript $S/PlotDiscordance.R --hist $W/Data/DiscordanceHistogram.csv \
   --summary $W/Data/DiscordanceSummary.csv --cells-per-site $W/Data/CellsPerSite.csv \
   --min-k 2 --max-k-plot 20 --min-sites 1000 \
   --output $W/Output/Discordance.png")
echo "   $J_P"
echo
echo "watch: squeue -u \$USER"
