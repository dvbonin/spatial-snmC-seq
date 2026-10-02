#!/bin/bash
# Spike-in controls across both FirstPrep plates:
#
#   FirstPrep_spikein_scatter.png        Lambda conversion vs pUC19 methylation, one
#                                        dot per well, coloured by plate
#   FirstPrep_Lambda_coverage_box.png    Lambda coverage, coloured by conversion
#   FirstPrep_pUC19_coverage_box.png     pUC19 coverage, coloured by methylation
#
# Inputs are the per-well spike-in VCFs (PP pipeline) and the Analysis well_stats.
# Independent of QC and of DataPrePrep: every well of every plate is plotted, since the
# point is the controls across the whole plate.
#
# This is a submitter -- pure sbatch calls, no compute.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/SpikeIn/Scripts
W=${SNMC_FIGS}/Preprint/SpikeIn
V=${SNMC_PIPE}/PP_USPrep_All/Workup/VCF
WS=${SNMC_PIPE}/Analysis_Combined/Results/WellStats
PQ_ENV=${SNMC_ENV_PQ}
TMP=${SNMC_TMP}
TABLE=$W/Data/spikein_per_well.tsv

mkdir -p "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <command>
  local name=$1 mem=$2 tl=$3 dep=$4; shift 4
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $PQ_ENV; export TMPDIR=$TMP; $*"
}

echo "== stage 1: per-well beta and coverage from both spike-in VCFs =="
J_TAB=$(sub collect 8G 01:00:00 - \
  "python3 $S/CollectSpikeIn.py --output $TABLE \
   --plates 'P3=$V/Lambda/FirstPrep/P3,$V/pUC19/FirstPrep/P3,$WS/FirstPrep_P3_well_stats.csv' \
            'P4=$V/Lambda/FirstPrep/P4,$V/pUC19/FirstPrep/P4,$WS/FirstPrep_P4_well_stats.csv'")
echo "   $J_TAB"

echo "== stage 2: figures =="
J1=$(sub plot_scatter 8G 00:20:00 "$J_TAB" \
  "python3 $S/PlotSpikeInScatter.py --table $TABLE \
   --output $W/Output/FirstPrep_spikein_scatter.png")
J2=$(sub plot_boxes 8G 00:20:00 "$J_TAB" \
  "python3 $S/PlotSpikeInCoverageBox.py --table $TABLE \
   --output-lambda $W/Output/FirstPrep_Lambda_coverage_box.png \
   --output-puc19 $W/Output/FirstPrep_pUC19_coverage_box.png")
echo "   $J1 $J2"
echo
echo "watch: squeue -u \$USER"
