#!/bin/bash
# Site-level pairwise methylation correlation between three independent sources --
# a bulk plate, pseudobulk from QC-passed single cells, and a public bulk reference --
# for two cell types. Four figures:
#
#   CorrHeatmap_Sites_AllSites   r over every CpG both profiles cover
#   CorrHeatmap_Sites_Top1pct    r over the top 1% differential sites
#   SiteCount_Sites_AllSites     how many sites each of those correlations used
#   SiteCount_Sites_Top1pct      likewise for the top 1%
#
# No binning and no coverage filter on the correlations: a pair uses the CpGs both
# profiles measured, full stop. The one depth requirement is on SELECTION of the top
# 1% (>= MIN_COV reads in both bulk profiles), because at one or two reads a beta is
# 0 or 1 outright and |diff| hits its ceiling by sampling accident.
#
# The five profiles this correlates are prepared by ../../DataPrePrep, which must have
# been run (pre-qc and post-qc, with ../../QC in between) before this script. The only
# input read directly here is the public reference in the datasets tree.
#
# This is a submitter -- pure sbatch calls, no compute. Stages are chained by job
# dependency.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/CorrHeatmaps/Scripts
DPP=${SNMC_FIGS}/Preprint/DataPrePrep/Data
W=${SNMC_FIGS}/Preprint/CorrHeatmaps
DS=${SNMC_DATA}/SNMC_B1
ENV=${SNMC_ENV_R}
TMP=${SNMC_TMP}
MIN_COV=5
TOP_PCT=1

mkdir -p "$W/Work/SiteCorr" "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <command>
  local name=$1 mem=$2 tl=$3 dep=$4; shift 4
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $ENV; export TMPDIR=$TMP; $*"
}

echo "== stage 1: top ${TOP_PCT}% differential sites (>= ${MIN_COV} reads in both bulk profiles) =="
KEYS=$W/Data/TopDiffSites_Top${TOP_PCT}pct.tsv
J_KEYS=$(sub topsites 48G 02:00:00 - \
  "Rscript $S/TopDiffSites.R --bulk-kz $DPP/Bulk_KZ.bed.gz --bulk-mcc $DPP/Bulk_MCC.bed.gz \
   --top-pct $TOP_PCT --min-coverage $MIN_COV --output-keys $KEYS")
echo "   $J_KEYS"

echo "== stage 2: one job per matrix cell =="
MAN=$W/Work/pairs_manifest.tsv
: > "$MAN"
LABELS=(Bulk_KZ Bulk_MCC Pseudobulk_KZ Pseudobulk_MCC Public_KZ)
declare -A PATHS=(
  [Bulk_KZ]=$DPP/Bulk_KZ.bed.gz                [Bulk_MCC]=$DPP/Bulk_MCC.bed.gz
  [Pseudobulk_KZ]=$DPP/Pseudobulk_KZ.bed.gz    [Pseudobulk_MCC]=$DPP/Pseudobulk_MCC.bed.gz
  [Public_KZ]=$DS/PubSeq/Keratinocyte_public.bed.gz )
for VARIANT in AllSites Top${TOP_PCT}pct; do
  K="-"; [ "$VARIANT" != AllSites ] && K="$KEYS"
  for i in "${!LABELS[@]}"; do for j in "${!LABELS[@]}"; do
    [ "$j" -lt "$i" ] && continue          # unordered pairs only, diagonal included
    A=${LABELS[$i]}; B=${LABELS[$j]}
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$VARIANT" "$A" "${PATHS[$A]}" "$B" "${PATHS[$B]}" "$K" >> "$MAN"
  done; done
done
N=$(wc -l < "$MAN")
J_ARR=$(sbatch --parsable --dependency=afterok:$J_KEYS --array=1-$N \
  "$S/CorrelateSitePair.sbatch")
echo "   $N cells (15 pairs x 2 variants), array $J_ARR"

echo "== stage 3: assemble the tables and draw the four figures =="
J_PLOT=$(sub collect 8G 00:30:00 "$J_ARR" \
  "for V in AllSites Top${TOP_PCT}pct; do \
     head -1 \$(ls $W/Work/SiteCorr/\${V}__*.csv | head -1) > $W/Data/SiteCorrelationPairs_\${V}.csv; \
     tail -n +2 -q $W/Work/SiteCorr/\${V}__*.csv >> $W/Data/SiteCorrelationPairs_\${V}.csv; \
   done; bash $S/PlotOnly.sh")
echo "   $J_PLOT"
echo
echo "watch: squeue -u \$USER"
