#!/bin/bash
# Our assay against published single-cell methylation methods (Fig. 1e, Extended Data
# Fig. 1k-l):
#
#   AssayComparison_CpGsPerCell.png       covered CpGs per cell
#   AssayComparison_CpGsPerReadPair.png   covered CpGs per raw read pair
#   AssayComparison_ReadsVsCpGs.png       raw read pairs vs covered CpGs, log-log fits
#   AssayComparison_summary.tsv           the numbers quoted in the text
#
# Stage 1 picks the published cells and stage 2 counts them, straight from GEO, SRA and
# NeMO; files are streamed and discarded, only the per-cell tables are kept, in Data/.
# Stage 2 downloads ~25 GB and is the slow part. sci-MET v3 needs its NovaSeq amethyst
# object (4.2 GB), which is kept in Data/ so the metadata can be re-extracted.
#
# Our own cells come from the DataPrePrep well summaries, so this runs after
# DataPrePrep (pre-qc) and QC -- the PassesQC column must be filled.
#
# This is a submitter -- pure sbatch calls, no compute.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

S=${SNMC_ROOT}/Scripts/Figures/Preprint/AssayComparison/Scripts
W=${SNMC_FIGS}/Preprint/AssayComparison
WSUM=${SNMC_FIGS}/Preprint/DataPrePrep/Data
PQ_ENV=${SNMC_ENV_PQ}
R_ENV=${SNMC_ENV_R}
TMP=${SNMC_TMP}
AMETHYST=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8433nnn/GSM8433386/suppl/GSM8433386_NovaSeq.amethyst.Rdata.gz

mkdir -p "$W/Data" "$W/Output" "$W/Logs"

sub() {  # sub <name> <mem> <time> <deps|-> <env> <command>
  local name=$1 mem=$2 tl=$3 dep=$4 env=$5; shift 5
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 4 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $env; export TMPDIR=$TMP; $*"
}

echo "== stage 1: published cells per method =="
J_MAN=$(sub manifests 4G 00:30:00 - "$PQ_ENV" \
  "python3 $S/BuildManifests.py --out-dir $W/Data")
echo "   $J_MAN"

echo "== stage 2: covered CpGs and raw read pairs per published cell =="
J_CNT=()
for m in scBS-seq snmC-seq2 snmC-seq3; do
  J_CNT+=($(sub "count_$m" 16G 12:00:00 "$J_MAN" "$PQ_ENV" \
    "python3 $S/CountPublishedCells.py --manifest $W/Data/manifest_$m.tsv \
     --output $W/Data/percell_$m.tsv --jobs 4 --tmp $TMP"))
done
J_SM=$(sub scimet 48G 04:00:00 - "$R_ENV" \
  "cd $W/Data && curl -sSL --retry 5 -C - -O $AMETHYST && \
   gunzip -kf GSM8433386_NovaSeq.amethyst.Rdata.gz && \
   Rscript $S/ExtractAmethystMetadata.R --input GSM8433386_NovaSeq.amethyst.Rdata \
     --output $W/Data/percell_sciMET-v3.tsv && \
   conda activate $PQ_ENV && \
   python3 $S/SciMETReadBound.py --metadata $W/Data/percell_sciMET-v3.tsv \
     --output $W/Data/scimet_read_bound.tsv")
echo "   ${J_CNT[*]} $J_SM"

echo "== stage 3: figures =="
DEPS=$(IFS=:; echo "${J_CNT[*]}:$J_SM")
J_PLOT=$(sub plot 8G 00:20:00 "$DEPS" "$PQ_ENV" \
  "python3 $S/PlotAssayComparison.py \
   --well-summary P3=$WSUM/FirstPrep_P3_well_summary.csv P4=$WSUM/FirstPrep_P4_well_summary.csv \
   --scbs $W/Data/percell_scBS-seq.tsv --snmc2 $W/Data/percell_snmC-seq2.tsv \
   --snmc3 $W/Data/percell_snmC-seq3.tsv \
   --scimet-cells $W/Data/percell_sciMET-v3.tsv --scimet-bound $W/Data/scimet_read_bound.tsv \
   --out-dir $W/Output")
echo "   $J_PLOT"
echo
echo "watch: squeue -u \$USER"
