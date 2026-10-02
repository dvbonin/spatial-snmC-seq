#!/bin/bash
# Everything the figure directories share, prepared once. Nothing here draws a figure.
#
# The run order for rebuilding the preprint from scratch is:
#
#   1.  RunAll.sh pre-qc     what QC needs, plus the inputs that need no QC at all
#   2.  ../../QC             Leiden QC decision -> writes passes_cpg_qc into the
#                            Analysis well_stats tables, and the QC-passed PCA
#   3.  RunAll.sh post-qc    the single-cell pseudobulks, which pool QC-passed wells
#   4.  the figure dirs      CorrHeatmaps, PubSeq, ExpressionData, MarkerGenes,
#                            Microarray, ...
#
# The split exists because QC is not purely downstream of the pipelines: it clusters
# DataPrePrep's merged MethSCAn matrix, so that matrix has to be built before QC runs,
# while the single-cell pseudobulks can only be built after it. The bulk plate, the
# microarray and the gene annotation have no QC dependency either way and sit in pre-qc
# simply to get them out of the way early. The gene annotation is here because two
# figures (ExpressionData, MarkerGenes) read the same gene bodies and TSS positions.
#
# This is a submitter -- pure sbatch calls, no compute.

# Sourced before `set -u`: /etc/bashrc reads an unset BASHRCSOURCED and would abort.
source ~/.bashrc
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

STAGE=${1:-}
case "$STAGE" in pre-qc|post-qc) ;; *) echo "usage: RunAll.sh pre-qc|post-qc" >&2; exit 2;; esac

S=${SNMC_ROOT}/Scripts/Figures/Preprint/DataPrePrep/Scripts
W=${SNMC_FIGS}/Preprint/DataPrePrep
DS=${SNMC_DATA}/SNMC_B1
PP=${SNMC_PIPE}/PP_USPrep_All/Workup/Beds/Main/FirstPrep
WS=${SNMC_PIPE}/Analysis_Combined/Results/WellStats
ENV=${SNMC_ENV_R}
PQ_ENV=${SNMC_ENV_PQ}
RSCRIPT_ARRAY=${SNMC_ENV_ARRAY}/bin/Rscript
LIFTOVER=${SNMC_ENV_LIFTOVER}/bin/liftOver
CHAIN=${SNMC_ROOT}/Refs/Liftover/hg19ToHg38.over.chain.gz
GTF=${SNMC_REFS}/Annotation/GENCODE_hg38/gencode.v50.annotation.gtf.gz
REF=${SNMC_REFS}/Builds/Biscuit_hg38/hg38.analysisSet.fa
TMP=${SNMC_TMP}
SRC="$PP/P3=$WS/FirstPrep_P3_well_stats.csv,$PP/P4=$WS/FirstPrep_P4_well_stats.csv"

mkdir -p "$W/Data" "$W/Work" "$W/Logs"

submit() {  # submit <env> <name> <mem> <time> <deps|-> <command>
  local env=$1 name=$2 mem=$3 tl=$4 dep=$5; shift 5
  local d=(); [ "$dep" != "-" ] && d=(--dependency=afterok:$dep)
  sbatch --parsable --job-name="$name" -p cpu_p --qos=cpu_normal -c 2 --mem="$mem" -t "$tl" \
    "${d[@]}" -o "$W/Logs/${name}_%j.log" -e "$W/Logs/${name}_%j.log" \
    --wrap "source ~/.bashrc; set -euo pipefail; conda activate $env; export TMPDIR=$TMP; $*"
}
sub()   { submit "$ENV" "$@"; }       # R
subpy() { submit "$PQ_ENV" "$@"; }    # python

if [ "$STAGE" = pre-qc ]; then
  echo "== merged MethSCAn matrix over all FirstPrep wells (QC still undecided) =="
  echo "   $(sbatch --parsable $S/BuildMergedFirstPrepMethScan.sbatch)"

  echo "== GENCODE gene bodies and TSS positions =="
  echo "   $(subpy annotations 8G 00:20:00 - \
    "python3 $S/BuildGeneAnnotations.py --gtf $GTF \
     --bodies $W/Data/gene_bodies.bed --tss $W/Data/gene_tss.bed")"

  # The bulk plate is bismark output with no strand column, so collapsing its CpG
  # dyads needs the reference CpG set; both pools wait on it.
  echo "== reference CpG dyad coordinates =="
  J_CPG=$(subpy cpg_sites 8G 01:00:00 - \
    "python3 $S/BuildCpGSites.py --fasta $REF --output $W/Data/CpGSites_hg38.bed.gz")
  echo "   $J_CPG"

  echo "== bulk plate CHI002: wells with >4 cells, per cell type =="
  echo "   $(sub bulk_kz 64G 06:00:00 "$J_CPG" \
    "Rscript $S/BuildBulkPseudobulk.R --annotation $DS/BulkPlates/CHI002/CHI002_well_stats.csv \
     --cov-dir $DS/BulkPlates/CHI002/context_files_merged --cell-type keratinocyte --min-cells 5 \
     --cpg-sites $W/Data/CpGSites_hg38.bed.gz \
     --output $W/Data/Bulk_KZ.bed.gz")"
  echo "   $(sub bulk_mcc 64G 06:00:00 "$J_CPG" \
    "Rscript $S/BuildBulkPseudobulk.R --annotation $DS/BulkPlates/CHI002/CHI002_well_stats.csv \
     --cov-dir $DS/BulkPlates/CHI002/context_files_merged --cell-type merkel --min-cells 5 \
     --cpg-sites $W/Data/CpGSites_hg38.bed.gz \
     --output $W/Data/Bulk_MCC.bed.gz")"

  echo "== EPIC IDATs -> hg19 betas -> hg38 =="
  # minfi and the EPIC annotation only exist in SNMC_Array, and liftOver in its own env,
  # so both are called by absolute path rather than activated.
  echo "   $(sub array 32G 02:00:00 - \
    "$RSCRIPT_ARRAY $S/ExtractArrayBetas.R --idat-dir $DS/Microarrays \
       --output $W/Work/ArrayProbes_hg19.bed; \
     $LIFTOVER $W/Work/ArrayProbes_hg19.bed $CHAIN \
       $W/Data/ArrayProbes_hg38.bed $W/Work/ArrayProbes_unmapped.bed; \
     echo \"   lifted: \$(wc -l < $W/Data/ArrayProbes_hg38.bed)\"; \
     echo \"   unmapped: \$(( \$(wc -l < $W/Work/ArrayProbes_unmapped.bed) / 2 ))\"")"

  echo
  echo "next: ../../QC/Scripts/RunAll.sbatch, then RunAll.sh post-qc"
else
  echo "== QC-passed FirstPrep single cells -> per-CpG pseudobulk, two weightings =="
  # Reads passes_cpg_qc from the well_stats tables, which QC must have written first.
  #   reads  every read counts equally -- what the correlation figures compare to bulk
  #   wells  every well counts equally -- what the metagene figures average over cells
  for CT in KZ:HealthyKeratinocytes MCC:MerkelCarcinoma; do
    TAG=${CT%%:*} TYPE=${CT##*:}
    echo "   $(sub pseudobulk_reads_$TAG 128G 06:00:00 - \
      "Rscript $S/BuildSingleCellPseudobulk.R --well-stats '$SRC' --weighting reads \
       --cell-type $TYPE --output $W/Data/Pseudobulk_$TAG.bed.gz")"
    echo "   $(sub pseudobulk_wells_$TAG 128G 06:00:00 - \
      "Rscript $S/BuildSingleCellPseudobulk.R --well-stats '$SRC' --weighting wells \
       --cell-type $TYPE --output $W/Data/PseudobulkWellMean_$TAG.bed.gz")"
  done
  echo
  echo "next: the figure dirs"
fi
echo "watch: squeue -u \$USER"
