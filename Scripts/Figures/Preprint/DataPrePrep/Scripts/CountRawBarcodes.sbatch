#!/bin/bash
#SBATCH -p cpu_p
#SBATCH --qos=cpu_normal
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH -o Logs/rawbc_%A_%a.log
#SBATCH -e Logs/rawbc_%A_%a.log
set -euo pipefail
: "${SNMC_ROOT:?source the repository config.sh first}"

# Raw read pairs per cell barcode, straight from the untrimmed FASTQ -- the one
# per-well number the pipeline cannot report, because barcodes are only assigned
# after alignment. One array task per sequencing chunk.
#
# R1 only: the barcode is part of the read name, which both mates share, so R1
# gives read PAIRS and halves the I/O. Raw records = 2x these counts.
#
# igzip (Intel ISA-L) rather than zcat -- a single gzip stream cannot be inflated
# in parallel, so decompression is one core per file and the parallelism has to
# come from running the chunks concurrently.

IN=${SNMC_DATA}/SNMC_B1/PrepUS
OUT=${SNMC_FIGS}/Preprint/DataPrePrep/Data/RawBarcodeCounts
IGZIP=${SNMC_ENV_BIO}/bin/igzip
mkdir -p "$OUT"

R1=$(ls "$IN"/FirstPrep_S*_P[34]_C*_R1.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
CHUNK=$(basename "$R1" _R1.fastq.gz)
echo "chunk $SLURM_ARRAY_TASK_ID: $CHUNK"

# Header lines are every 4th. The barcode is the last ':'-separated token of the
# first whitespace field, always 8 bp, so take it by offset.
"$IGZIP" -dc "$R1" \
| gawk 'NR % 4 == 1 { c[substr($1, length($1) - 7)]++ }
        END { for (b in c) print b "\t" c[b] }' \
> "$OUT/${CHUNK}.tsv"

echo "distinct barcodes: $(wc -l < "$OUT/${CHUNK}.tsv")"
echo "read pairs:        $(awk -F'\t' '{s+=$2} END{print s}' "$OUT/${CHUNK}.tsv")"
