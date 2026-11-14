# Spatial snmC-seq — analysis code

Analysis code for the manuscript "Spatially resolved single-nucleus DNA methylation profiling in FFPE tissue" (von Bonin et al.).

    git clone https://github.com/dvbonin/spatial-snmC-seq

Single-nucleus methylome sequencing (snmC-seq3) of laser-capture-microdissected
Merkel cell carcinoma and keratinocyte nuclei, from raw reads through to the
published figures.

## What is here

| Path | Contents |
|---|---|
| `Pipelines/B1_PP_USPrep/` | Preprocessing pipeline (Snakemake): trimming, alignment, barcode correction, deduplication, per-well demultiplexing, methylation calling. |
| `Pipelines/B1_Analysis_Combined/` | Analysis pipeline (Snakemake): per-well statistics, MethSCAn region calling, methylation matrices, PCA, clustering. Configured for the two plates used in this paper. |
| `Scripts/Python/`, `Scripts/R/` | Scripts called by the two pipelines. |
| `Scripts/Figures/Preprint/` | One directory per figure or figure group, each with its own `RunAll` submitter. |
| `Envs/` | Conda environment definitions (`*.yaml`) and the exact resolved package versions used (`resolved_*.txt`). |

## Order of execution

The figure directories depend on the pipelines, and the quality-control decision
gates most of them, so the order matters:

1. `Pipelines/B1_PP_USPrep` — preprocessing, per plate
2. `Pipelines/B1_Analysis_Combined` — per-well statistics and methylation matrices
3. `Scripts/Figures/Preprint/DataPrePrep/Scripts/RunAll.sh pre-qc` — shared inputs: pooled MethSCAn matrix, gene annotations, reference CpG coordinates, bulk-plate and array data, raw read pairs per well
4. `Scripts/Figures/Preprint/QC/Scripts/RunAll.sbatch` — the quality-control decision, written back into the per-well tables
5. `Scripts/Figures/Preprint/DataPrePrep/Scripts/RunAll.sh post-qc` — single-cell pseudobulks, which pool only QC-passed wells, and the per-plate well summaries
6. The remaining figure directories, in any order

## Figures

| Directory | Figure content |
|---|---|
| `AssayComparison` | CpG coverage and read efficiency against published single-cell methylation assays |
| `CellTypeCounts` | CpG coverage per cell, by cell type |
| `CorrHeatmaps` | Correlation between single-cell pseudobulk, bulk plate, and public data |
| `CovBoxplots` | CpG coverage by number of cells per well |
| `ExpressionData` | Methylation metagenes by expression quartile |
| `MarkerGenes` | Promoter and gene-body methylation at marker genes |
| `Microarray` | Per-cell correlation to EPIC array references |
| `PubSeq` | Pseudobulk versus published bulk keratinocyte methylome |
| `QC` | Quality-control panels and threshold sweeps |
| `ReadCountHistogram` | Read counts per well |
| `SelfCorrelation` | Pairwise agreement between cells of the same type |
| `SpikeIn` | Lambda and pUC19 spike-in conversion controls |
| `UsableReads` | Read survival through the pipeline |

## Configuration

All machine-specific paths live in one file. Copy the template, set the values for
your system, and source it before running anything:

    cp config.sh.example config.sh
    # edit config.sh
    source config.sh

It defines where the code lives, where input data is read from, where output is
written, the reference and annotation locations, scratch space, and the conda
environments. SLURM exports these to the jobs, so sourcing it once per shell is
enough.

The Snakemake pipelines take their own configuration the same way:

    cp Pipelines/B1_PP_USPrep/Configs/Config_FullPlate.yaml.example \
       Pipelines/B1_PP_USPrep/Configs/Config_FullPlate.yaml
    cp Pipelines/B1_Analysis_Combined/Configs/Config_Analysis.yaml.example \
       Pipelines/B1_Analysis_Combined/Configs/Config_Analysis.yaml

Nothing outside `config.sh` and those two files needs editing.

Environments are defined in `Envs/*.yaml`; the `resolved_*.txt` files alongside
them record the exact package versions used, should an exact rebuild be needed.

Jobs are submitted with `sbatch`; running under a different scheduler means
adapting those calls. The code is published as the record of what produced the
figures rather than as a general-purpose tool.

## Published methods in the assay comparison

`AssayComparison` downloads and counts the published data itself (GEO, SRA/ENA and
the NeMO Archive); nothing has to be staged by hand. Per method:

| Method | Source | Cells |
|---|---|---|
| scBS-seq (Smallwood et al. 2014) | GSE56879 CpG reports, SRP041257 | all single cells; controls, bulk samples and the two resequenced "(deeper)" libraries excluded |
| snmC-seq2 (Luo et al. 2018) | GSE112471 allc files, SRP136737 | every 4th sample in GSM order |
| snmC-seq3 (Tian et al. 2023) | NeMO `Ecker_Human_Methylome_2023` allc and FASTQ | 60 cells, seeded random sample |
| sci-MET v3 (Nichols et al. 2025) | GSE273592 NovaSeq amethyst object, ENA run counts | all 44,840 cells |

CpGs are strand-merged dinucleotides and reads are raw read pairs for every method.
sci-MET v3 deposits no per-cell read counts; it enters the read-efficiency panels as a
bound (total read pairs over cells), not as data.

## Data

Sequencing data: [accession, to be added]. Public datasets used:
GSE223275 and GSE107871 (bulk RNA-seq), and the published bulk keratinocyte
methylome and EPIC array data cited in the manuscript.

## Citation

[citation, to be added once the preprint is posted]

Repository: https://github.com/dvbonin/spatial-snmC-seq
