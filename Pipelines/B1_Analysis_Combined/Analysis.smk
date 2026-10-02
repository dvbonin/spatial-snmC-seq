##############################################################################################
#----------------------------------  Analysis (second half)  --------------------------------#
##############################################################################################
# Analysis.smk -- shared by both first-half variants, and combined across however many of them
# are configured. Starts from the merged per-well BAMs / wells.tsv / read-counts produced by
# Pipelines/B1_PP_USPrep/PP_FullPlate.smk, and optionally the quadrant-barcoded variant (both
# emit the same 384-well output format). config["Sources"] names each first-half OutDir once and
# config["Samples"] says, per sample, which of them produced it, which plates to run and which
# Mapping CSV(s) annotate its wells. All PP paths are built from that by pp().
# All outputs from this file go to the single shared config["OutDir"], independent of which
# source a given (sample, plate) came from (see Configs/Config_Analysis.yaml).

import os
import csv

##############################################################################################
### Configs

WorkingDir = config["WorkingDir"]
OutDir     = config["OutDir"]

_spikes      = config.get("SpikeIns", {})
Spike_UnMeth = _spikes.get("Unmethylated", None)
Spike_Meth   = _spikes.get("Methylated",   None)
SpikeInRef   = _spikes.get("CombinedRef",  None)
UseSpikes    = bool(Spike_UnMeth and Spike_Meth and SpikeInRef)

PythonScripts = WorkingDir + "/Scripts/Python"
RScripts      = WorkingDir + "/Scripts/R"

# This half never parses raw fastq, so (sample, plate) pairs can't be discovered on their own
# -- config["Samples"] lists them explicitly, along with which first-half run produced each
# sample and which Mapping CSV(s) annotate its wells.
SAMPLE_PLATES = [
    (sample, f"P{n}")
    for sample, cfg in config["Samples"].items()
    for n in cfg["Plates"]
]

def _mappings_for(wildcards):
    return config["Samples"][wildcards.sample]["Mappings"]


def pp(sample, plate, what, type="Main"):
    """Path to a first-half (PP) artifact -- the one place that knows their layout.

    Only the root is configurable (config["Sources"], selected per sample); everything
    below it is fixed by PP_FullPlate.smk / PP_Quadrants.smk, so it is defined here once
    instead of being spelled out at every input. Four families of artifact cross the
    pipeline boundary: per-stream read counts, the per-context methylation summary, the
    per-well CpG beds, and the spike-in VCFs. Main-genome VCFs never cross.
    """
    root = config["Sources"][config["Samples"][sample]["Source"]]
    return {
        "read_counts":  f"{root}/Results/Split/{sample}/{plate}/{type}_read_counts.csv",
        "meth_summary": f"{root}/Results/MethStats/{sample}_{plate}_MethStats_Summary.csv",
        "beds_done":    f"{root}/Workup/Beds/Main/{sample}/{plate}/.BedsDone",
        "beds_dir":     f"{root}/Workup/Beds/Main/{sample}/{plate}",
        "spikein_done": f"{root}/Workup/VCF/{type}/{sample}/{plate}/.VCFDone",
        "spikein_dir":  f"{root}/Workup/VCF/{type}/{sample}/{plate}",
    }[what]

wildcard_constraints:
    sample      = r"[A-Za-z0-9][A-Za-z0-9-]*",
    plate       = r"P\d+",
    well        = r"[A-P]\d{2}",
    type        = r"Main|Lambda|pUC19",
    matrix_type = r"methscan|bin100k",


##############################################################################################
### Outputs

_SI_OUTPUTS = [
    ("Lambda",   "well_reads_heatmap.png"),
    ("Lambda",   "well_methylation_heatmap.png"),
    ("pUC19",    "well_reads_heatmap.png"),
    ("pUC19",    "well_methylation_heatmap.png"),
    ("Combined", "lambda_conversion_all.png"),
    ("Combined", "puc19_methylation_cpg.png"),
    ("Combined", "lambdaALL_vs_puc19CPG.png"),
]

MATRIX_TYPES = ["methscan", "bin100k"]

OUT_clustering_panels = [
    OutDir + f"/Plots/Clustering/{s}/{p}/{matrix_type}_panel.png"
    for s, p in SAMPLE_PLATES for matrix_type in MATRIX_TYPES
]

OUT_leiden = (
    [OutDir + f"/Plots/Leiden/{s}/{p}/{mt}_umap_clusters.png"
     for s, p in SAMPLE_PLATES for mt in MATRIX_TYPES]
    + [OutDir + f"/Plots/Leiden/{s}/{p}/{mt}_cpg_coverage_by_cluster.png"
       for s, p in SAMPLE_PLATES for mt in MATRIX_TYPES]
)

OUT_cpg_coverage = (
    [OutDir + f"/Plots/CpG_Coverage/{s}_{p}_boxplot_cpg_coverage.png" for s, p in SAMPLE_PLATES]
    + [OutDir + f"/Plots/MethylationRates/{s}_{p}_boxplot_methylation_rates.png" for s, p in SAMPLE_PLATES]
)

OUT_spikein = (
    [OutDir + f"/Plots/SpikeIn/{s}/{p}/{subdir}/{fname}"
     for s, p in SAMPLE_PLATES for subdir, fname in _SI_OUTPUTS]
    if UseSpikes else []
)


##############################################################################################
### All rule

rule all:
    input:
        OUT_clustering_panels
        + OUT_leiden
        + OUT_cpg_coverage
        + OUT_spikein


##############################################################################################

rule build_well_stats:
    # Consolidated per-well stats table -- every downstream rule reads this. The
    # methylation columns come straight from the first-half pipeline's
    # MethStats_Summary.csv (WellMethStats.py); this half no longer rescans the
    # per-well VCFs to recompute numbers that already exist.
    input:
        mappings      = _mappings_for,
        main_counts   = lambda wc: pp(wc.sample, wc.plate, "read_counts", "Main"),
        alt_counts    = lambda wc: pp(wc.sample, wc.plate, "read_counts", "Alt"),
        lambda_counts = (lambda wc: pp(wc.sample, wc.plate, "read_counts", "Lambda")) if UseSpikes else [],
        puc19_counts  = (lambda wc: pp(wc.sample, wc.plate, "read_counts", "pUC19")) if UseSpikes else [],
        meth_summary  = lambda wc: pp(wc.sample, wc.plate, "meth_summary"),
    output:
        OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv"
    log:
        OutDir + "/Logs/WellStats/{sample}/{plate}_well_stats.log"
    params:
        mappings   = lambda wc, input: ",".join(input.mappings),
        spike_args = lambda wc, input: (
            f'--lambda-counts "{input.lambda_counts}" --puc19-counts "{input.puc19_counts}"'
        ) if UseSpikes else "",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/BuildWellStats.R" \
            --mappings     "{params.mappings}" \
            --main-counts  "{input.main_counts}" \
            --alt-counts   "{input.alt_counts}" \
            --meth-summary "{input.meth_summary}" \
            {params.spike_args} \
            --output       "{output}" >> {log} 2>&1
        """

# Distribution of per-well CpG coverage by well type, showing how far single-cell wells
# separate from empty ones.
rule plot_cpg_coverage:
    input:
        well_stats = OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv",
    output:
        OutDir + "/Plots/CpG_Coverage/{sample}_{plate}_boxplot_cpg_coverage.png",
    log:
        OutDir + "/Logs/Beds/Main/{sample}/{plate}_cpg_coverage.log"
    params:
        prefix = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/PlotCpGCoverage.R" \
            --well-stats "{input.well_stats}" \
            --prefix     "{params.prefix}" \
            --output     "{output}" >> {log} 2>&1
        """

# Mean per-well methylation in each cytosine context, pooled across the plate. CpG should be
# high; CHG and CHH near zero unless bisulfite conversion failed.
rule plot_methylation_rates:
    input:
        well_stats = OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv",
    output:
        OutDir + "/Plots/MethylationRates/{sample}_{plate}_boxplot_methylation_rates.png",
    log:
        OutDir + "/Logs/Beds/Main/{sample}/{plate}_methylation_rates.log"
    params:
        prefix = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/PlotMethylationRates.R" \
            --well-stats "{input.well_stats}" \
            --prefix     "{params.prefix}" \
            --output     "{output}" >> {log} 2>&1
        """


##############################################################################################
### MethScan

# The plate's per-well CpG beds (first-half output, located via pp()) are symlinked
# into a staging dir and converted to MethSCAn's internal format.
rule methscan_prepare:
    input:
        beds_done = lambda wc: pp(wc.sample, wc.plate, "beds_done")
    output:
        done = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared/.PrepareDone"
    log:
        OutDir + "/Logs/MethScan/{sample}/{plate}_prepare.log"
    params:
        beds_dir     = lambda wc: pp(wc.sample, wc.plate, "beds_dir"),
        staging_dir  = OutDir + "/Workup/MethScan/{sample}/{plate}/Staged",
        prepared_dir = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared",
        prefix       = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_MethScan.yaml"
    shell:
        r"""
        set -euo pipefail
        rm -rf "{params.staging_dir}" "{params.prepared_dir}"
        mkdir -p "{params.staging_dir}"
        for bed in {params.beds_dir}/*_MethCPG.bed.gz; do
            well=$(basename "$bed" _MethCPG.bed.gz)
            ln -s "$bed" "{params.staging_dir}/{params.prefix}_${{well}}.bed.gz"
        done
        methscan prepare \
            --input-format biscuit_short \
            {params.staging_dir}/*.bed.gz \
            "{params.prepared_dir}" \
            > {log} 2>&1
        touch {output.done}
        """

# A plate-wide smoothed methylation baseline is estimated, which the variance scan measures
# individual cells against.
rule methscan_smooth:
    input:
        OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared/.PrepareDone"
    output:
        OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared/.SmoothDone"
    log:
        OutDir + "/Logs/MethScan/{sample}/{plate}_smooth.log"
    params:
        prepared_dir = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared",
    conda:
        WorkingDir + "/Envs/SNMC_MethScan.yaml"
    shell:
        r"""
        set -euo pipefail
        methscan smooth "{params.prepared_dir}" > {log} 2>&1
        touch {output}
        """

# The genome is scanned for variably methylated regions -- windows where cells disagree more
# than the plate baseline predicts. These regions become the features for clustering.
rule methscan_scan:
    input:
        OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared/.SmoothDone"
    output:
        OutDir + "/Results/MethScan/{sample}/{plate}_VMRs.bed"
    log:
        OutDir + "/Logs/MethScan/{sample}/{plate}_scan.log"
    params:
        prepared_dir = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared",
    threads:
        8
    conda:
        WorkingDir + "/Envs/SNMC_MethScan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output}")"
        methscan scan --threads {threads} "{params.prepared_dir}" "{output}" > {log} 2>&1
        """

# Each cell's methylation at every called region is extracted into a cell-by-region matrix of
# shrunken residuals, the input to compute_pca.
rule methscan_matrix:
    input:
        vmrs  = OutDir + "/Results/MethScan/{sample}/{plate}_VMRs.bed",
        ready = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared/.SmoothDone",
    output:
        residuals = OutDir + "/Results/MethScan/{sample}/{plate}/mean_shrunken_residuals.csv.gz",
        meth      = OutDir + "/Results/MethScan/{sample}/{plate}/methylated_sites.csv.gz",
        total     = OutDir + "/Results/MethScan/{sample}/{plate}/total_sites.csv.gz",
        fractions = OutDir + "/Results/MethScan/{sample}/{plate}/methylation_fractions.csv.gz",
    log:
        OutDir + "/Logs/MethScan/{sample}/{plate}_matrix.log"
    params:
        prepared_dir = OutDir + "/Workup/MethScan/{sample}/{plate}/Prepared",
        matrix_dir   = OutDir + "/Results/MethScan/{sample}/{plate}",
    threads:
        8
    conda:
        WorkingDir + "/Envs/SNMC_MethScan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.matrix_dir}"
        methscan matrix --threads {threads} \
            "{input.vmrs}" "{params.prepared_dir}" "{params.matrix_dir}" \
            > {log} 2>&1
        """


##############################################################################################
### 100kb-bin methylation matrix -- alternative to MethSCAn's VMR-based residuals matrix,
### built directly from the per-well BED files (coverage-weighted mean beta per bin).

# A second cell-by-region matrix from the same beds, on fixed 100 kb bins rather than
# data-driven regions -- a check that structure is not an artefact of region calling. Feeds
# its own compute_pca / leiden_cluster branch in parallel with the MethSCAn one.
rule bin100k_matrix:
    input:
        beds_done = lambda wc: pp(wc.sample, wc.plate, "beds_done")
    output:
        matrix = OutDir + "/Results/Bin100k/{sample}/{plate}/mean_meth_100k.csv.gz"
    log:
        OutDir + "/Logs/Bin100k/{sample}/{plate}_matrix.log"
    params:
        beds_pattern = lambda wc: pp(wc.sample, wc.plate, "beds_dir") + "/*_MethCPG.bed.gz",
        prefix       = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output.matrix}")"
        Rscript "{RScripts}/Compute100kBinMatrix.R" \
            --beds-pattern "{params.beds_pattern}" \
            --prefix       "{params.prefix}" \
            --output       "{output.matrix}" >> {log} 2>&1
        """


##############################################################################################
### PCA -> Leiden clustering -> plots, run independently on the MethSCAn residuals matrix and
### the 100kb-bin matrix. One persisted PCA per (plate, matrix type) is shared by the
### clustering, the UMAP scatter and the clustering panel. All wells are included --
### Empty, single-cell and multi-cell.

# The cell-by-region matrix is imputed, scaled and reduced to its leading principal
# components. One persisted embedding per (plate, matrix type), shared by leiden_cluster,
# plot_leiden_umap and plot_clustering_panel so they cannot disagree about the space.
rule compute_pca:
    input:
        matrix = lambda wc: (
            OutDir + f"/Results/MethScan/{wc.sample}/{wc.plate}/mean_shrunken_residuals.csv.gz"
            if wc.matrix_type == "methscan" else
            OutDir + f"/Results/Bin100k/{wc.sample}/{wc.plate}/mean_meth_100k.csv.gz"
        ),
    output:
        pca = OutDir + "/Results/PCA/{sample}/{plate}/{matrix_type}_pca.csv.gz"
    log:
        OutDir + "/Logs/PCA/{sample}/{plate}_{matrix_type}_pca.log"
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output.pca}")"
        Rscript "{RScripts}/ComputePCA.R" \
            --matrix  "{input.matrix}" \
            --pca-out "{output.pca}" >> {log} 2>&1
        """

# Cells are linked into a shared-nearest-neighbour graph in that reduced space and
# partitioned into clusters. Labels only -- the plots are separate rules.
rule leiden_cluster:
    input:
        pca = OutDir + "/Results/PCA/{sample}/{plate}/{matrix_type}_pca.csv.gz"
    output:
        clusters = OutDir + "/Results/Leiden/{sample}/{plate}/{matrix_type}_clusters.csv"
    log:
        OutDir + "/Logs/Leiden/{sample}/{plate}_{matrix_type}_cluster.log"
    params:
        prefix = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output.clusters}")"
        Rscript "{RScripts}/LeidenCluster.R" \
            --pca          "{input.pca}" \
            --prefix       "{params.prefix}" \
            --clusters-out "{output.clusters}" >> {log} 2>&1
        """

# Two-dimensional view of the same space coloured by cluster, with each cluster's empty-well
# content in the legend.
rule plot_leiden_umap:
    input:
        pca        = OutDir + "/Results/PCA/{sample}/{plate}/{matrix_type}_pca.csv.gz",
        clusters   = OutDir + "/Results/Leiden/{sample}/{plate}/{matrix_type}_clusters.csv",
        well_stats = OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv",
    output:
        OutDir + "/Plots/Leiden/{sample}/{plate}/{matrix_type}_umap_clusters.png"
    log:
        OutDir + "/Logs/Leiden/{sample}/{plate}_{matrix_type}_umap.log"
    params:
        prefix = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/PlotLeidenUMAP.R" \
            --pca        "{input.pca}" \
            --clusters   "{input.clusters}" \
            --well-stats "{input.well_stats}" \
            --prefix     "{params.prefix}" \
            --output     "{output}" >> {log} 2>&1
        """

# Per-cluster CpG coverage, which is what distinguishes a cluster of real cells from one of
# low-coverage empty wells.
rule plot_cpg_by_cluster:
    input:
        clusters   = OutDir + "/Results/Leiden/{sample}/{plate}/{matrix_type}_clusters.csv",
        well_stats = OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv",
    output:
        OutDir + "/Plots/Leiden/{sample}/{plate}/{matrix_type}_cpg_coverage_by_cluster.png"
    log:
        OutDir + "/Logs/Leiden/{sample}/{plate}_{matrix_type}_cpg_boxplot.log"
    params:
        prefix = "{sample}_{plate}",
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/PlotCpGByCluster.R" \
            --clusters   "{input.clusters}" \
            --well-stats "{input.well_stats}" \
            --prefix     "{params.prefix}" \
            --output     "{output}" >> {log} 2>&1
        """

# The same embedding coloured four ways -- cell type, cell count, coverage and mean CpG
# methylation -- to see which of them the structure actually tracks.
rule plot_clustering_panel:
    input:
        pca        = OutDir + "/Results/PCA/{sample}/{plate}/{matrix_type}_pca.csv.gz",
        well_stats = OutDir + "/Results/WellStats/{sample}_{plate}_well_stats.csv",
    output:
        OutDir + "/Plots/Clustering/{sample}/{plate}/{matrix_type}_panel.png"
    log:
        OutDir + "/Logs/Clustering/{sample}/{plate}_{matrix_type}_panel.log"
    conda:
        WorkingDir + "/Envs/SNMC_CGI.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        Rscript "{RScripts}/PlotClusteringPanel.R" \
            --pca        "{input.pca}" \
            --well-stats "{input.well_stats}" \
            --sample     "{wildcards.sample}" \
            --plate      "{wildcards.plate}" \
            --output     "{output}" >> {log} 2>&1
        """


##############################################################################################
### SpikeIn well heatmaps (Lambda + pUC19, abundance + methylation)

# Plate-layout view of spike-in abundance and methylation per well, read from the first-half
# pipeline's per-well spike-in VCFs.
rule plot_spikein_well_heatmap:
    input:
        vcf_done = lambda wc: pp(wc.sample, wc.plate, "spikein_done", wc.type)
    output:
        reads_heatmap = OutDir + "/Plots/SpikeIn/{sample}/{plate}/{type}/well_reads_heatmap.png",
        meth_heatmap  = OutDir + "/Plots/SpikeIn/{sample}/{plate}/{type}/well_methylation_heatmap.png",
    log:
        OutDir + "/Logs/SpikeIn/{sample}/{plate}_{type}_heatmap.log"
    params:
        vcf_dir     = lambda wc: pp(wc.sample, wc.plate, "spikein_dir", wc.type),
        cx           = lambda wc: "CG" if wc.type == "pUC19" else "all",
        invert_flag  = lambda wc: "--invert" if wc.type == "Lambda" else "",
        title        = lambda wc: f"{wc.sample} {wc.plate} {wc.type}",
        meth_label   = lambda wc: "Conversion" if wc.type == "Lambda" else "Methylation",
        min_cov      = 0,
    wildcard_constraints:
        type = "Lambda|pUC19"
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output.reads_heatmap}")"
        python3 "{PythonScripts}/PlotSpikeInWells.py" \
            --vcf-dir "{params.vcf_dir}" \
            --output  "{output.reads_heatmap}" \
            --metric  coverage \
            --cx      all \
            --min-cov {params.min_cov} \
            --title   "{params.title} — Coverage" >> {log} 2>&1
        python3 "{PythonScripts}/PlotSpikeInWells.py" \
            --vcf-dir "{params.vcf_dir}" \
            --output  "{output.meth_heatmap}" \
            --metric  methylation \
            --cx      {params.cx} \
            --min-cov {params.min_cov} \
            {params.invert_flag} \
            --title   "{params.title} — {params.meth_label}" >> {log} 2>&1
        """


##############################################################################################
### SpikeIn QC plots (per sample/plate)

# Bisulfite conversion efficiency (lambda) against false-negative rate (pUC19) for the plate:
# the primary chemistry QC for the experiment.
rule spike_qc_plots:
    input:
        lambda_done = lambda wc: pp(wc.sample, wc.plate, "spikein_done", "Lambda"),
        puc19_done  = lambda wc: pp(wc.sample, wc.plate, "spikein_done", "pUC19"),
    output:
        summary = OutDir + "/Workup/SpikeIn/{sample}_{plate}_spike_qc_summary.tsv",
        p1      = OutDir + "/Plots/SpikeIn/{sample}/{plate}/Combined/lambda_conversion_all.png",
        p2      = OutDir + "/Plots/SpikeIn/{sample}/{plate}/Combined/puc19_methylation_cpg.png",
        p3      = OutDir + "/Plots/SpikeIn/{sample}/{plate}/Combined/lambdaALL_vs_puc19CPG.png",
    log:
        OutDir + "/Logs/SpikeIn/{sample}/{plate}_qc_plots.log"
    params:
        lambda_vcf_dir = lambda wc: pp(wc.sample, wc.plate, "spikein_dir", "Lambda"),
        puc19_vcf_dir  = lambda wc: pp(wc.sample, wc.plate, "spikein_dir", "pUC19"),
        outdir         = OutDir + "/Plots/SpikeIn/{sample}/{plate}/Combined",
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname "{output.summary}")"
        python3 "{PythonScripts}/PlotSpikeIn.py" \
            --lambda-dir  "{params.lambda_vcf_dir}" \
            --puc19-dir   "{params.puc19_vcf_dir}" \
            --outdir      "{params.outdir}" \
            --sample      "{wildcards.sample}" \
            --plate       "{wildcards.plate}" \
            --summary-out "{output.summary}" >> {log} 2>&1
        """

