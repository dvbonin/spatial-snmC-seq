import os
import re
import glob as _glob
from collections import defaultdict


##############################################################################################
#------------------------------  First half (FullPlate barcoding) ---------------------------#
##############################################################################################
# Raw fastq -> deduped per-type BAMs -> split-by-well -> merged per (sample, plate, type).
# Produces a 384-well output format shared with the quadrant-barcoded variant of this
# pipeline -- both
# variants feed the shared Pipelines/B1_Analysis_Combined/Analysis.smk the same way. Use this
# variant for plates barcoded as one full 384-well set (single Mapping_384.csv covering all
# wells).

##############################################################################################
### Configs

### Directories
WorkingDir     = config["WorkingDir"]
OutDir         = config["OutDir"]
DatasetConfigs = config["ConfigDir"]
InDir          = config["InDir"]

# References
MainRef = config["RefIndex"][config["Assembly"]]

# Spike-in controls
_spikes      = config.get("SpikeIns", {})
Spike_UnMeth = _spikes.get("Unmethylated", None)
Spike_Meth   = _spikes.get("Methylated",   None)
SpikeInRef   = _spikes.get("CombinedRef",  None)
UseSpikes    = bool(Spike_UnMeth and Spike_Meth and SpikeInRef)

# Alignment streams that exist as their own deduped chunk BAM. "Alt" is NOT one of
# these -- it is carved out of the Main stream at the well split, so it must never
# reach TagBarcodes/dupsifter (whose input lambdas have no Alt branch and would
# silently fall through to the pUC19 BAM).
STREAM_TYPES = ["Main"] + (["Lambda", "pUC19"] if UseSpikes else [])

# Per-well output types -- what split_chunk_by_well / merge_wells produce.
WELL_TYPES = ["Main", "Alt"] + (["Lambda", "pUC19"] if UseSpikes else [])

# Contig names for BAM splitting
_LAMBDA_CONTIG = _spikes.get("LambdaContig", None)
_PUC19_CONTIG  = _spikes.get("pUC19Contig",  None)


##############################################################################################
### Paths

PythonScripts = WorkingDir + "/Scripts/Python"
BashScripts   = WorkingDir + "/Scripts/Bash"
RScripts      = WorkingDir + "/Scripts/R"


##############################################################################################
### Variables

# config["InFiles"] maps sample name -> list of plate numbers to process.
# Chunks are already split on disk (SampleName_S<n>_P<plate>_C<chunk>_R1/R2.fastq.gz);
# every chunk found for a requested (sample, plate) is picked up automatically.
_RUN_TO_R1  = {}
_SP_TO_RUNS = defaultdict(list)

for _sample, _plates in config["InFiles"].items():
    for _plate_num in _plates:
        _plate = f"P{_plate_num}"
        for _r1 in sorted(_glob.glob(f"{InDir}/{_sample}_S*_{_plate}_C*_R1.fastq.gz")):
            _chunk = os.path.basename(_r1)[: -len("_R1.fastq.gz")]
            _RUN_TO_R1[_chunk] = _r1
            _SP_TO_RUNS[(_sample, _plate)].append(_chunk)

CHUNKS        = sorted(_RUN_TO_R1)
SAMPLE_PLATES = sorted(_SP_TO_RUNS)

def _chunk_dones_for_sp_type(wildcards):
    return [
        OutDir + f"/Workup/Split/Chunks/{chunk}/{wildcards.type}/.done"
        for chunk in _SP_TO_RUNS[(wildcards.sample, wildcards.plate)]
    ]

# QC only needs to sanity-check that a plate sequenced and aligned reasonably
# -- one representative chunk per (sample, plate) is enough, not every chunk.
QC_CHUNKS = sorted(sorted(_chunks)[0] for _chunks in _SP_TO_RUNS.values())


wildcard_constraints:
    chunk  = r"[A-Za-z0-9][A-Za-z0-9-]*_S\d+_P\d+_C\d+",
    sample = r"[A-Za-z0-9][A-Za-z0-9-]*",
    plate  = r"P\d+",
    well   = r"[A-P]\d{2}",
    type   = r"Main|Alt|Lambda|pUC19",


##############################################################################################
### Outputs

OUT_deduped = [
    OutDir + f"/Workup/Tag/{chunk}_{t}_Deduped.bam"
    for chunk in CHUNKS for t in STREAM_TYPES
]

OUT_split = [
    OutDir + f"/Results/Split/{s}/{p}/{t}_wells.tsv"
    for s, p in SAMPLE_PLATES for t in WELL_TYPES
]

OUT_pileup_beds = [
    OutDir + f"/Workup/Beds/Main/{s}/{p}/.BedsDone"
    for s, p in SAMPLE_PLATES
]

OUT_pileup_spikein = (
    [OutDir + f"/Workup/VCF/{t}/{s}/{p}/.VCFDone"
     for s, p in SAMPLE_PLATES for t in ("Lambda", "pUC19")]
    if UseSpikes else []
)

OUT_read_count_plots = (
    [OutDir + f"/Plots/ReadTracking/Main/{s}_{p}_Main_read_count_tracking.png"
     for s, p in SAMPLE_PLATES]
    + ([OutDir + f"/Plots/ReadTracking/Lambda/{s}_{p}_Lambda_read_count_tracking.png"
        for s, p in SAMPLE_PLATES]
       + [OutDir + f"/Plots/ReadTracking/pUC19/{s}_{p}_pUC19_read_count_tracking.png"
          for s, p in SAMPLE_PLATES]
       if UseSpikes else [])
)

OUT_contig_plots = [
    OutDir + f"/Plots/ContigAlignment/{s}_{p}_contig_alignment_share.png"
    for s, p in SAMPLE_PLATES
]

OUT_welltype_plots = [
    OutDir + f"/Plots/WellTypeCounts/{s}_{p}_{t}_well_counts.png"
    for s, p in SAMPLE_PLATES for t in WELL_TYPES
]

OUT_qc = (
    expand(OutDir + "/Results/FastQC/Raw/{chunk}_{read}_fastqc.html",
           chunk=QC_CHUNKS, read=["R1", "R2"])
    + expand(OutDir + "/Results/Trim_Fastp/{chunk}.html", chunk=QC_CHUNKS)
    + expand(OutDir + "/Results/SamtoolsStats/{chunk}.stats.txt", chunk=QC_CHUNKS)
    + [OutDir + "/Results/MultiQC/multiqc_report.html"]
)


##############################################################################################
### All rule

rule all:
    input:
        OUT_deduped
        + OUT_split
        + OUT_pileup_beds
        + OUT_pileup_spikein
        + OUT_read_count_plots
        + OUT_contig_plots
        + OUT_welltype_plots
        + OUT_qc


##############################################################################################
#--------------------------------------  Run Pipeline  --------------------------------------#
##############################################################################################

# Raw paired FASTQ chunks are adapter- and quality-trimmed and their poly-G/poly-X tails
# removed, so that nothing downstream has to reason about sequencing artefacts.
rule fastp_trim:
    input:
        r1 = lambda wc: _RUN_TO_R1[wc.chunk],
        r2 = lambda wc: _RUN_TO_R1[wc.chunk].replace("_R1.fastq.gz", "_R2.fastq.gz")
    output:
        r1   = OutDir + "/Workup/Trim_Fastp/{chunk}_R1_Trimmed.fastq.gz",
        r2   = OutDir + "/Workup/Trim_Fastp/{chunk}_R2_Trimmed.fastq.gz",
        html = OutDir + "/Results/Trim_Fastp/{chunk}.html",
        json = OutDir + "/Results/Trim_Fastp/{chunk}.json"
    threads: 8
    log:
        OutDir + "/Logs/Trim_Fastp/{chunk}.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe \
          --trim_poly_g \
          --trim_poly_x \
          --cut_tail \
          --cut_window_size 4 \
          --cut_mean_quality 15 \
          --length_required 20 \
          --n_base_limit 5 \
          --thread {threads} \
          --html {output.html} \
          --json {output.json} \
          > {log} 2>&1
        """

# Trimmed reads are aligned to the bisulfite-converted human genome. Confidently placed reads
# (MAPQ>=10) form the main-genome BAM; everything unmapped or ambiguously placed is handed
# back as FASTQ so it can be tried against the spike-in references instead.
rule biscuit_align:
    input:
        r1 = OutDir + "/Workup/Trim_Fastp/{chunk}_R1_Trimmed.fastq.gz",
        r2 = OutDir + "/Workup/Trim_Fastp/{chunk}_R2_Trimmed.fastq.gz",
    output:
        bam      = OutDir + "/Workup/Align_Biscuit/{chunk}_hg38.bam",
        r1       = OutDir + "/Workup/Align_Biscuit/{chunk}_leftover_R1.fastq.gz",
        r2       = OutDir + "/Workup/Align_Biscuit/{chunk}_leftover_R2.fastq.gz",
    threads:
        24
    resources:
        mem_mb = 200000
    log:
        OutDir + "/Logs/Align_Biscuit/{chunk}.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail

        tmpdir="{resources.tmpdir}/biscuit_sort_tmp/{wildcards.chunk}.${{SLURM_JOB_ID:-$$}}"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        biscuit align \
          -@ 16 \
          -R "@RG\\tID:{wildcards.chunk}\\tSM:{wildcards.chunk}\\tPL:ILLUMINA" \
          {MainRef} {input.r1} {input.r2} 2> {log} \
        | samtools view -u -F 2304 - \
        | samtools sort -@ 4 -m 2G -T "$tmpdir/{wildcards.chunk}" -O BAM \
            -o "$tmpdir/full.bam" >> {log} 2>&1
        # Primary mapped reads with MAPQ>=10 → hg38 BAM
        # Everything else (unmapped or low-MAPQ) → FASTQ for spike-in re-alignment
        samtools view -bF 4 -e 'mapq >= 10' \
            -U "$tmpdir/low.bam" -o {output.bam} "$tmpdir/full.bam" >> {log} 2>&1
        samtools sort -n -@ 4 "$tmpdir/low.bam" \
        | samtools fastq -@ 4 \
            -1 {output.r1} -2 {output.r2} \
            -0 /dev/null -s /dev/null >> {log} 2>&1
        """

# Reads that could not be placed on the human genome are realigned against the combined
# lambda/pUC19 reference, recovering the spike-in controls that were mixed into the library.
rule biscuit_align_spikein:
    input:
        r1 = OutDir + "/Workup/Align_Biscuit/{chunk}_leftover_R1.fastq.gz",
        r2 = OutDir + "/Workup/Align_Biscuit/{chunk}_leftover_R2.fastq.gz",
    output:
        bam = OutDir + "/Workup/Align_Biscuit/{chunk}_spikein.bam",
    params:
        ref = SpikeInRef
    threads:
        16
    resources:
        mem_mb = 50000
    log:
        OutDir + "/Logs/Align_Biscuit/{chunk}_spikein.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail

        tmpdir="{resources.tmpdir}/biscuit_spikein_tmp/{wildcards.chunk}.${{SLURM_JOB_ID:-$$}}"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        biscuit align \
          -@ 12 \
          -R "@RG\\tID:{wildcards.chunk}\\tSM:{wildcards.chunk}\\tPL:ILLUMINA" \
          {params.ref} {input.r1} {input.r2} 2> {log} \
        | samtools view -u -F 2308 - \
        | samtools sort -@ 4 -m 2G -T "$tmpdir/{wildcards.chunk}" -O BAM \
            -o {output.bam} >> {log} 2>&1
        """

# The spike-in alignment is separated by contig into its two control streams: lambda
# (unmethylated) and pUC19 (fully methylated).
rule split_spikein_bam:
    input:
        bam = OutDir + "/Workup/Align_Biscuit/{chunk}_spikein.bam"
    output:
        lambda_bam = OutDir + "/Workup/Align_Biscuit/{chunk}_lambda.bam",
        puc19_bam  = OutDir + "/Workup/Align_Biscuit/{chunk}_puc19.bam",
    params:
        lambda_contig = _LAMBDA_CONTIG,
        puc19_contig  = _PUC19_CONTIG,
    log:
        OutDir + "/Logs/Align_Biscuit/{chunk}_spikein_split.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools index {input.bam} 2>> {log}
        samtools view -bh {input.bam} {params.lambda_contig} \
            -o {output.lambda_bam} 2>> {log}
        samtools view -bh {input.bam} {params.puc19_contig} \
            -o {output.puc19_bam}  2>> {log}
        """

# Reads from the main-genome and spike-in alignments have their cell barcode checked and
# corrected against the plate whitelist, are grouped by read name, lose any read whose mate
# did not survive alignment filtering, and are finally deduplicated.
#
# The mate-pair filter is what defines "usable reads": dupsifter needs both mates and
# biscuit_align's MAPQ>=10 cut can drop one of them. Its count goes to PairFilter.txt.
# Pairing is by adjacent identical QNAME, so -F 0x900 keeps secondary and supplementary
# alignments out. The trailing `samtools view -b` is required -- dupsifter emits SAM.
rule dedup_reads:
    input:
        bam = lambda wc: (
            OutDir + f"/Workup/Align_Biscuit/{wc.chunk}_hg38.bam"   if wc.type == "Main"   else
            OutDir + f"/Workup/Align_Biscuit/{wc.chunk}_lambda.bam"  if wc.type == "Lambda" else
            OutDir + f"/Workup/Align_Biscuit/{wc.chunk}_puc19.bam"
        ),
        mapping = DatasetConfigs + "/Mapping_384.csv"
    output:
        bam         = OutDir + "/Workup/Tag/{chunk}_{type}_Deduped.bam",
        tag_stats   = OutDir + "/Results/Tag/{chunk}_{type}_CorrectionStats.txt",
        pair_stats  = OutDir + "/Results/Tag/{chunk}_{type}_PairFilter.txt",
    wildcard_constraints:
        type = "Main|Lambda|pUC19"
    params:
        ref = lambda wc: MainRef if wc.type == "Main" else SpikeInRef
    threads:
        8
    resources:
        mem_mb = lambda wc: 50000 if wc.type == "Main" else 20000
    log:
        OutDir + "/Logs/Dupsifter/{chunk}_{type}_Dedup.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")" "$(dirname "{output.tag_stats}")"

        tmpdir="{resources.tmpdir}/dedup_tmp/{wildcards.chunk}_{wildcards.type}.${{SLURM_JOB_ID:-$$}}"
        mkdir -p "$tmpdir"
        trap 'rm -rf "$tmpdir"' EXIT

        samtools view -h -F 0x900 {input.bam} \
        | awk 'BEGIN{{OFS="\t"}} /^@/{{print; next}}
               {{n=split($1,a,":"); $0=$0 "\tCB:Z:" a[n]; print}}' \
        | python3 {PythonScripts}/CorrectCBCTags.py {input.mapping} \
            --drop-ambiguous --drop-unmatched --max-distance 1 \
            --stats-out {output.tag_stats} \
        | samtools sort -n -@ {threads} -O sam -T "$tmpdir/ns" - \
        | awk -v stats={output.pair_stats} \
              '/^@/{{print; next}}
               prev==$1{{print prev_line; print; kept+=2; prev=""; next}}
               {{if(prev!="") dropped++; prev=$1; prev_line=$0}}
               END{{if(prev!="") dropped++;
                    print "kept_reads\t" kept+0 "\ndropped_singletons\t" dropped+0 > stats}}' \
        | dupsifter -r -B {params.ref} /dev/stdin 2> {log} \
        | samtools view -b -@ {threads} -o {output.bam}
        """


##############################################################################################
### Split by well: one job per chunk, then merge across chunks

# One sequencing chunk's main-genome reads are distributed into per-well BAMs by their
# corrected barcode, and at the same time separated by alignment contig: the primary
# chromosomes, which carry all downstream methylation analysis, and an Alt stream holding chrM
# plus the unplaced and decoy scaffolds.
rule split_chunk_by_well_main:
    input:
        bam = OutDir + "/Workup/Tag/{chunk}_Main_Deduped.bam"
    output:
        main_done   = temp(OutDir + "/Workup/Split/Chunks/{chunk}/Main/.done"),
        main_wells  = OutDir + "/Workup/Split/Chunks/{chunk}/Main/wells.tsv",
        alt_done    = temp(OutDir + "/Workup/Split/Chunks/{chunk}/Alt/.done"),
        alt_wells   = OutDir + "/Workup/Split/Chunks/{chunk}/Alt/wells.tsv",
        contigs     = OutDir + "/Workup/Split/Chunks/{chunk}/contig_counts.tsv",
    log:
        OutDir + "/Logs/Split/Chunks/{chunk}_Main.log"
    params:
        outdir  = OutDir + "/Workup/Split/Chunks/{chunk}",
        writers = 800,
    threads: 2
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        rm -rf "{params.outdir}/Main" "{params.outdir}/Alt"
        mkdir -p "{params.outdir}"
        python3 {PythonScripts}/GroupReadsToCells.py \
          --in                {input.bam} \
          --tag               XP \
          --outdir            "{params.outdir}" \
          --max-open          {params.writers} \
          --require-name-sorted \
          --contig-split \
          --contig-counts-out {output.contigs} \
          >> {log} 2>&1
        touch {output.main_done} {output.alt_done}
        """

# The same per-well distribution for the two spike-in streams, which need no contig
# separation.
rule split_chunk_by_well_spikein:
    input:
        bam = OutDir + "/Workup/Tag/{chunk}_{type}_Deduped.bam"
    output:
        done  = temp(OutDir + "/Workup/Split/Chunks/{chunk}/{type}/.done"),
        wells = OutDir + "/Workup/Split/Chunks/{chunk}/{type}/wells.tsv",
    wildcard_constraints:
        type = "Lambda|pUC19"
    log:
        OutDir + "/Logs/Split/Chunks/{chunk}_{type}.log"
    params:
        outdir  = OutDir + "/Workup/Split/Chunks/{chunk}/{type}",
        writers = 400,
    threads: 2
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        rm -rf "{params.outdir}"
        mkdir -p "{params.outdir}"
        python3 {PythonScripts}/GroupReadsToCells.py \
          --in         {input.bam} \
          --tag        XP \
          --outdir     "{params.outdir}" \
          --max-open   {params.writers} \
          --require-name-sorted \
          --wells-out  {output.wells} \
          >> {log} 2>&1
        touch {output.done}
        """


# A well's reads are scattered across every sequencing chunk of its plate; here they are
# concatenated back into one BAM per well and the per-well read totals are tallied.
checkpoint merge_wells:
    input:
        _chunk_dones_for_sp_type
    output:
        done   = OutDir + "/Workup/Merged/{sample}/{plate}/{type}/.SplitDone",
        wells  = OutDir + "/Results/Split/{sample}/{plate}/{type}_wells.tsv",
        counts = OutDir + "/Results/Split/{sample}/{plate}/{type}_read_counts.csv"
    log:
        OutDir + "/Logs/Split/{sample}/{plate}_{type}_merge.log"
    params:
        chunk_dirs = lambda wc: [
            OutDir + f"/Workup/Split/Chunks/{chunk}/{wc.type}"
            for chunk in _SP_TO_RUNS[(wc.sample, wc.plate)]
        ],
        outdir = OutDir + "/Workup/Merged/{sample}/{plate}/{type}"
    threads: 8
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        mkdir -p "$(dirname "{output.wells}")"
        python3 {PythonScripts}/MergeWellChunks.py \
          --chunk-dirs {params.chunk_dirs} \
          --outdir     "{params.outdir}" \
          --wells-out  {output.wells} \
          --counts-out {output.counts} \
          --threads    {threads} \
          >> {log} 2>&1
        touch {output.done}
        
        rm -rf {params.chunk_dirs}
        """


##############################################################################################
### Post-split: Pileup

# Main
def _main_wells_for(wildcards):
    # Dynamic well list for a (sample, plate), read from merge_wells' checkpoint
    # output -- lets pileup run as one job per well instead of one job per plate.
    ck = checkpoints.merge_wells.get(sample=wildcards.sample, plate=wildcards.plate, type="Main")
    wells = []
    with open(ck.output.wells) as f:
        next(f)
        for line in f:
            w = line.split("\t")[0].strip()
            if w:
                wells.append(w)
    return wells

# Every cytosine covered by a well's reads is called against the reference, giving methylated
# and unmethylated read counts at each genomic position.
#
# The pileup itself is the expensive, stable step; the bed and the per-context
# stats derived from its VCF are cheap and get redefined often. Keeping them in
# one rule meant any change to either forced a re-pileup of every well, so they
# are three rules over the same VCF.
rule pileup_well:
    input:
        bam = OutDir + "/Workup/Merged/{sample}/{plate}/Main/{well}.bam"
    output:
        vcf = OutDir + "/Workup/VCF/Main/{sample}/{plate}/{well}_VCF.gz",
        tbi = OutDir + "/Workup/VCF/Main/{sample}/{plate}/{well}_VCF.gz.tbi",
    log:
        OutDir + "/Logs/VCF/Main/{sample}/{plate}/{well}.log"
    params:
        tmp_dir = OutDir + "/Workup/Tmp/VCF/Main/{sample}/{plate}/{well}",
        ref     = MainRef
    threads:
        2
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        rm -f {output.vcf} {output.tbi}
        python3 "{PythonScripts}/PileupWell.py" \
            --well    "{wildcards.well}" \
            --bam     "{input.bam}" \
            --vcf-out "{output.vcf}" \
            --ref     "{params.ref}" \
            --tmp-dir "{params.tmp_dir}" \
            --threads {threads} >> {log} 2>&1
        """

# The per-cytosine calls are reduced to a CpG table, with the two cytosines of each CpG
# dyad merged into one site (biscuit mergecg) -- they carry the same biological mark, and
# keeping them apart split the coverage of every site in two. A "CpG" therefore means a
# dyad from here on, and rows span two bases.
#
# This is the sole entry point for every clustering result: Analysis.smk's methscan_prepare
# and bin100k_matrix both build their cell-by-region matrix from these beds, as do the
# merged-plate MethSCAn run and the pseudobulk script in Scripts/Figures. All of them read
# beta from column 4 and coverage from column 5 (MethSCAn: --input-format biscuit_short).
rule well_cpg_bed:
    input:
        vcf = OutDir + "/Workup/VCF/Main/{sample}/{plate}/{well}_VCF.gz",
        ref = MainRef
    output:
        bed = OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethCPG.bed.gz",
        tbi = OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethCPG.bed.gz.tbi",
    log:
        OutDir + "/Logs/Beds/Main/{sample}/{plate}/{well}_bed.log"
    params:
        tmp_dir = OutDir + "/Workup/Tmp/Beds/Main/{sample}/{plate}/{well}"
    threads:
        2
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        rm -f {output.bed} {output.tbi}
        python3 "{PythonScripts}/WellCpGBed.py" \
            --well    "{wildcards.well}" \
            --vcf     "{input.vcf}" \
            --bed-out "{output.bed}" \
            --ref     "{input.ref}" \
            --tmp-dir "{params.tmp_dir}" \
            --threads {threads} >> {log} 2>&1
        """

# The same calls are summarised per well into coverage and mean methylation for each cytosine
# context (CpG, CHG, CHH). CHG and CHH come from the VCF, which alone still holds them; CpG
# comes from the merged bed, so that CpGsCovered counts dyads and cannot drift from the bed
# every other consumer reads. Reaches the analysis pipeline via the plate summary below.
rule well_meth_stats:
    input:
        vcf = OutDir + "/Workup/VCF/Main/{sample}/{plate}/{well}_VCF.gz",
        bed = OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethCPG.bed.gz"
    output:
        stats = OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethStats.tsv"
    log:
        OutDir + "/Logs/Beds/Main/{sample}/{plate}/{well}_stats.log"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        rm -f {output.stats}
        python3 "{PythonScripts}/WellMethStats.py" \
            --well       "{wildcards.well}" \
            --vcf        "{input.vcf}" \
            --bed        "{input.bed}" \
            --stats-out  "{output.stats}" >> {log} 2>&1
        """

# A plate's per-well summaries are gathered into one table and the plate is marked complete.
# Both outputs cross into Analysis.smk: build_well_stats joins the summary into well_stats,
# and the sentinel is what gates methscan_prepare and bin100k_matrix on the plate being done.
rule collect_plate_meth_stats:
    input:
        wells_tsv = OutDir + "/Results/Split/{sample}/{plate}/Main_wells.tsv",
        beds      = lambda wc: expand(
            OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethCPG.bed.gz",
            sample=wc.sample, plate=wc.plate, well=_main_wells_for(wc),
        ),
        stats     = lambda wc: expand(
            OutDir + "/Workup/Beds/Main/{sample}/{plate}/{well}_MethStats.tsv",
            sample=wc.sample, plate=wc.plate, well=_main_wells_for(wc),
        ),
    output:
        sentinel = OutDir + "/Workup/Beds/Main/{sample}/{plate}/.BedsDone",
        summary  = OutDir + "/Results/MethStats/{sample}_{plate}_MethStats_Summary.csv",
    log:
        OutDir + "/Logs/Beds/Main/{sample}/{plate}_aggregate.log"
    params:
        beds_dir = OutDir + "/Workup/Beds/Main/{sample}/{plate}"
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.summary}")"
        python3 "{PythonScripts}/AggregateMethStats.py" \
            --wells-tsv   "{input.wells_tsv}" \
            --beds-dir    "{params.beds_dir}" \
            --summary-out "{output.summary}" \
            --sentinel    "{output.sentinel}" >> {log} 2>&1
        """

# Cytosine calling for the spike-in controls. Lambda reports how completely bisulfite
# conversion worked, pUC19 how often true methylation is missed.
#
# Lambda & pUC19
rule pileup_spikein_plate:
    input:
        wells_tsv = OutDir + "/Results/Split/{sample}/{plate}/{type}_wells.tsv"
    output:
        OutDir + "/Workup/VCF/{type}/{sample}/{plate}/.VCFDone"
    log:
        OutDir + "/Logs/VCF/{type}/{sample}/{plate}.log"
    params:
        split_dir = OutDir + "/Workup/Merged/{sample}/{plate}/{type}",
        vcf_dir   = OutDir + "/Workup/VCF/{type}/{sample}/{plate}",
        tmp_dir   = OutDir + "/Workup/Tmp/VCF/{type}/{sample}/{plate}",
        ref       = SpikeInRef
    wildcard_constraints:
        type = "Lambda|pUC19"
    threads:
        16
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output}")"
        python3 "{PythonScripts}/ProcessSpikeInCells.py" \
            --wells-tsv "{input.wells_tsv}" \
            --split-dir "{params.split_dir}" \
            --vcf-dir   "{params.vcf_dir}" \
            --spike-ref "{params.ref}" \
            --tmp-dir   "{params.tmp_dir}" \
            --threads   {threads} \
            --sentinel  "{output}" >> {log} 2>&1
        """


##############################################################################################
#----------------------------------  Read Count Tracking  -----------------------------------#
##############################################################################################

# Read survival is traced through every processing stage, from raw FASTQ to deduplicated
# alignments, so any loss can be attributed to the step that caused it.
rule collect_read_counts:
    input:
        done = OutDir + "/Workup/Merged/{sample}/{plate}/{type}/.SplitDone"
    output:
        OutDir + "/Results/ReadCounts/{sample}/{plate}/{type}_read_counts.csv"
    wildcard_constraints:
        type = "Main|Lambda|pUC19"
    params:
        runs        = lambda wc: ",".join(_SP_TO_RUNS[(wc.sample, wc.plate)]),
        r1_files    = lambda wc: ",".join(
                          _RUN_TO_R1[r] for r in _SP_TO_RUNS[(wc.sample, wc.plate)]),
        trimmed_dir = OutDir + "/Workup/Trim_Fastp",
        align_dir   = OutDir + "/Workup/Align_Biscuit",
        tag_dir       = OutDir + "/Workup/Tag",
        tag_stats_dir = OutDir + "/Results/Tag",
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        python3 "{PythonScripts}/CollectReadCounts.py" \
            --runs        "{params.runs}" \
            --r1-files    "{params.r1_files}" \
            --type        "{wildcards.type}" \
            --trimmed-dir "{params.trimmed_dir}" \
            --align-dir   "{params.align_dir}" \
            --tag-dir     "{params.tag_dir}" \
            --tag-stats-dir "{params.tag_stats_dir}" \
            --output      "{output}"
        """

# How the plate's reads distribute across the chromosomes, separating genuine nuclear signal
# from mitochondrial reads and from scaffold mismapping.
rule plot_contig_alignment_share:
    input:
        contigs = lambda wc: expand(
            OutDir + "/Workup/Split/Chunks/{chunk}/contig_counts.tsv",
            chunk=_SP_TO_RUNS[(wc.sample, wc.plate)],
        )
    output:
        OutDir + "/Plots/ContigAlignment/{sample}_{plate}_contig_alignment_share.png"
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        python3 "{PythonScripts}/PlotContigAlignmentShare.py" \
            --contig-counts {input.contigs} \
            --title   "{wildcards.sample} {wildcards.plate}" \
            --output  "{output}"
        """

# Plate-layout view of how many reads each well contributed to each stream, which makes
# spatial artefacts and failed wells immediately visible.
#
# One PNG per stream. All streams are still read together in a single job: the
# bracketed percentage on each cell is that stream's share of the well's total
# across all streams, so a per-stream job could not compute it.
rule plot_well_type_heatmaps:
    input:
        counts = lambda wc: expand(
            OutDir + "/Results/Split/{sample}/{plate}/{type}_read_counts.csv",
            sample=wc.sample, plate=wc.plate, type=WELL_TYPES,
        )
    output:
        expand(OutDir + "/Plots/WellTypeCounts/{{sample}}_{{plate}}_{type}_well_counts.png",
               type=WELL_TYPES)
    params:
        outdir = OutDir + "/Plots/WellTypeCounts",
        specs  = lambda wc: " ".join(
            f"{t}={OutDir}/Results/Split/{wc.sample}/{wc.plate}/{t}_read_counts.csv"
            for t in WELL_TYPES
        ),
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "{params.outdir}"
        python3 "{PythonScripts}/PlotWellTypeHeatmaps.py" \
            --counts {params.specs} \
            --outdir "{params.outdir}" \
            --prefix "{wildcards.sample}_{wildcards.plate}" \
            --title  "{wildcards.sample} {wildcards.plate}"
        """

# The stage-by-stage read survival for the main-genome stream, drawn as a funnel.
rule plot_main_read_count_tracking:
    input:
        summary = OutDir + "/Results/ReadCounts/{sample}/{plate}/Main_read_counts.csv"
    output:
        OutDir + "/Plots/ReadTracking/Main/{sample}_{plate}_Main_read_count_tracking.png"
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        python3 "{PythonScripts}/PlotReadCountTracking.py" \
            --summary "{input.summary}" \
            --title   "{wildcards.sample} {wildcards.plate} Main" \
            --output  "{output}"
        """

# The same survival funnel for one spike-in stream.
rule plot_spikein_read_count_tracking:
    input:
        csv = OutDir + "/Results/ReadCounts/{sample}/{plate}/{type}_read_counts.csv",
    output:
        OutDir + "/Plots/ReadTracking/{type}/{sample}_{plate}_{type}_read_count_tracking.png"
    wildcard_constraints:
        type = "Lambda|pUC19"
    conda:
        WorkingDir + "/Envs/SNMC_Analysis.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        python3 "{PythonScripts}/PlotReadCountTracking.py" \
            --summary    "{input.csv}" \
            --title      "{wildcards.sample} {wildcards.plate} {wildcards.type}" \
            --drop-first \
            --output     "{output}"
        """


##############################################################################################
#---------------------------------------  QC  -----------------------------------------------#
##############################################################################################

# Standard sequencing-quality report on the raw reads of one representative chunk, as a sanity
# check on the run itself.
rule fastqc_raw:
    input:
        r1 = lambda wc: _RUN_TO_R1[wc.chunk],
        r2 = lambda wc: _RUN_TO_R1[wc.chunk].replace("_R1.fastq.gz", "_R2.fastq.gz")
    output:
        html_r1 = OutDir + "/Results/FastQC/Raw/{chunk}_R1_fastqc.html",
        zip_r1  = OutDir + "/Results/FastQC/Raw/{chunk}_R1_fastqc.zip",
        html_r2 = OutDir + "/Results/FastQC/Raw/{chunk}_R2_fastqc.html",
        zip_r2  = OutDir + "/Results/FastQC/Raw/{chunk}_R2_fastqc.zip",
    params:
        outdir = OutDir + "/Results/FastQC/Raw"
    log:
        OutDir + "/Logs/FastQC/Raw/{chunk}.log"
    threads:
        2
    conda:
        WorkingDir + "/Envs/SNMC_QC.yaml"
    shell:
        r"""
        mkdir -p "{params.outdir}" "{resources.tmpdir}"
        export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir={resources.tmpdir}"
        fastqc -o "{params.outdir}" -t {threads} {input.r1} {input.r2} > {log} 2>&1
        """

# Alignment-level summary statistics for a representative chunk.
rule samtools_stats:
    input:
        bam = OutDir + "/Workup/Tag/{chunk}_Main_Deduped.bam"
    output:
        OutDir + "/Results/SamtoolsStats/{chunk}.stats.txt"
    log:
        OutDir + "/Logs/SamtoolsStats/{chunk}.log"
    threads:
        4
    conda:
        WorkingDir + "/Envs/SNMC_Bio.yaml"
    shell:
        r"""
        mkdir -p "$(dirname "{output}")"
        samtools stats -@ {threads} "{input.bam}" > "{output}" 2> {log}
        """

# The trimming, sequencing-quality and alignment reports are aggregated into one browsable
# page for the whole dataset.
rule multiqc:
    input:
        expand(OutDir + "/Results/Trim_Fastp/{chunk}.json", chunk=QC_CHUNKS),
        expand(OutDir + "/Results/FastQC/Raw/{chunk}_{read}_fastqc.zip",
               chunk=QC_CHUNKS, read=["R1", "R2"]),
        expand(OutDir + "/Results/SamtoolsStats/{chunk}.stats.txt", chunk=QC_CHUNKS),
    output:
        OutDir + "/Results/MultiQC/multiqc_report.html"
    params:
        search_dirs = " ".join([
            OutDir + "/Results/Trim_Fastp",
            OutDir + "/Results/FastQC/Raw",
            OutDir + "/Results/SamtoolsStats",
        ]),
        outdir = OutDir + "/Results/MultiQC"
    log:
        OutDir + "/Logs/MultiQC/multiqc.log"
    conda:
        WorkingDir + "/Envs/SNMC_QC.yaml"
    shell:
        r"""
        mkdir -p "{params.outdir}"
        multiqc {params.search_dirs} \
            -o "{params.outdir}" \
            --force \
            > {log} 2>&1
        """
