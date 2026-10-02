import argparse
import csv
import os
import re
import subprocess


def _count_fastq_pairs(fastq_gz):
    r = subprocess.run(f"zcat {fastq_gz} | wc -l", shell=True,
                       capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        return 0
    return int(r.stdout.strip()) // 4


def _count_bam_pairs(bam):
    r = subprocess.run(f"samtools view -c -F 2304 -f 64 {bam}",
                       shell=True, capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        return 0
    return int(r.stdout.strip())


_BAM_SUFFIX = {"Main": "_hg38.bam", "Lambda": "_lambda.bam", "pUC19": "_puc19.bam"}


def _assigned_pairs(stats_path):
    """Reads surviving barcode correction, from CorrectCBCTags.py's summary.

    Tag correction no longer has a BAM of its own -- it is streamed inside
    dedup_reads -- so the count comes from its stats file. The file counts SAM
    records (mates); every other column here is pairs, so halve it.
    """
    try:
        with open(stats_path) as f:
            for line in f:
                if "assigned XP tag" in line:
                    return int(re.search(r"([\d,]+)", line.split("assigned XP tag")[1])
                               .group(1).replace(",", "")) // 2
    except (OSError, AttributeError):
        pass
    return 0


def _pair_filtered_pairs(stats_path):
    """Pairs surviving the mate-pair filter, from dedup_reads' PairFilter.txt."""
    try:
        with open(stats_path) as f:
            for line in f:
                if line.startswith("kept_reads"):
                    return int(line.split("\t")[1]) // 2
    except (OSError, IndexError, ValueError):
        pass
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs",        required=True)
    ap.add_argument("--r1-files",    required=True)
    ap.add_argument("--type",        required=True, choices=["Main", "Lambda", "pUC19"])
    ap.add_argument("--trimmed-dir", required=True)
    ap.add_argument("--align-dir",   required=True)
    ap.add_argument("--tag-dir",     required=True)
    ap.add_argument("--tag-stats-dir", required=True,
                    help="Results/Tag -- CorrectionStats.txt / PairFilter.txt from dedup_reads")
    ap.add_argument("--output",      required=True)
    args = ap.parse_args()

    runs       = args.runs.split(",")
    r1_files   = args.r1_files.split(",")
    t          = args.type
    is_spikein = t in ("Lambda", "pUC19")
    bam_suffix = _BAM_SUFFIX[t]

    if is_spikein:
        fieldnames = ["run", "sample", "plate", "chunk",
                      "leftover", "aligned", "tagged", "paired", "deduped"]
    else:
        fieldnames = ["run", "sample", "plate", "chunk",
                      "raw", "trimmed", "aligned", "tagged", "paired", "deduped"]

    # FullPlate runs (e.g. "FirstPrep_S1_P3_C2"): chunk = numeric chunk index.
    # Quadrant runs (e.g. "Spatial5MC_S5_P2_Q1_C1"): chunk keeps the quadrant
    # prefix (e.g. "Q1_C1") so it stays unique within a plate.
    _RUN_RE_PLAIN    = re.compile(r"^([A-Za-z0-9-]+)_S\d+_(P\d+)_C(\d+)$")
    _RUN_RE_QUADRANT = re.compile(r"^([A-Za-z0-9-]+)_S\d+_(P\d+)_(Q[1-4]_C\d+)$")

    def _parse_run_name(run):
        m = _RUN_RE_QUADRANT.match(run)
        if m:
            return m.group(1), m.group(2), m.group(3)
        m = _RUN_RE_PLAIN.match(run)
        if m:
            return m.group(1), m.group(2), m.group(3)
        raise ValueError(f"Run name {run!r} does not match expected pattern")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for run, r1 in zip(runs, r1_files):
            sample, plate, chunk = _parse_run_name(run)

            if is_spikein:
                row = {
                    "run":      run,
                    "sample":   sample,
                    "plate":    plate,
                    "chunk":    chunk,
                    "leftover": _count_fastq_pairs(
                        os.path.join(args.align_dir, f"{run}_leftover_R1.fastq.gz")),
                    "aligned":  _count_bam_pairs(
                        os.path.join(args.align_dir, f"{run}{bam_suffix}")),
                    "tagged":   _assigned_pairs(
                        os.path.join(args.tag_stats_dir, f"{run}_{t}_CorrectionStats.txt")),
                    "paired":   _pair_filtered_pairs(
                        os.path.join(args.tag_stats_dir, f"{run}_{t}_PairFilter.txt")),
                    "deduped":  _count_bam_pairs(
                        os.path.join(args.tag_dir, f"{run}_{t}_Deduped.bam")),
                }
            else:
                row = {
                    "run":      run,
                    "sample":   sample,
                    "plate":    plate,
                    "chunk":    chunk,
                    "raw":      _count_fastq_pairs(r1),
                    "trimmed":  _count_fastq_pairs(
                        os.path.join(args.trimmed_dir, f"{run}_R1_Trimmed.fastq.gz")),
                    "aligned":  _count_bam_pairs(
                        os.path.join(args.align_dir, f"{run}{bam_suffix}")),
                    "tagged":   _assigned_pairs(
                        os.path.join(args.tag_stats_dir, f"{run}_{t}_CorrectionStats.txt")),
                    "paired":   _pair_filtered_pairs(
                        os.path.join(args.tag_stats_dir, f"{run}_{t}_PairFilter.txt")),
                    "deduped":  _count_bam_pairs(
                        os.path.join(args.tag_dir, f"{run}_{t}_Deduped.bam")),
                }
            writer.writerow(row)
            print(f"  {run}: done")


if __name__ == "__main__":
    main()
