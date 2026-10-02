#!/usr/bin/env python3
"""
MergeWellChunks.py

Merge per-well BAMs from multiple chunk split directories into a single output
directory using samtools cat, and aggregate well read counts from per-chunk
wells.tsv files.
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--chunk-dirs", nargs="+", required=True,
                   help="Per-chunk split directories (each contains per-well BAMs and wells.tsv)")
    p.add_argument("--outdir", required=True,
                   help="Output directory for merged per-well BAMs")
    p.add_argument("--wells-out", required=True,
                   help="Output TSV: well<TAB>records_written")
    p.add_argument("--counts-out", required=True,
                   help="Output CSV: WellPosition;ReadCounts")
    p.add_argument("--threads", type=int, default=4,
                   help="Number of parallel samtools cat merges")
    return p.parse_args()


def merge_well(well, chunk_dirs, outdir):
    parts = [
        os.path.join(d, f"{well}.bam")
        for d in chunk_dirs
        if os.path.exists(os.path.join(d, f"{well}.bam"))
    ]
    out_bam = os.path.join(outdir, f"{well}.bam")
    subprocess.run(
        ["samtools", "cat", "-o", out_bam] + parts,
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Aggregate per-well counts from per-chunk wells.tsv files
    per_key = defaultdict(int)
    for d in args.chunk_dirs:
        tsv = os.path.join(d, "wells.tsv")
        if not os.path.exists(tsv):
            sys.stderr.write(f"# WARNING: no wells.tsv in {d}\n")
            continue
        with open(tsv) as f:
            next(f)  # skip header
            for line in f:
                well, count = line.rstrip("\n").split("\t")
                per_key[well] += int(count)

    all_wells = sorted(per_key)
    sys.stderr.write(f"# merging {len(all_wells)} wells from {len(args.chunk_dirs)} chunks\n")

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(merge_well, well, args.chunk_dirs, args.outdir): well
            for well in all_wells
        }
        for f in as_completed(futures):
            f.result()

    os.makedirs(os.path.dirname(args.wells_out), exist_ok=True)
    with open(args.wells_out, "w") as out:
        out.write("well\trecords_written\n")
        for k in all_wells:
            out.write(f"{k}\t{per_key[k]}\n")

    os.makedirs(os.path.dirname(args.counts_out), exist_ok=True)
    with open(args.counts_out, "w") as out:
        out.write("WellPosition;ReadCounts\n")
        for k in all_wells:
            out.write(f"{k};{per_key[k]}\n")

    sys.stderr.write(
        f"# merge complete\n"
        f"wells_merged\t{len(all_wells)}\n"
        f"total_records\t{sum(per_key.values())}\n"
    )


if __name__ == "__main__":
    main()
