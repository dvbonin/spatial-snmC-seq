#!/usr/bin/env python3
"""
WellCpGBed.py

Per-well CpG bed from the pileup VCF -- the file MethSCAn and the 100kb-bin
matrix read. Split out of the old ProcessCells.py so bed conventions can change
without re-running biscuit pileup.

Sites are CpG dyads, not single cytosines: `biscuit mergecg` sums the calls from
the C on the plus strand and the C on the minus strand (reference G, one base
further along), which are the same biological mark. The output therefore spans
two bases per row, and a CpG count means dyads everywhere downstream.

Note the flags: mergecg wants beta in column 4 and coverage in column 5, so
vcf2bed must run WITHOUT -e -- that flag inserts four context columns and pushes
them to 8 and 9. Dropping it also means consumers read columns 4/5, and MethSCAn
needs --input-format biscuit_short rather than biscuit.

No contig filtering happens here any more: the well BAM reaching pileup is
already primary-chromosome-only, because split_chunk_by_well_main routes
chrM/EBV/decoy/unplaced alignments into the separate Alt stream.
"""

import argparse
import os
import shutil
import subprocess


def run(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\nSTDERR:\n{result.stderr}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--well",    required=True)
    p.add_argument("--vcf",     required=True)
    p.add_argument("--bed-out", required=True)
    p.add_argument("--ref",     required=True, help="Reference FASTA, for mergecg")
    p.add_argument("--tmp-dir", required=True)
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.bed_out), exist_ok=True)
    tmp = args.tmp_dir
    os.makedirs(tmp, exist_ok=True)

    run(f"biscuit vcf2bed -k 1 {args.vcf} > {tmp}/percyt.bed")
    run(f"biscuit mergecg {args.ref} {tmp}/percyt.bed > {tmp}/meth.bed")
    run(f"bgzip -@ {args.threads} -c {tmp}/meth.bed > {args.bed_out}")
    run(f"tabix -p bed {args.bed_out}")

    shutil.rmtree(tmp, ignore_errors=True)
    print(f"Done: {args.well}")


if __name__ == "__main__":
    main()
