#!/usr/bin/env python3
"""
PileupWell.py

Per-well biscuit pileup: coordinate-sort the merged well BAM and call every
cytosine genome-wide into a bgzipped, tabixed VCF.

This is the expensive, stable half of what used to be ProcessCells.py. The
cheap, frequently-tweaked halves -- the per-CpG bed (WellCpGBed.py) and the
per-context stats (WellMethStats.py) -- are separate rules reading this VCF, so
changing how a bed or a statistic is defined no longer forces a re-pileup of
every well.
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
    p.add_argument("--bam",     required=True)
    p.add_argument("--vcf-out", required=True)
    p.add_argument("--ref",     required=True)
    p.add_argument("--tmp-dir", required=True)
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.vcf_out), exist_ok=True)
    tmp = args.tmp_dir
    os.makedirs(tmp, exist_ok=True)
    t = args.threads

    run(f"samtools sort -@ {t} -o {tmp}/coord.bam {args.bam}")
    run(f"samtools index -@ {t} {tmp}/coord.bam")
    run(f"biscuit pileup -@ {t} -o {tmp}/pileup.vcf {args.ref} {tmp}/coord.bam")
    run(f"bgzip -@ {t} -c {tmp}/pileup.vcf > {args.vcf_out}")
    run(f"tabix -p vcf {args.vcf_out}")

    shutil.rmtree(tmp, ignore_errors=True)
    print(f"Done: {args.well}")


if __name__ == "__main__":
    main()
