#!/usr/bin/env python3
"""
WellMethStats.py

Per-well, per-context methylation summary.

CG comes from the merged-CpG bed (one row per dyad); CHG and CHH come from the
pileup VCF. So n_sites for CG -- which becomes CpGsCovered downstream -- counts
CpG dyads, not single cytosines.

One `biscuit vcf2bed -t c` pass covers CG, CHG and CHH together, replacing the
two passes (-t cg, -t ch) the old ProcessCells.py made plus the third pass the
Analysis pipeline's ExtractMethStats.py made independently.

The reported methylation is `mean_beta`: each covered cytosine's own
methylation fraction M/(M+U) is computed first, then averaged over sites. That
answers "what share of covered sites is methylated", giving every site equal
weight. It is deliberately NOT the read-pooled ratio sum(M)/sum(M+U), which
weights a site by its coverage -- in this data sites are covered ~2x on average
(diploid alleles plus overlapping mates), so the two differ by ~0.03 and the
gap tracks coverage.

n_calls is kept alongside n_sites so per-well mean coverage (n_calls/n_sites)
is available downstream.
"""

import argparse
import gzip
import os
import subprocess

CONTEXTS = ("CG", "CHG", "CHH")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--well",       required=True)
    p.add_argument("--vcf",        required=True)
    p.add_argument("--bed",        required=True,
                   help="Merged-CpG bed from WellCpGBed.py -- the CG row is counted from "
                        "this rather than from the VCF, so a CpG always means a dyad and "
                        "CpGsCovered cannot drift from the bed everything else reads.")
    p.add_argument("--stats-out",  required=True)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.stats_out), exist_ok=True)

    totals = {c: {"n_sites": 0, "n_calls": 0, "n_methylated": 0, "sum_beta": 0.0}
              for c in CONTEXTS}

    # CG: one row per dyad, beta in column 4 and coverage in column 5.
    cg = totals["CG"]
    with gzip.open(args.bed, "rt") as bed:
        for line in bed:
            f = line.rstrip("\n").split("\t")
            beta, cov = float(f[3]), int(f[4])
            cg["n_sites"] += 1
            cg["n_calls"] += cov
            cg["n_methylated"] += round(beta * cov)
            cg["sum_beta"] += beta

    # CHG/CHH from the VCF; CG rows are skipped, having just been counted above.
    proc = subprocess.Popen(
        ["biscuit", "vcf2bed", "-t", "c", "-e", "-c", "-k", "1", args.vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        ctx = f[4]
        if ctx == "CG":
            continue
        t = totals.get(ctx)
        if t is None:
            continue
        m, u = int(f[8]), int(f[9])
        cov = m + u
        if cov == 0:
            continue
        t["n_sites"] += 1
        t["n_calls"] += cov
        t["n_methylated"] += m
        t["sum_beta"] += m / cov

    stderr = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"vcf2bed failed for {args.vcf}:\nSTDERR:\n{stderr}")

    with open(args.stats_out, "w") as out:
        out.write("context\tn_sites\tn_calls\tn_methylated\tmean_beta\n")
        for c in CONTEXTS:
            t = totals[c]
            mean_beta = t["sum_beta"] / t["n_sites"] if t["n_sites"] else float("nan")
            out.write(f"{c}\t{t['n_sites']}\t{t['n_calls']}\t{t['n_methylated']}\t{mean_beta}\n")

    print(f"Done: {args.well}")


if __name__ == "__main__":
    main()
