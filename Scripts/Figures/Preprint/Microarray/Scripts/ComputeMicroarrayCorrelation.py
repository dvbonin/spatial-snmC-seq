#!/usr/bin/env python3
"""
Per-well Pearson correlation of single-cell CpG methylation against the two EPIC
reference profiles, over the array's most differentially methylated sites.

For each well: pull only the top-diff positions out of the well's CpG VCF, convert
to per-site beta, and correlate against beta_mc and beta_kc. Restricting to sites
that separate the two references is what makes this work at single-cell coverage --
a genome-wide correlation is dominated by sites where both references agree and
carries almost no cell-type signal.

Both cytosines of each CpG dyad are used. An array probe names one cytosine, but
the mark is symmetric and either strand may have been the one sequenced, so the
tabix window is widened by a base to catch the partner and the two are summed
into one site. Querying the probe base alone halved the usable sites: on the
best-covered cell, 55,847 probes were hit directly and a further 55,772 only via
the partner.

The dyad partner is found from the reference base rather than by adjacency --
in a CGCG stretch positions 2 and 3 are both CpG cytosines but belong to
different dyads. Base "C" pairs with the following position, "G" with the
preceding one. That is the same arithmetic biscuit mergecg does, done here to
avoid 16 worker threads each loading the hg38 FASTA.

The site list comes from ../../Microarray, which builds it from the raw IDATs; the
VCFs come from the PP pipeline.
"""
import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

BIN = os.path.join(os.environ["SNMC_ENV_BIO"], "bin")
BISCUIT, TABIX = f"{BIN}/biscuit", f"{BIN}/tabix"

# Below this many shared sites the correlation is noise rather than a measurement.
MIN_SITES = 10


def well_correlation(vcf_dir, well, region_bed, array_df):
    vcf = f"{vcf_dir}/{well}_VCF.gz"
    tabix_out = subprocess.run([TABIX, "-h", "-R", region_bed, vcf],
                               capture_output=True, text=True).stdout
    if not tabix_out:
        return well, 0, np.nan, np.nan

    # -c gives M and U per cytosine, which is what lets the two strands be summed.
    bed = subprocess.run([BISCUIT, "vcf2bed", "-t", "cg", "-e", "-c", "-k", "1", "/dev/stdin"],
                         input=tabix_out, capture_output=True, text=True).stdout

    dyads = {}
    for line in bed.splitlines():
        f_ = line.split("\t")
        if len(f_) < 10:
            continue
        chrom, pos, base, m, u = f_[0], int(f_[1]), f_[3], f_[8], f_[9]
        if m == "." or u == ".":
            continue
        m, u = int(m), int(u)
        if m + u == 0:
            continue
        start = pos if base == "C" else pos - 1
        acc = dyads.setdefault((chrom, start), [0, 0])
        acc[0] += m
        acc[1] += u

    # A probe may name either base of the dyad, so offer the merged beta under both.
    rows = []
    for (chrom, start), (m, u) in dyads.items():
        beta = m / (m + u)
        rows.append((f"{chrom}:{start}", beta))
        rows.append((f"{chrom}:{start + 1}", beta))
    if len(dyads) < MIN_SITES:
        return well, len(dyads), np.nan, np.nan

    well_df = pd.DataFrame(rows, columns=["key", "beta"]).drop_duplicates(subset="key")
    merged = well_df.merge(array_df, on="key")
    if len(merged) < MIN_SITES:
        return well, len(merged), np.nan, np.nan

    return (well, len(merged),
            merged["beta"].corr(merged["beta_mc"]),
            merged["beta"].corr(merged["beta_kc"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf-dir", required=True, help="PP pipeline Workup/VCF/Main/<sample>/<plate>")
    ap.add_argument("--wells", required=True, help="comma-separated well IDs")
    ap.add_argument("--topdiff-bed", required=True,
                    help="TopDiffSites_array.bed.gz from ../../Microarray")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    array_df = pd.read_csv(args.topdiff_bed, sep="\t", header=None,
                           names=["chrom", "start", "end", "beta_mc", "beta_kc", "diff", "key"])
    print(f"{len(array_df):,} array sites", flush=True)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    # tabix -R needs a plain file, so the regions are written out alongside the result.
    region_bed = args.output.replace(".csv", "_regions.bed")
    # Widened by one base either side so the dyad partner of each probe is fetched too.
    regions = array_df[["chrom", "start", "end"]].copy()
    regions["start"] = (regions["start"] - 1).clip(lower=0)
    regions["end"] = regions["end"] + 1
    regions.to_csv(region_bed, sep="\t", header=False, index=False)

    wells = args.wells.split(",")
    print(f"{len(wells)} wells from {args.vcf_dir}", flush=True)

    rows = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = [ex.submit(well_correlation, args.vcf_dir, w, region_bed, array_df)
                   for w in wells]
        for i, fut in enumerate(as_completed(futures), 1):
            rows.append(fut.result())
            if i % 50 == 0:
                print(f"  {i}/{len(wells)} done", flush=True)

    out = pd.DataFrame(rows, columns=["well", "n_sites_covered", "corr_to_mcc", "corr_to_kz"])
    out.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({out['corr_to_mcc'].notna().sum()} wells with a correlation; "
          f"median {out['n_sites_covered'].median():.0f} shared sites)", flush=True)


if __name__ == "__main__":
    main()
