#!/usr/bin/env python3
"""
Per-CpG pseudobulk -> per (gene, bin) methylation sums, for one cell type and one
region type.

  --mode genebody   gene body scaled to 100 bins, 0 = TSS, 1 = TES
  --mode promoter   TSS +/- flank bp in fixed bins of --bin-width

A CpG is credited to EVERY gene whose region contains it. That matters: gene bodies
nest (31,825 of the 77,081 bodies sit entirely inside a longer one) and TSS windows
overlap heavily, so one CpG legitimately belongs to several genes' profiles. The
regions are therefore walked one at a time, each taking the slice of CpG positions it
spans, rather than each CpG being matched to a single region.

Output rows carry mean_beta, the plain unweighted mean of the CpG betas that fell in
that gene's bin, so nothing downstream is weighted by coverage: a CpG seen in one cell
counts the same as one seen in forty. n_cpg is reported alongside for reference only.
"""
import argparse

import numpy as np
import pandas as pd

PRIMARY = ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pseudobulk", required=True, help="chrom, start, end, beta, coverage")
    p.add_argument("--regions", required=True, help="bodies bed (5 col) or tss bed (4 col)")
    p.add_argument("--mode", required=True, choices=["genebody", "promoter"])
    p.add_argument("--n-bins", type=int, default=100, help="genebody only")
    p.add_argument("--flank", type=int, default=5000, help="promoter only")
    p.add_argument("--bin-width", type=int, default=40, help="promoter only")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    if args.mode == "genebody":
        r = pd.read_csv(args.regions, sep="\t", header=None,
                        names=["chrom", "start", "end", "gene", "strand"])
        n_bins = args.n_bins
    else:
        r = pd.read_csv(args.regions, sep="\t", header=None,
                        names=["chrom", "tss", "gene", "strand"])
        r["start"] = np.maximum(r["tss"] - args.flank, 0)
        r["end"] = r["tss"] + args.flank + 1
        n_bins = 2 * args.flank // args.bin_width

    cpg = pd.read_csv(args.pseudobulk, sep="\t", header=None, usecols=[0, 1, 3],
                      names=["chrom", "pos", "beta"])
    cpg = cpg[cpg["chrom"].isin(PRIMARY)]
    print(f"{len(cpg)} CpGs, {len(r)} regions", flush=True)

    num = np.zeros((len(r), n_bins))
    cnt = np.zeros((len(r), n_bins), dtype=np.int32)

    for chrom, cs in cpg.groupby("chrom", sort=False):
        rs = r[r["chrom"] == chrom]
        if len(rs) == 0:
            continue
        cs = cs.sort_values("pos")
        pos = cs["pos"].to_numpy()
        beta = cs["beta"].to_numpy()

        lo = np.searchsorted(pos, rs["start"].to_numpy(), "left")
        hi = np.searchsorted(pos, rs["end"].to_numpy(), "left")
        rows = rs.index.to_numpy()
        starts, ends = rs["start"].to_numpy(), rs["end"].to_numpy()
        plus = (rs["strand"] == "+").to_numpy()
        tss = rs["tss"].to_numpy() if args.mode == "promoter" else None

        for j in np.flatnonzero(hi > lo):
            sl = slice(lo[j], hi[j])
            p_ = pos[sl]
            if args.mode == "genebody":
                L = ends[j] - starts[j]
                frac = (p_ - starts[j]) / L if plus[j] else (ends[j] - 1 - p_) / L
                b = (frac * n_bins).astype(np.int64)
            else:
                rel = p_ - tss[j] if plus[j] else tss[j] - p_
                b = np.clip((rel + args.flank) // args.bin_width, 0, n_bins - 1).astype(np.int64)
            i = rows[j]
            num[i] += np.bincount(b, weights=beta[sl], minlength=n_bins)
            cnt[i] += np.bincount(b, minlength=n_bins).astype(np.int32)

    ri, bi = np.nonzero(cnt)
    if args.mode == "genebody":
        coord = (bi + 0.5) / n_bins
    else:
        coord = -args.flank + (bi + 0.5) * args.bin_width
    out = pd.DataFrame({"gene": r["gene"].to_numpy()[ri], "bin": coord,
                        "mean_beta": num[ri, bi] / cnt[ri, bi], "n_cpg": cnt[ri, bi]})
    out.to_parquet(args.output, index=False)
    print(f"Wrote {args.output}: {len(out)} (gene, bin) rows, {out.gene.nunique()} genes, "
          f"{int(out.n_cpg.sum())} (CpG, gene) contributions", flush=True)


if __name__ == "__main__":
    main()
