#!/usr/bin/env python3
"""
Every CpG observation of every QC-passed single cell across the marker gene loci, from
--upstream bp before the TSS through to the TES, as one long table.

The range covers exactly what the two violin figures need: the promoter window is
TSS +/- 2 kb and the body is the whole gene, so 2 kb upstream through the TES suffices.

This is per-cell data, so it cannot come from a pseudobulk -- it is read straight out of
the per-well beds. The loci are tiny (five genes, ~58 kb in total), so the beds are
queried by region with tabix rather than read whole: five regions per well instead of
50 million rows per well.

Coordinates are oriented per gene: pos is signed distance from the TSS, running from
-upstream through 0 at the TSS to the gene length at the TES, on both strands. beta is
the well's own value at that CpG, so a cell with one read there contributes a 0 or a 1.
"""
import os
import argparse
import io
import subprocess

import pandas as pd

CELL_TYPES = ["HealthyKeratinocytes", "MerkelCarcinoma"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--well-stats", required=True,
                   help="Comma-separated beds_dir=well_stats_csv pairs, one per plate")
    p.add_argument("--bodies", required=True, help="gene_bodies.bed: chrom,start,end,gene,strand")
    p.add_argument("--genes", required=True, help="Comma-separated gene symbols")
    p.add_argument("--upstream", type=int, default=2000)
    p.add_argument("--tabix", default=os.environ.get("SNMC_TABIX", "tabix"))
    p.add_argument("--output", required=True)
    args = p.parse_args()

    want = args.genes.split(",")
    bodies = pd.read_csv(args.bodies, sep="\t", header=None,
                         names=["chrom", "start", "end", "gene", "strand"])
    g = bodies[bodies["gene"].isin(want)].set_index("gene").loc[want].reset_index()

    # TSS-upstream .. TES, in genomic orientation
    g["lo"] = (g["start"] - args.upstream).where(g["strand"] == "+", g["start"])
    g["hi"] = g["end"].where(g["strand"] == "+", g["end"] + args.upstream)
    g["tss"] = g["start"].where(g["strand"] == "+", g["end"] - 1)
    regions = [f"{r.chrom}:{r.lo + 1}-{r.hi}" for r in g.itertuples()]
    for r in g.itertuples():
        print(f"  {r.gene:7s} {r.chrom}:{r.lo}-{r.hi}  strand {r.strand}  "
              f"body {r.end - r.start} bp", flush=True)

    rows = []
    for spec in args.well_stats.split(","):
        beds_dir, ws_path = spec.split("=")
        plate = beds_dir.rstrip("/").split("/")[-1]
        ws = pd.read_csv(ws_path)
        sel = ws[ws["well_class"].isin(CELL_TYPES) & (ws["passes_cpg_qc"] == True)]
        print(f"{ws_path.split('/')[-1]}: {len(sel)} QC-passed wells", flush=True)

        for w in sel.itertuples():
            out = subprocess.run([args.tabix, f"{beds_dir}/{w.WellPosition}_MethCPG.bed.gz",
                                  *regions], capture_output=True, text=True, check=True).stdout
            if not out:
                continue
            # Merged-CpG bed: chrom, start, end, beta, coverage, detail. Rows are CpG
            # dyads and vcf2bed emits CG only, so there is no context column to filter.
            d = pd.read_csv(io.StringIO(out), sep="\t", header=None,
                            usecols=[0, 1, 3, 4], names=["chrom", "pos", "beta", "cov"])
            for r in g.itertuples():
                m = d[(d["chrom"] == r.chrom) & (d["pos"] >= r.lo) & (d["pos"] < r.hi)]
                if m.empty:
                    continue
                rel = m["pos"] - r.tss if r.strand == "+" else r.tss - m["pos"]
                rows.append(pd.DataFrame({
                    "gene": r.gene, "cell_type": w.well_class, "plate": plate,
                    "well": w.WellPosition, "pos": rel.to_numpy(),
                    "beta": m["beta"].to_numpy(), "cov": m["cov"].to_numpy()}))

    d = pd.concat(rows, ignore_index=True)
    d.to_parquet(args.output, index=False)
    print(f"\nWrote {args.output}: {len(d)} (cell, CpG) observations", flush=True)
    print(d.groupby(["gene", "cell_type"]).agg(cells=("well", "nunique"),
                                               cpgs=("pos", "nunique"),
                                               obs=("beta", "size")).to_string(), flush=True)


if __name__ == "__main__":
    main()
