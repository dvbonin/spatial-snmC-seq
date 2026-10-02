#!/usr/bin/env python3
"""
Per-well spike-in performance for both controls, as one table.

  Lambda  unmethylated control -> conversion is 1 - beta
  pUC19   CpG-methylated control -> methylation is beta

Each control is measured in the context set it was built for, which is not the same set
for the two. The asymmetry is the point, not an oversight:

  Lambda  ALL cytosines. It is unmethylated everywhere, so every C reports on conversion.
          Restricting it to CpGs raises beta from 0.066 to 0.238 -- apparent conversion
          93% -> 76% -- because human reads that cross-map onto lambda carry ~80% CpG
          methylation but ~1% CH methylation, so a CpG-only estimate measures conversion
          plus contamination.

  pUC19   CpG only. M.SssI methylates it at CpGs and nowhere else, so its CH sites are
          unmethylated by construction and including them dilutes the control: beta falls
          from 0.977 to 0.368 (in one well, 316 CpGs at beta 0.961 against 1031 CH sites
          at 0.132).

Beta is the coverage-weighted mean of the VCF's per-site BT over sites with CV > 0.

cov is the summed per-cytosine coverage and is what the box plots show. reads is the
Analysis pipeline's own spike-in read count, carried through for reference but not
plotted -- note the two are NOT interchangeable: on FirstPrep P3 they differ by a median
factor of 14 and range from 0.2x to 21x, because duplication scales them differently.

Every well of every plate is included -- all 384, empty and multi-cell alongside the
single cells -- because the point is to see the controls across the whole plate. well_class
and passes_cpg_qc are carried through so any subset can be taken later without recomputing.

There is no cache: reading 1536 small VCFs takes a couple of minutes, which is not worth
the staleness risk a cache carries.
"""
import argparse
import glob
import gzip
import os

import pandas as pd

SPIKES = ["Lambda", "pUC19"]
READ_COL = {"Lambda": "ReadCounts_Lambda", "pUC19": "ReadCounts_pUC19"}
# None means every context; see the module docstring for why the two differ.
CONTEXT = {"Lambda": None, "pUC19": "CX=CG"}


def spike_beta(vcf_path, context):
    """Coverage-weighted mean beta over the covered sites of one context, with count and depth."""
    total_cov, weighted_bt, n = 0, 0.0, 0
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or (context is not None and context not in fields[7]):
                continue
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            if "CV" not in fmt or "BT" not in fmt:
                continue
            cv = int(fmt["CV"])
            if cv > 0:
                total_cov += cv
                weighted_bt += float(fmt["BT"]) * cv
                n += 1
    return (weighted_bt / total_cov, n, total_cov) if total_cov else (None, 0, 0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--plates", nargs="+", required=True,
                   help="PLATE=lambda_dir,puc19_dir,well_stats_csv")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = []
    for spec in args.plates:
        plate, rest = spec.split("=", 1)
        lambda_dir, puc19_dir, ws_path = rest.split(",")
        ws = pd.read_csv(ws_path).set_index("WellPosition")
        for spike, d in zip(SPIKES, [lambda_dir, puc19_dir]):
            for vcf in sorted(glob.glob(os.path.join(d, "*_VCF.gz"))):
                well = os.path.basename(vcf).replace("_VCF.gz", "")
                beta, n, cov = spike_beta(vcf, CONTEXT[spike])
                if beta is None:
                    continue
                w = ws.loc[well]
                rows.append((plate, spike, well, w["well_class"], bool(w["passes_cpg_qc"]),
                             beta, n, cov, int(w[READ_COL[spike]])))
            print(f"  {plate} {spike}: {sum(r[0] == plate and r[1] == spike for r in rows)} wells",
                  flush=True)

    d = pd.DataFrame(rows, columns=["plate", "spike_type", "well", "well_class",
                                    "passes_cpg_qc", "beta", "n_sites", "cov", "reads"])
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    d.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output}: {len(d)} rows", flush=True)
    print(d.groupby(["plate", "spike_type"]).agg(
        wells=("well", "size"), median_beta=("beta", "median"),
        median_sites=("n_sites", "median"), median_cov=("cov", "median"),
        median_reads=("reads", "median")).round(4).to_string(),
        flush=True)


if __name__ == "__main__":
    main()
