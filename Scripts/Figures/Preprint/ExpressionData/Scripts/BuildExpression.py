#!/usr/bin/env python3
"""
Two public bulk RNA-seq datasets -> one table of mean log-normalised expression per gene
and cell type, which is all the metagene plots need (they only rank genes into quartiles).

  cMCC  GSE223275, all 8 Merkel carcinoma samples of the raw count matrix
  bK    GSE107871, the four monolayer keratinocyte samples from healthy donors
        (series titles "KC NN 5".."KC NN 8", GSM2882398-2882401). The other 20 samples
        of that series are psoriasis patients or whole skin biopsies, so they are out.

Counts are turned into CPM, log1p'd, then averaged across samples -- log after
normalising, so a gene that is high in one sample and absent in another does not get its
mean dominated by the one large value.
"""
import argparse

import numpy as np
import pandas as pd

KZ_SAMPLES = ["GSM2882398_Sample9", "GSM2882399_Sample10",
              "GSM2882400_Sample11", "GSM2882401_Sample12"]


def mean_log_cpm(counts, subtype):
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    v = np.log1p(cpm).mean(axis=1)
    return pd.DataFrame({"gene": v.index, "subtype": subtype, "mean_lognorm": v.to_numpy()})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mcc-counts", required=True)
    p.add_argument("--kz-dir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    mcc = pd.read_csv(args.mcc_counts, sep="\t", index_col=0)
    mcc = mcc.loc[:, ~mcc.isna().all()]          # the header carries a trailing tab
    mcc = mcc.groupby(mcc.index).sum()

    kz = pd.DataFrame({s: pd.read_csv(f"{args.kz_dir}/{s}.txt.gz", sep=" ", quotechar='"')
                        .groupby("gene_symbol")["count"].sum() for s in KZ_SAMPLES})

    out = pd.concat([mean_log_cpm(mcc, "cMCC"), mean_log_cpm(kz, "bK")], ignore_index=True)
    out.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")
    print(f"  cMCC: {mcc.shape[1]} samples, {mcc.shape[0]} genes")
    print(f"  bK:   {kz.shape[1]} samples, {kz.shape[0]} genes")


if __name__ == "__main__":
    main()
