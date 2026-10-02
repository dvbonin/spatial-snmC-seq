#!/usr/bin/env python3
"""
GENCODE gene records -> two region files the metagene profiles are built on:

  gene_bodies.bed   chrom, start, end, gene, strand   (0-based, half-open)
  gene_tss.bed      chrom, tss, gene, strand          (0-based TSS position)

A gene_name can appear on several GTF records (PAR copies, scaffold copies, plain
duplicates). Expression is keyed on the symbol, so exactly one locus per symbol has to
be chosen; the longest span wins, with GTF order breaking exact ties (stable sort), which
is deterministic and picks the real locus rather than whatever the file lists first. Only chr1-22, X, Y are kept -- the
sequencing data has no usable coverage on the scaffolds anyway.

All gene_types are kept, not just protein_coding, so that the expression tables (which
carry lncRNAs and pseudogenes) can be matched on symbol without silently dropping rows.
"""
import argparse
import gzip
import re

import pandas as pd

PRIMARY = ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--bodies", required=True)
    p.add_argument("--tss", required=True)
    args = p.parse_args()

    rows = []
    with gzip.open(args.gtf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            fl = line.rstrip("\n").split("\t")
            if fl[2] != "gene" or fl[0] not in PRIMARY:
                continue
            rows.append((fl[0], int(fl[3]) - 1, int(fl[4]), fl[6],
                         re.search(r'gene_name "([^"]+)"', fl[8]).group(1)))

    g = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"])
    n_records = len(g)
    g = (g.assign(length=g["end"] - g["start"])
          .sort_values(["gene", "length"], ascending=[True, False], kind="mergesort")
          .drop_duplicates("gene", keep="first"))

    g[["chrom", "start", "end", "gene", "strand"]].to_csv(args.bodies, sep="\t", header=False, index=False)
    print(f"Wrote {args.bodies}: {len(g)} genes from {n_records} primary-chromosome records")

    # 0-based TSS: the body start on +, the last base of the body on -
    g["tss"] = g["start"].where(g["strand"] == "+", g["end"] - 1)
    g[["chrom", "tss", "gene", "strand"]].to_csv(args.tss, sep="\t", header=False, index=False)
    print(f"Wrote {args.tss}: {len(g)} genes")


if __name__ == "__main__":
    main()
