#!/usr/bin/env python3
"""
Bar chart of the share of reads landing on each primary chromosome, summed
across one (sample, plate)'s chunks. chrM gets its own bar -- it carries ~10-20x
more reads than all other non-primary contigs combined, so lumping it in would
hide both signals. Everything else non-primary (EBV, decoy, unplaced and
unlocalized scaffolds) goes into a single "Other" bar; otherwise ~80 tiny bars
swamp the plot without being individually informative.

Input is the per-chunk contig_counts.tsv written by GroupReadsToCells.py, i.e.
counted off the deduped chunk BAM that feeds the well split -- post-MAPQ-filter,
post-barcode-correction, post-dedup. These numbers therefore reconcile with the
per-well stream heatmaps, and are NOT comparable to `samtools idxstats` on the
raw alignment.
"""
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRIMARY_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
MITO = "chrM"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contig-counts", nargs="+", required=True,
                    help="One or more contig_counts.tsv files (contig, records)")
    ap.add_argument("--title", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    counts = defaultdict(int)
    for path in args.contig_counts:
        with open(path) as f:
            next(f)  # header
            for line in f:
                contig, n = line.rstrip("\n").split("\t")
                counts[contig] += int(n)

    total = sum(counts.values())
    other = sum(n for c, n in counts.items()
                if c not in PRIMARY_ORDER and c != MITO)

    labels = PRIMARY_ORDER + [MITO, "Other"]
    shares = ([100 * counts.get(c, 0) / total for c in PRIMARY_ORDER]
              + [100 * counts.get(MITO, 0) / total, 100 * other / total])
    colors = ["#4393c3"] * len(PRIMARY_ORDER) + ["#d95f02", "#999999"]

    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.bar(labels, shares, color=colors)
    ax.set_ylabel("Share of split-input reads (%)")
    ax.set_title(f"Read alignment share by contig{' -- ' + args.title if args.title else ''}")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output} ({total} reads; "
          f"{100 * counts.get(MITO, 0) / total:.3g}% chrM, {100 * other / total:.3g}% other)",
          flush=True)


if __name__ == "__main__":
    main()
