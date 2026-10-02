#!/usr/bin/env python3
"""
Per-cell mean methylation at the marker genes, one violin per gene and cell type.

Two region definitions, selected with --region:

  body      the whole gene, [start, end)
  promoter  TSS - --promoter-upstream to TSS + --promoter-downstream, asymmetric because
            a promoter is not centred on the TSS: the regulatory sequence sits upstream
            and only the first few hundred bases of the transcript belong to it. The
            default 2 kb / 500 bp is the common convention.

The original figure used a symmetric TSS +/- 2 kb, recovered by testing candidate windows
against its surviving per-cell tables (r = 0.94 for the promoter betas, r = 0.91 for the
body), so numbers here will not match it exactly.

A cell's value is the plain mean of its own CpG betas in the region. At single-cell
depth that is a handful of CpGs, so a large share of the values are exactly 0 or 1 --
the log reports how many. --min-cpgs drops the thinnest cells; it defaults to off, which
is what the original did.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import sys
sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import FONT_STACK, apply_style, save_figure
apply_style()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLORS = {"HealthyKeratinocytes": "#2a78d6", "MerkelCarcinoma": "#eb6834"}
DISPLAY = {"HealthyKeratinocytes": "Keratinocytes", "MerkelCarcinoma": "Merkel Cell Carcinoma"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--observations", required=True)
    p.add_argument("--bodies", required=True)
    p.add_argument("--region", required=True, choices=["body", "promoter"])
    p.add_argument("--promoter-upstream", type=int, default=2000)
    p.add_argument("--promoter-downstream", type=int, default=500)
    p.add_argument("--left-genes", required=True, help="first group, e.g. the MCC markers")
    p.add_argument("--right-genes", required=True, help="second group, e.g. the keratinocyte markers")
    p.add_argument("--min-cpgs", type=int, default=1)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    genes = args.left_genes.split(",") + args.right_genes.split(",")
    d = pd.read_parquet(args.observations)
    bodies = pd.read_csv(args.bodies, sep="\t", header=None,
                         names=["chrom", "start", "end", "gene", "strand"]).set_index("gene")
    d["len"] = d["gene"].map(bodies["end"] - bodies["start"])

    if args.region == "body":
        d = d[(d["pos"] >= 0) & (d["pos"] < d["len"])]
    else:
        d = d[(d["pos"] >= -args.promoter_upstream) & (d["pos"] <= args.promoter_downstream)]

    per_cell = (d.groupby(["gene", "cell_type", "plate", "well"])["beta"]
                 .agg(mean_beta="mean", n_covered="size").reset_index())
    per_cell = per_cell[per_cell["n_covered"] >= args.min_cpgs]
    sat = ((per_cell.mean_beta == 0) | (per_cell.mean_beta == 1)).mean()
    print(f"{args.region}: {len(per_cell)} (gene, cell) values, median {int(per_cell.n_covered.median())} "
          f"CpGs per cell, {100 * sat:.1f}% exactly 0 or 1", flush=True)
    per_cell.to_csv(args.output.replace(".png", "_percell.csv"), index=False)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    positions, data, colors, ticks = [], [], [], []
    pos = 0
    for gene in genes:
        group = []
        for ct in ["HealthyKeratinocytes", "MerkelCarcinoma"]:
            v = per_cell.loc[(per_cell.gene == gene) & (per_cell.cell_type == ct), "mean_beta"].to_numpy()
            data.append(v); colors.append(COLORS[ct]); positions.append(pos); group.append(pos)
            ax.text(pos, 1.08, str(len(v)), ha="center", va="bottom", fontsize=9, color=COLORS[ct])
            pos += 1
        ticks.append(np.mean(group))
        pos += 1

    vp = ax.violinplot(data, positions=positions, widths=0.8, showmedians=True, showextrema=False)
    for body, c in zip(vp["bodies"], colors):
        body.set_facecolor(c); body.set_alpha(0.6); body.set_edgecolor("none")
    vp["cmedians"].set_color("black")

    # divider sits between the two marker groups, derived from --left-genes
    n_left = len(args.left_genes.split(","))
    ax.axvline((positions[2 * n_left - 1] + positions[2 * n_left]) / 2,
               color="grey", linestyle=":", linewidth=1)
    ax.set_xticks(ticks); ax.set_xticklabels(genes)
    ax.set_ylabel("Mean methylation (per cell)")
    ax.set_ylim(-0.05, 1.15); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    label = "Gene Body" if args.region == "body" else "Promoter"
    ax.set_title(f"{label} Methylation at Marker Genes")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.6) for c in COLORS.values()]
    ax.legend(handles, [DISPLAY[c] for c in COLORS], loc="upper center",
              bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    save_figure(fig, args.output, dpi=150)
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
