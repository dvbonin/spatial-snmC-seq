#!/usr/bin/env python3
"""
Per-well spike-in coverage, one PNG per control, one box per plate, with every well
overlaid as a jittered point coloured red-to-green by how well that well performed:
conversion (1 - beta) for Lambda, methylation (beta) for pUC19. A well that worked is
green in both plots.

The y axis is summed per-site coverage over the control's genome, across whichever
context that control is measured in, which is not the same thing as a read count -- see
CollectSpikeIn.py.
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

CMAP = LinearSegmentedColormap.from_list("red_green", ["#D9776B", "#6FA287"])
PERFORMANCE = {"Lambda": ("Conversion", lambda b: 1 - b), "pUC19": ("Methylation", lambda b: b)}
# The context each control is measured in, mirroring CONTEXT in CollectSpikeIn.py.
CONTEXT_LABEL = {"Lambda": "All Cytosines", "pUC19": "CpG Only"}
DISPLAY_LABEL = {"P3": "P1", "P4": "P2"}


def fmt_count(x, _pos=None):
    for div, suf in ((1e9, "B"), (1e6, "M"), (1e3, "k")):
        if x >= div:
            return f"{x / div:g}{suf}"
    return f"{x:g}"


def plot_one(d, spike, output):
    label, fn = PERFORMANCE[spike]
    plates = sorted(d["plate"].unique())
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    rng = np.random.default_rng(0)

    ax.boxplot([d.loc[(d.plate == p) & (d.spike_type == spike), "cov"] for p in plates],
               positions=range(len(plates)), widths=0.5, showfliers=False,
               boxprops=dict(color="black"), whiskerprops=dict(color="black"),
               capprops=dict(color="black"), medianprops=dict(color="black", linewidth=1.5))

    sc = None
    for i, p in enumerate(plates):
        s = d[(d.plate == p) & (d.spike_type == spike)]
        sc = ax.scatter(i + rng.uniform(-0.15, 0.15, len(s)), s["cov"],
                        c=fn(s["beta"].to_numpy()), cmap=CMAP, vmin=0, vmax=1,
                        s=30, alpha=0.8, edgecolors="none")

    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_count))
    ax.set_xticks(range(len(plates)))
    ax.set_xticklabels([DISPLAY_LABEL.get(p, p) for p in plates])
    ax.set_xlim(-0.5, len(plates) - 0.5)
    ax.set_ylabel(f"Coverage ({CONTEXT_LABEL[spike]})")
    ax.set_title(f"{spike} Spike-In Coverage", loc="left", fontsize=14)
    for spine in ax.spines.values():
        spine.set_color("grey")
    ax.grid(axis="y", color="lightgrey", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.colorbar(sc, ax=ax, label=label, shrink=0.8)

    os.makedirs(os.path.dirname(output), exist_ok=True)
    save_figure(fig, output, dpi=150, bbox_inches="tight")
    print(f"Wrote {output}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--table", required=True)
    p.add_argument("--output-lambda", required=True)
    p.add_argument("--output-puc19", required=True)
    args = p.parse_args()
    d = pd.read_csv(args.table, sep="\t")
    plot_one(d, "Lambda", args.output_lambda)
    plot_one(d, "pUC19", args.output_puc19)


if __name__ == "__main__":
    main()
