#!/usr/bin/env python3
"""
Every well as one point: CpG coverage against its mean fraction-empty-in-own-cluster
across the Leiden sweep, coloured by cell count.

The two QC signals plotted against each other. Coverage is the direct, per-well
measurement; the empty fraction is the clustering-derived one that the pass/fail
call actually uses. Wells that disagree between them -- high coverage but
clustering with empties, or vice versa -- are the interesting ones, and they are
only visible in this view.

Python rather than R because the shared R conda env cannot compile source
packages (see PlotEmptyFractionBeeswarm.py for the same reason).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import sys
sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import FONT_STACK, apply_style, save_figure
apply_style()
matplotlib.rcParams["font.family"] = FONT_STACK[0]
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--well-empty-frac-csv", required=True,
                    help="FirstPrep_Merged_well_empty_frac.csv from SweepLeidenClusterUMAP.R")
parser.add_argument("--well-stats-p3", required=True)
parser.add_argument("--well-stats-p4", required=True)
parser.add_argument("--threshold", type=float, default=None,
                    help="Draw the QC empty-fraction cutoff as a horizontal line")
parser.add_argument("--n-sweeps", type=int, required=True,
                    help="len(n_neighbors) x len(resolution) from the Leiden sweep")
parser.add_argument("--output", required=True)
args = parser.parse_args()

ef = pd.read_csv(args.well_empty_frac_csv)
ws = pd.concat([pd.read_csv(args.well_stats_p3).assign(Plate="P3"),
                pd.read_csv(args.well_stats_p4).assign(Plate="P4")])
d = ef.merge(ws[["Plate", "WellPosition", "CpGsCovered"]],
             left_on=["Plate", "well"], right_on=["Plate", "WellPosition"])

cell_count_levels = sorted(set([0, 1, 2, 5, 10]) | set(d["CellCount"].unique()))
non_zero_levels = [c for c in cell_count_levels if c != 0]
ramp = matplotlib.colors.LinearSegmentedColormap.from_list("ramp", ["#89C2D9", "#2F6B5E"])
colors = {0: "#E8998D"}
for i, c in enumerate(non_zero_levels):
    colors[c] = matplotlib.colors.to_hex(ramp(i / max(len(non_zero_levels) - 1, 1)))

fig, ax = plt.subplots(figsize=(9, 7), dpi=150)
# Plot the deepest wells first so the sparse, low-coverage ones stay visible on top.
for c in sorted(cell_count_levels, reverse=True):
    sub = d[d["CellCount"] == c]
    if not len(sub):
        continue
    ax.scatter(sub["CpGsCovered"], sub["mean_empty_frac"], s=26, alpha=0.75,
               c=colors[c], edgecolors="none", label=str(c))

if args.threshold is not None:
    ax.axhline(args.threshold, linestyle="--", linewidth=1, color="black")
    # Sits just above the line so the dashes do not strike through the text, and
    # right-aligned at the axes edge so the label ends exactly where the line
    # does. get_yaxis_transform() blends axes-fraction x with data-coordinate y,
    # so this is unaffected by the later switch to a log x scale.
    ax.text(1.0, args.threshold + 0.008, f"QC cutoff {args.threshold:g}",
            transform=ax.get_yaxis_transform(),
            va="bottom", ha="right", fontsize=9)

ax.set_xscale("log")
ax.set_xlabel("CpGs covered (log scale)")
ax.set_ylabel(f"Mean Fraction Empty in Own Cluster ({args.n_sweeps} Sweeps)")
ax.set_title("CpG coverage vs. own-cluster empty fraction, per well "
             f"(n={len(d)})")
ax.spines[["top", "right"]].set_visible(False)

handles = [mlines.Line2D([], [], color=colors[c], marker="o", linestyle="none", markersize=7)
           for c in reversed(cell_count_levels)]
ax.legend(handles, [str(c) for c in reversed(cell_count_levels)],
          title="Cell Count", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

fig.tight_layout()
save_figure(fig, args.output, dpi=150)
print(f"Wrote: {args.output}")
