#!/usr/bin/env python3
import argparse
import os

import matplotlib
matplotlib.rcParams["font.family"] = "Nimbus Sans"
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Beeswarm-style version of the per-well "mean fraction empty in own cluster"
# distribution from SweepLeidenClusterUMAP.R -- ggplot2/ggbeeswarm isn't
# installable in the shared R conda env (its Makeconf can't compile C code),
# so this one plot is done in Python instead, reading the cached CSV.
#
# The data is heavily tied (only ~25 distinct values across ~768 wells, one
# value shared by nearly a quarter of them), so a circle-packing beeswarm
# either overlaps badly or needs huge marker-to-width ratios. Instead, each
# well is drawn as a thin vertical tick, and wells sharing the exact same
# value are spread deterministically (not randomly) across a row -- ticks
# need far less horizontal room per point than circles, so far more of them
# fit without overlapping, and row width still reflects how many wells share
# that value.

parser = argparse.ArgumentParser()
parser.add_argument("--well-empty-frac-csv", required=True)
parser.add_argument("--n-sweeps", type=int, required=True, help="len(n_neighbors) x len(resolution) from the sweep")
parser.add_argument("--output", required=True)
args = parser.parse_args()

well_rank = pd.read_csv(args.well_empty_frac_csv)

cell_count_levels = sorted(set([0, 1, 2, 5, 10]) | set(well_rank["CellCount"].unique()))
non_zero_levels = [c for c in cell_count_levels if c != 0]
cell_count_colors = {0: "#E8998D"}
ramp = matplotlib.colors.LinearSegmentedColormap.from_list("ramp", ["#89C2D9", "#2F6B5E"])
for i, c in enumerate(non_zero_levels):
    cell_count_colors[c] = matplotlib.colors.to_hex(ramp(i / max(len(non_zero_levels) - 1, 1)))

max_frac = well_rank["mean_empty_frac"].max()
y_max = np.ceil(max_frac / 0.05) * 0.05
y_breaks = np.arange(0, y_max + 1e-9, 0.05)

# ---- Deterministic per-row spread: same step size per tick regardless of ----
# ---- row size, so a row's width still reflects how many wells tie there ----
max_group_n = well_rank.groupby("mean_empty_frac").size().max()
step = 1.6 / max_group_n

well_rank = well_rank.sort_values(["mean_empty_frac", "CellCount"]).reset_index(drop=True)
xs = np.empty(len(well_rank))
for _, idx in well_rank.groupby("mean_empty_frac").indices.items():
    n = len(idx)
    xs[idx] = (np.arange(n) - (n - 1) / 2) * step

colors = well_rank["CellCount"].map(cell_count_colors).to_numpy()

fig, ax = plt.subplots(figsize=(7, 9), dpi=150)
ax.scatter(xs, well_rank["mean_empty_frac"], marker="|", s=45, linewidths=0.8, c=colors, alpha=0.85)
half_width = (max_group_n - 1) / 2 * step
ax.set_xlim(-half_width - 0.1, half_width + 0.1)
ax.set_xticks([])
ax.set_ylim(-0.01, y_max)
ax.set_yticks(y_breaks)
ax.set_xlabel(None)
ax.set_ylabel(f"Mean Fraction Empty in Own Cluster ({args.n_sweeps} Sweeps)")
ax.set_title(f"Own-Cluster Empty Fraction per Well (n={len(well_rank)})")
ax.spines[["top", "right"]].set_visible(False)

legend_handles = [
    mlines.Line2D([], [], color=cell_count_colors[c], marker="|", linestyle="none", markersize=10, markeredgewidth=1.5)
    for c in reversed(cell_count_levels)
]
ax.legend(legend_handles, [str(c) for c in reversed(cell_count_levels)],
          title="Cell Count", bbox_to_anchor=(1.02, 1), loc="upper left")

fig.tight_layout()
fig.savefig(args.output, dpi=150)
# SVG alongside the PNG, off unless SNMC_SVG is set in the environment. Vector output
# is wanted only occasionally, and for dense scatters it can run to hundreds of MB.
if os.environ.get("SNMC_SVG"):
    fig.savefig( args.output.replace(".png", ".svg"), dpi=150)
print(f"Wrote: {args.output}")
