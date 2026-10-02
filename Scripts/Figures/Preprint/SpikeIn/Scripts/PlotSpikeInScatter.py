#!/usr/bin/env python3
"""
One dot per well: Lambda conversion against pUC19 methylation, coloured by plate.

Both controls should be near 1 in a well that worked -- lambda fully converted and pUC19
recovered as methylated -- so the cloud piles into the top right corner and the plates sit
on top of each other. Dots are therefore drawn strongly transparent, so overlap builds up
rather than one plate hiding the other, and a well that failed either control falls out of
the corner along the axis it failed.

Plate colours are the ones the CovBoxplots and CellTypeCounts figures use, defined in
Shared/FirstPrepCoverage.R; they are repeated here because that file is R.
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
import pandas as pd

PLATE_COLORS = {"P3": "#8B7CC0", "P4": "#D9A441"}
DISPLAY_LABEL = {"P3": "P1", "P4": "P2"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--table", required=True)
    p.add_argument("--alpha", type=float, default=0.35)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = pd.read_csv(args.table, sep="\t")
    w = d.pivot_table(index=["plate", "well"], columns="spike_type", values="beta").reset_index()
    w = w.dropna(subset=["Lambda", "pUC19"])
    w["conversion"] = 1 - w["Lambda"]

    fig, ax = plt.subplots(figsize=(6, 6))
    for plate, col in PLATE_COLORS.items():
        s = w[w["plate"] == plate]
        ax.scatter(s["conversion"], s["pUC19"], s=42, color=col, alpha=args.alpha,
                   edgecolors="none", label=DISPLAY_LABEL.get(plate, plate))
        print(f"{plate}: {len(s)} wells, median conversion {s.conversion.median():.4f}, "
              f"median pUC19 methylation {s['pUC19'].median():.4f}", flush=True)

    ax.set_xlabel("Lambda conversion  (1 - beta, all cytosines)")
    ax.set_ylabel("pUC19 methylation  (beta, CpG only)")
    ax.set_title("Spike-in performance per well", loc="left", fontsize=13)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_aspect("equal")
    ax.grid(color="lightgrey", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("grey")
    # matplotlib centres a legend title over the handles by default, which leaves it
    # floating away from the entries; align it with them.
    leg = ax.legend(title="Plate", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    leg._legend_box.align = "left"

    save_figure(fig, args.output, dpi=150, bbox_inches="tight")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
