#!/usr/bin/env python3
"""
Read-yield funnel per plate, from Results/ReadCounts/{sample}/{plate}/Main_read_counts.csv
(one row per sequencing chunk, written by CollectReadCounts.py), summed per plate.

Two variants, selected with --stages:

  main  Sequenced / Aligned / Deduplicated, drawn as overlaid bars. The headline
        figure. Note that the Aligned -> Deduplicated step folds in two smaller
        losses besides duplicates -- reads whose barcode did not match the
        whitelist, and reads whose mate did not survive the MAPQ filter -- each
        around 1% or less.

  all   Every stage the pipeline records, drawn as grouped bars with each bar's
        share of the sequenced total underneath. Overlaid bars cannot work here:
        consecutive stages differ by well under 1% of the axis, so their bands
        would be invisible. This variant is what shows where reads are actually
        lost -- trimming in particular can remove far more than alignment does.
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import sys
sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import FONT_STACK, apply_style, save_figure
apply_style()
matplotlib.rcParams["font.family"] = FONT_STACK[0]
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

STAGE_LABELS = {
    "raw":     "Sequenced",
    "trimmed": "Trimmed",
    "aligned": "Aligned",
    "tagged":  "Barcode-assigned",
    "paired":  "Mate-pair complete",
    "deduped": "Deduplicated (Unique)",
}
STAGE_SETS = {
    "main": ["raw", "aligned", "deduped"],
    "all":  ["raw", "trimmed", "aligned", "tagged", "paired", "deduped"],
}
MAIN_COLORS = ["#BFD7EA", "#5B9BD5", "#1F4E79"]
ALL_COLORS = ["#DCE9F2", "#BFD7EA", "#8CB8DC", "#5B9BD5", "#31719E", "#1F4E79"]

# Display-only rename on the x axis: plate 3 -> "P1", plate 4 -> "P2".
DISPLAY_LABEL = {"P3": "P1", "P4": "P2"}


def plate_totals(csv_path, stages):
    """Sum each stage over the plate's sequencing chunks. The file is tab-separated
    despite the .csv name."""
    totals = {s: 0 for s in stages}
    with open(csv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for s in stages:
                totals[s] += int(row[s])
    return totals


def plot_overlaid(ax, plates, stages, colors):
    """Concentric bars, tallest stage behind. Values are labelled above each
    band's top, except the final stage which is labelled inside it in white."""
    spacing, width = 0.45, 0.3
    x = [i * spacing for i in range(len(plates))]
    max_raw = max(t[stages[0]] for _, t in plates)

    for pos, (label, totals) in zip(x, plates):
        for j, (stage, color) in enumerate(zip(stages, colors)):
            ax.bar(pos, totals[stage], width=width, color=color, zorder=j + 1,
                   label=STAGE_LABELS[stage] if pos == x[0] else None)

        first, mid, last = (totals[s] for s in stages)
        ax.text(pos, first + max_raw * 0.015, f"{first / 1e6:.0f}M",
                ha="center", va="bottom", fontsize=8, color="#333333")
        # If the middle and final bands are too close for two separate labels,
        # combine them onto one line above the bar.
        if (mid - last) < max_raw * 0.05:
            ax.text(pos, mid + max_raw * 0.015, f"{mid/1e6:.0f}M / {last/1e6:.0f}M",
                    ha="center", va="bottom", fontsize=8, color="#333333")
        else:
            ax.text(pos, mid + max_raw * 0.015, f"{mid / 1e6:.0f}M",
                    ha="center", va="bottom", fontsize=8, color="#333333")
            ax.text(pos, last - max_raw * 0.015, f"{last / 1e6:.0f}M",
                    ha="center", va="top", fontsize=8, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY_LABEL.get(l, l) for l, _ in plates])
    # Half a bar-width of margin either side: with only a couple of plates a
    # wider axis leaves the bars floating in white space.
    ax.set_xlim(x[0] - width * 0.75, x[-1] + width * 0.75)
    return 3


def plot_grouped(ax, plates, stages, colors):
    """One bar per stage per plate, with each stage's share of the sequenced
    total written underneath its value."""
    n = len(stages)
    bar_w, group_gap = 0.8 / n, 1.0
    centres = [i * group_gap for i in range(len(plates))]
    ymax = max(t[stages[0]] for _, t in plates)

    for gi, (label, totals) in enumerate(zip([p[0] for p in plates], [p[1] for p in plates])):
        base = centres[gi] - 0.4 + bar_w / 2
        for j, (stage, color) in enumerate(zip(stages, colors)):
            pos = base + j * bar_w
            ax.bar(pos, totals[stage], width=bar_w * 0.9, color=color, zorder=2,
                   label=STAGE_LABELS[stage] if gi == 0 else None)
            pct = 100 * totals[stage] / totals[stages[0]]
            ax.text(pos, totals[stage] + ymax * 0.012,
                    f"{totals[stage]/1e6:.0f}M\n{pct:.0f}%",
                    ha="center", va="bottom", fontsize=7, color="#333333",
                    linespacing=1.15)

    ax.set_xticks(centres)
    ax.set_xticklabels([DISPLAY_LABEL.get(l, l) for l, _ in plates])
    ax.set_xlim(-0.5, centres[-1] + 0.5)
    ax.set_ylim(0, ymax * 1.22)
    return 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read-counts", nargs="+", required=True,
                    help="plate_label=path/to/Main_read_counts.csv, one per plate")
    ap.add_argument("--stages", choices=list(STAGE_SETS), default="main")
    ap.add_argument("--title", default="Read Yield")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    stages = STAGE_SETS[args.stages]
    plates = []
    for spec in args.read_counts:
        label, path = spec.split("=", 1)
        plates.append((label, plate_totals(path, stages)))

    # Width follows the plate count so the bars fill the panel; the constant
    # covers the y axis and the legend's widest entry.
    if args.stages == "all":
        fig, ax = plt.subplots(figsize=(3.6 + 1.7 * len(plates), 6.2))
        ncol = plot_grouped(ax, plates, stages, ALL_COLORS)
    else:
        fig, ax = plt.subplots(figsize=(3.6 + 1.1 * len(plates), 6))
        ncol = plot_overlaid(ax, plates, stages, MAIN_COLORS)

    ax.set_ylabel("Read pairs")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e6:.0f}M"))
    ax.set_title(args.title)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=ncol)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    save_figure(fig, args.output, dpi=150)
    for label, totals in plates:
        drops = "  ".join(f"{s}={totals[s]/1e6:.1f}M" for s in stages)
        print(f"  {label}: {drops}")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
