#!/usr/bin/env python3
"""
Plate heatmaps of per-well read counts, one PNG per read-class stream
(Main / Alt / Lambda / pUC19), from the per-type `{type}_read_counts.csv` files
written by MergeWellChunks.py.

Each cell shows the read count and, in brackets, that count as a percentage of
the well's total across all streams -- so the streams sum to 100% per well by
construction. All streams must therefore be passed together even though each is
plotted separately. Reads that entered the split but were written to no well
(no XP tag, dropped singletons, tag/contig-inconsistent pairs) carry no well
assignment and cannot appear here; they are reported per chunk in the
GroupReadsToCells.py summary instead.

Wells absent from a stream's CSV genuinely had zero reads of that class (a
stream only lists wells it observed), so the grid is left-joined and 0-filled
rather than shown as missing.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Nimbus Sans"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROWS = [chr(c) for c in range(ord("A"), ord("P") + 1)]   # A-P
COLS = list(range(1, 25))                                 # 1-24
WELLS = [f"{r}{c:02d}" for r in ROWS for c in COLS]

# Tops out at a medium green rather than a dark one: every cell label is drawn in
# black, so the darkest end of the scale still has to keep black text readable.
GREENS = matplotlib.colors.LinearSegmentedColormap.from_list(
    "well_greens", ["#ffffff", "#e4f3e7", "#b7e0c2", "#84cb9c", "#57b77e"])


def compact(n):
    if n < 1000:
        return str(int(n))
    if n < 1_000_000:
        return f"{n / 1e3:.0f}k" if n >= 10_000 else f"{n / 1e3:.1f}k"
    return f"{n / 1e6:.1f}M"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", nargs="+", required=True,
                    help="One or more Name=path pairs, e.g. Main=/path/Main_read_counts.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True,
                    help="Output basename prefix, e.g. FirstPrep_P3")
    ap.add_argument("--suffix", default="well_counts",
                    help="Output basename suffix [default: %(default)s]")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    streams = []
    grid = pd.DataFrame(index=WELLS)
    for spec in args.counts:
        name, path = spec.split("=", 1)
        df = pd.read_csv(path, sep=";")
        counts = df.set_index("WellPosition")["ReadCounts"]
        grid[name] = counts.reindex(WELLS).fillna(0).astype("int64")
        streams.append(name)

    total = grid[streams].sum(axis=1)
    os.makedirs(args.outdir, exist_ok=True)

    for name in streams:
        counts = grid[name]
        mat = counts.to_numpy().reshape(len(ROWS), len(COLS)).astype(float)
        pct = np.where(total.to_numpy() > 0, 100 * counts.to_numpy() / total.to_numpy(), 0.0)
        pct = pct.reshape(len(ROWS), len(COLS))

        fig, ax = plt.subplots(figsize=(28, 19), dpi=150)

        # log colour scale: Main is ~10^5 reads/well while Alt is ~10^1, so a
        # linear scale would render every stream but Main as a flat block.
        # White -> medium green, deliberately stopping short of a dark green so
        # black cell labels stay legible across the whole range (the previous
        # dark-ended map forced a white/black text switch that was unreadable
        # in the mid tones).
        im = ax.imshow(np.log10(mat + 1), cmap=GREENS, aspect="auto")
        cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
        cb.set_label(r"$\log_{10}$(reads + 1)", fontsize=12)
        cb.ax.tick_params(labelsize=10)

        for i in range(len(ROWS)):
            for j in range(len(COLS)):
                if mat[i, j] == 0:
                    ax.text(j, i, "0", ha="center", va="center", fontsize=17, color="black")
                    continue
                ax.text(j, i - 0.17, compact(mat[i, j]), ha="center", va="center",
                        fontsize=17, color="black")
                ax.text(j, i + 0.23, f"({pct[i, j]:.1f}%)", ha="center", va="center",
                        fontsize=14, color="black")

        ax.set_xticks(range(len(COLS)), [str(c) for c in COLS], fontsize=12)
        ax.set_yticks(range(len(ROWS)), ROWS, fontsize=12)
        ax.set_xticks(np.arange(-0.5, len(COLS), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(ROWS), 1), minor=True)
        ax.grid(which="minor", color="#b0b0b0", linewidth=0.6)
        ax.tick_params(which="minor", length=0)

        subtitle = (f"{int(counts.sum()):,} reads "
                    f"({100 * counts.sum() / max(total.sum(), 1):.2f}% of plate)")
        ax.set_title(f"{name} per-well read counts"
                     f"{'  --  ' + args.title if args.title else ''}\n{subtitle}"
                     "\n(bracketed value = share of that well's assigned reads)",
                     fontsize=15)

        out = os.path.join(args.outdir, f"{args.prefix}_{name}_{args.suffix}.png")
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"Wrote {out} ({int(counts.sum()):,} reads)")


if __name__ == "__main__":
    main()
