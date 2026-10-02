#!/usr/bin/env python3
"""
Plot read-count retention across pipeline steps.

Single stream: --summary only — one line per run coloured by quadrant.
Two streams:   --summary + --summary2 — solid vs dashed lines, labelled by
               --label1 / --label2 (e.g. Lambda vs pUC19).
"""
import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

METADATA      = {"run", "sample", "plate", "chunk"}
CHUNK_COLORS  = ["#e15759", "#4e79a7", "#59a14f", "#f28e2b", "#b07aa1", "#ff9da7",
                 "#9c755f", "#bab0ac"]


def _fmt_count(x, _):
    if x >= 1e9: return f"{x/1e9:.3g}B"
    if x >= 1e6: return f"{x/1e6:.3g}M"
    if x >= 1e3: return f"{x/1e3:.3g}k"
    return f"{x:.3g}"


def _add_lines(axes, df, col_names, step_labels, linestyle, label_prefix, seen,
               drop_first=False):
    min_frac = float("inf")
    max_frac = 0.0
    for i, (_, row) in enumerate(df.iterrows()):
        q        = str(row["chunk"])
        col      = CHUNK_COLORS[i % len(CHUNK_COLORS)]
        counts   = [row[c] for c in col_names]
        baseline = counts[0] if counts[0] > 0 else 1
        fracs    = [c / baseline for c in counts]
        if drop_first:
            counts      = counts[1:]
            fracs       = fracs[1:]
            xlabels     = step_labels[1:]
        else:
            xlabels     = step_labels
        if fracs:
            min_frac = min(min_frac, min(fracs))
            max_frac = max(max_frac, max(fracs))
        key      = f"{label_prefix} {q}" if label_prefix else q
        label    = key if key not in seen else None
        seen.add(key)
        kw = dict(color=col, alpha=0.75, linewidth=1.4,
                  marker="o", markersize=4, linestyle=linestyle)
        axes[0].plot(xlabels, counts, label=label, **kw)
        axes[1].plot(xlabels, fracs,  label=label, **kw)
    if min_frac == float("inf"):
        min_frac = 0.0
    return min_frac, max_frac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary",    required=True)
    ap.add_argument("--summary2",   default=None)
    ap.add_argument("--label1",     default="")
    ap.add_argument("--label2",     default="")
    ap.add_argument("--output",     required=True)
    ap.add_argument("--title",      default="")
    ap.add_argument("--drop-first", action="store_true",
                    help="Use first column as retention baseline but exclude from x-axis")
    args = ap.parse_args()

    df1   = pd.read_csv(args.summary, sep="\t")
    cols1 = [c for c in df1.columns if c not in METADATA]
    labs1 = [c.replace("_", " ").title() for c in cols1]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    seen               = set()
    min_frac, max_frac = _add_lines(axes, df1, cols1, labs1, "-", args.label1, seen,
                                    drop_first=args.drop_first)

    if args.summary2:
        df2   = pd.read_csv(args.summary2, sep="\t")
        cols2 = [c for c in df2.columns if c not in METADATA]
        labs2 = [c.replace("_", " ").title() for c in cols2]
        min_f2, max_f2 = _add_lines(axes, df2, cols2, labs2, "--", args.label2,
                                     seen, drop_first=args.drop_first)
        min_frac = min(min_frac, min_f2)
        max_frac = max(max_frac, max_f2)

    axes[0].yaxis.set_major_formatter(mticker.FuncFormatter(_fmt_count))
    axes[0].set_ylabel("Read pairs")
    axes[0].set_title("Absolute read counts")

    axes[1].yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{x*100:.3g}%"
    ))
    axes[1].set_ylim(max(0.0, min_frac * 0.8), max(max_frac, 0.01) * 1.2)
    axes[1].set_ylabel("Fraction of Leftover reads" if args.drop_first
                       else "Fraction of starting reads retained")
    axes[1].set_title("Read retention")

    for ax in axes:
        ax.set_xlabel(None)
        ax.tick_params(axis="x", rotation=30)
        ax.grid(axis="y", linestyle=":", alpha=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Chunk",
               loc="upper right", bbox_to_anchor=(0.99, 0.95))

    title = "Read count tracking"
    if args.title:
        title += f" — {args.title}"
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 0.88, 1])

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
