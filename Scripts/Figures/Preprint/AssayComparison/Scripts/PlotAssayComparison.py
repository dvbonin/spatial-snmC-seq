#!/usr/bin/env python3
"""
Our assay against published single-cell methylation methods, on matched conventions.

  AssayComparison_CpGsPerCell       covered CpGs per cell (Extended Data Fig. 1k)
  AssayComparison_CpGsPerReadPair   covered CpGs per raw read pair (Extended Data Fig. 1l)
  AssayComparison_ReadsVsCpGs       raw read pairs against covered CpGs per cell, log-log,
                                    with a linear fit per method (Fig. 1e)
  AssayComparison_summary.tsv       medians, cell numbers and fitted slopes behind the text

CpGs are strand-merged dinucleotides throughout (CountPublishedCells.py for the published
methods; mergecg dyads in our pipeline). Reads are raw read PAIRS throughout. Our
well_summary tables carry RawReads = 2 x read pairs (MergeRawBarcodeCounts.py), so they
are halved here -- comparing them unhalved against SRA spots counts every one of our
pairs twice and halves our apparent efficiency. Summaries carrying RawReadPairs are read
directly.

Our assay enters as QC-passed single-cell wells only (Cells == 1 & PassesQC).

sci-MET v3 has no per-cell read counts (SciMETReadBound.py), only an upper bound on reads
per cell. It is drawn as a bound, not as data: in the read-pair panel as a lower bound on
CpGs per pair (marker with an upward arrow), in the scatter as a marker at the upper bound
on reads with the interquartile range of its CpGs (the legend says "read bound"). Its CpG
distribution is complete (44,840 cells); 350 random cells are drawn as points for
legibility, the box uses all of them.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import apply_style, save_figure
apply_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

OURS = "spatial snmC-seq"
# label, year, colour; ordered by median CpGs per cell
METHODS = [
    ("scBS-seq",   2014, "#eb6834"),
    ("snmC-seq2",  2018, "#1baf7a"),
    (OURS,         2026, "#2a78d6"),
    ("snmC-seq3",  2023, "#4a3aa7"),
    ("sci-MET v3", 2025, "#6b6a66"),
]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e6e2"
MAX_POINTS = 350
FS = 7


def load_ours(specs):
    frames = []
    for spec in specs:
        plate, path = spec.split("=", 1)
        d = pd.read_csv(path)
        d["Plate"] = plate
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    passed = d["PassesQC"].astype(str).str.upper() == "TRUE"
    d = d[(d["Cells"] == 1) & passed].copy()
    if "RawReadPairs" in d:
        d["raw_read_pairs"] = d["RawReadPairs"]
    else:  # summaries written before RawReadPairs was added carry RawReads only
        if (d["RawReads"] % 2).any():
            raise SystemExit("RawReads is expected to be 2 x read pairs (even); found odd values")
        d["raw_read_pairs"] = d["RawReads"] // 2
    return pd.DataFrame({"cpgs": d["CpGs"].to_numpy(), "pairs": d["raw_read_pairs"].to_numpy()})


def load_published(path):
    d = pd.read_csv(path, sep="\t")
    return pd.DataFrame({"cpgs": d["covered_dinucleotides"].to_numpy(),
                         "pairs": d["raw_read_pairs"].to_numpy()})


def style(ax, grid_axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=grid_axis, color=GRID, lw=0.5, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.6)
    ax.tick_params(length=2.5, width=0.6, colors=INK2, labelsize=FS)


def box(ax, i, vals, color, rng):
    bp = ax.boxplot([vals], positions=[i], widths=0.62, showfliers=False, patch_artist=True,
                    medianprops=dict(color=INK, lw=1.2),
                    whiskerprops=dict(color=INK2, lw=0.6), capprops=dict(color=INK2, lw=0.6),
                    boxprops=dict(lw=0.6))
    bp["boxes"][0].set_facecolor(color)
    bp["boxes"][0].set_alpha(0.28)
    bp["boxes"][0].set_edgecolor(color)
    v = vals if len(vals) <= MAX_POINTS else rng.choice(vals, MAX_POINTS, replace=False)
    ax.scatter(rng.normal(i, 0.055, len(v)), v, s=2.5, color=color, alpha=0.45,
               edgecolor="none", zorder=3)
    return max(bp["caps"][1].get_ydata()[0], np.max(v))


def label_above(ax, i, top, text, **kw):
    """Label centred above a box, clear of its whisker and points."""
    style = dict(fontsize=FS - 1, fontweight="bold", color=INK)
    style.update(kw)
    ax.annotate(text, xy=(i, top), xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", **style)


def xticks(ax, ns=None):
    """Method and year per tick; with ns, the cell number on a third line."""
    labels = [f"{m.replace(' ', chr(10), 1) if m == OURS else m}\n({y})" for m, y, _ in METHODS]
    if ns is not None:
        labels = [f"{l}\nn={n:,}" for l, n in zip(labels, ns)]
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(labels, fontsize=FS - 1)
    ax.set_xlim(-0.6, len(METHODS) - 0.4)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--well-summary", nargs="+", required=True,
                   help="PLATE=path to FirstPrep_<plate>_well_summary.csv")
    p.add_argument("--scbs", required=True)
    p.add_argument("--snmc2", required=True)
    p.add_argument("--snmc3", required=True)
    p.add_argument("--scimet-cells", required=True)
    p.add_argument("--scimet-bound", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    data = {
        "scBS-seq": load_published(args.scbs),
        "snmC-seq2": load_published(args.snmc2),
        OURS: load_ours(args.well_summary),
        "snmC-seq3": load_published(args.snmc3),
    }
    scimet = pd.read_csv(args.scimet_cells, sep="\t")["cg_cov"].to_numpy()
    bound = pd.read_csv(args.scimet_bound, sep="\t").iloc[0]
    sm_pairs_hi = float(bound["read_pairs_per_cell_upper"])
    sm_cpr_lo = float(np.median(scimet)) / sm_pairs_hi
    rng = np.random.default_rng(0)
    os.makedirs(args.out_dir, exist_ok=True)
    summary = []

    # ---- covered CpGs per cell ----
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    ns = []
    for i, (m, _, c) in enumerate(METHODS):
        vals = scimet if m == "sci-MET v3" else data[m]["cpgs"].to_numpy()
        top = box(ax, i, vals, c, rng)
        label_above(ax, i, top, f"{np.median(vals):,.0f}")
        ns.append(len(vals))
    ax.set_yscale("log")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.8)            # headroom for the labels
    xticks(ax, ns)
    ax.set_ylabel("CpGs per cell (merged dinucleotides)", fontsize=FS)
    ax.set_title("(1) CpGs per cell", fontsize=FS + 1)
    style(ax)
    fig.tight_layout()
    save_figure(fig, os.path.join(args.out_dir, "AssayComparison_CpGsPerCell.png"), dpi=300)
    plt.close(fig)

    # ---- covered CpGs per raw read pair ----
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    for i, (m, _, c) in enumerate(METHODS):
        if m == "sci-MET v3":
            ax.scatter([i], [sm_cpr_lo], s=22, facecolors="none", edgecolors=c,
                       linewidths=1.0, zorder=5)
            ax.annotate("", xy=(i, sm_cpr_lo + 0.22), xytext=(i, sm_cpr_lo + 0.03),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9,
                                        mutation_scale=6))
            label_above(ax, i, sm_cpr_lo + 0.22, f"≥{sm_cpr_lo:.2f}\n(bound)",
                        fontweight="normal", fontstyle="italic", color=INK2)
            continue
        cpr = (data[m]["cpgs"] / data[m]["pairs"]).to_numpy()
        top = box(ax, i, cpr, c, rng)
        label_above(ax, i, top, f"{np.median(cpr):.2f}")
    xticks(ax)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
    ax.set_ylabel("CpGs per raw read pair", fontsize=FS)
    ax.set_title("(2) CpGs per raw read pair", fontsize=FS + 1)
    style(ax)
    fig.tight_layout()
    save_figure(fig, os.path.join(args.out_dir, "AssayComparison_CpGsPerReadPair.png"), dpi=300)
    plt.close(fig)

    # ---- raw read pairs vs covered CpGs ----
    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    handles = []
    for m, _, c in METHODS:
        if m == "sci-MET v3":
            lo, med, hi = np.percentile(scimet, [25, 50, 75])
            ax.plot([sm_pairs_hi, sm_pairs_hi], [lo, hi], color=c, lw=0.9, zorder=4)
            ax.scatter([sm_pairs_hi], [med], s=22, facecolors="none", edgecolors=c,
                       linewidths=1.0, zorder=5)
            handles.append(Line2D([0], [0], marker="o", ls="", mfc="none", mec=c, ms=4,
                                  label=f"{m} (read bound)"))
            summary.append(dict(method=m, n=len(scimet), median_cpgs=med,
                                median_read_pairs=f"<={sm_pairs_hi:.0f}",
                                median_cpg_per_pair=f">={sm_cpr_lo:.3f}", loglog_slope=""))
            continue
        d = data[m]
        d = d[(d["cpgs"] > 0) & (d["pairs"] > 0)]
        x, y = d["pairs"].to_numpy(), d["cpgs"].to_numpy()
        b, a = np.polyfit(np.log10(x), np.log10(y), 1)
        show = rng.choice(len(x), MAX_POINTS, replace=False) if len(x) > MAX_POINTS else slice(None)
        ax.scatter(x[show], y[show], s=3, color=c, alpha=0.45, edgecolor="none", zorder=3)
        gx = np.linspace(np.log10(x).min(), np.log10(x).max(), 50)
        ax.plot(10 ** gx, 10 ** (a + b * gx), color=c, lw=1.2, zorder=6, solid_capstyle="round")
        handles.append(Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=4, label=m))
        summary.append(dict(method=m, n=len(x), median_cpgs=np.median(y),
                            median_read_pairs=np.median(x),
                            median_cpg_per_pair=round(float(np.median(y / x)), 3),
                            loglog_slope=round(float(b), 3)))
    handles.append(Line2D([0], [0], color=INK2, lw=1.2, label="Linear fit (log-log)"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(left=sm_pairs_hi / 1.6)                # keep the bound marker off the axis
    ax.set_xlabel("Raw read pairs per cell", fontsize=FS)
    ax.set_ylabel("CpGs per cell (merged dinucleotides)", fontsize=FS)
    ax.set_title("(3) Raw read pairs vs. CpGs per cell", fontsize=FS + 1)
    ax.legend(handles=handles, fontsize=FS - 1.5, frameon=False, loc="lower right",
              handletextpad=0.3, borderaxespad=0.2)
    style(ax, grid_axis="both")
    fig.tight_layout()
    save_figure(fig, os.path.join(args.out_dir, "AssayComparison_ReadsVsCpGs.png"), dpi=300)
    plt.close(fig)

    s = pd.DataFrame(summary)
    s.to_csv(os.path.join(args.out_dir, "AssayComparison_summary.tsv"), sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
