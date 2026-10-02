#!/usr/bin/env python3
"""
Per-well spike-in heatmap for a single (sample, plate).

--metric methylation : coverage-weighted mean BT (or 1-BT if --invert)
--metric coverage    : total sum of CV across selected sites
--cx CG|all          : restrict to CpG sites or use all cytosines
--invert             : plot 1 - value (for lambda conversion rate)
"""
import argparse
import glob
import gzip
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

ROWS = list("ABCDEFGHIJKLMNOP")
COLS = list(range(1, 25))


def _well_to_rc(well):
    r = ord(well[0].upper()) - ord("A")
    c = int(well[1:]) - 1
    return r, c


def _well_value(vcf_path, metric, cx, invert, min_cov=0):
    total_cov, weighted = 0, 0.0
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue
            if cx == "CG":
                info = {}
                for item in fields[7].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v
                if info.get("CX") != "CG":
                    continue
            fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
            try:
                cv = int(fmt["CV"])
                bt = float(fmt["BT"])
                if cv > 0:
                    total_cov += cv
                    weighted  += bt * cv
            except (KeyError, ValueError):
                continue
    if total_cov == 0:
        return None
    if min_cov > 0 and total_cov < min_cov:
        return None
    if metric == "coverage":
        return float(total_cov)
    val = weighted / total_cov
    return (1.0 - val) if invert else val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf-dir",  required=True)
    ap.add_argument("--output",   required=True)
    ap.add_argument("--title",    default="")
    ap.add_argument("--metric",   choices=["methylation", "coverage"], default="methylation")
    ap.add_argument("--cx",       choices=["CG", "all"], default="CG")
    ap.add_argument("--invert",   action="store_true",
                    help="Plot 1 - value (e.g. lambda conversion rate)")
    ap.add_argument("--min-cov",  type=int, default=0,
                    help="Exclude wells with total cytosine coverage below this threshold")
    args = ap.parse_args()

    matrix = np.full((16, 24), np.nan)
    for path in glob.glob(os.path.join(args.vcf_dir, "*_VCF.gz")):
        well = os.path.basename(path).replace("_VCF.gz", "")
        if len(well) < 2:
            continue
        try:
            r, c = _well_to_rc(well)
        except (ValueError, IndexError):
            continue
        if not (0 <= r < 16 and 0 <= c < 24):
            continue
        val = _well_value(path, args.metric, args.cx, args.invert, args.min_cov)
        if val is not None:
            matrix[r, c] = val

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    if args.metric == "coverage":
        vals = matrix[~np.isnan(matrix)]
        vmax = float(np.percentile(vals, 95)) if len(vals) > 0 else 1.0
        im = ax.imshow(matrix, aspect="auto", cmap="plasma",
                       vmin=0, vmax=vmax, interpolation="none")
        cbar_label = "Total cytosine coverage"
    else:
        im = ax.imshow(matrix, aspect="auto", cmap="viridis",
                       vmin=0, vmax=1, interpolation="none")
        if args.invert:
            cbar_label = ("CpG conversion rate" if args.cx == "CG"
                          else "All-context conversion rate")
        else:
            cbar_label = ("CpG methylation rate" if args.cx == "CG"
                          else "All-context methylation rate")

    ax.set_xticks(range(24))
    ax.set_xticklabels([str(c) for c in COLS], fontsize=7)
    ax.set_yticks(range(16))
    ax.set_yticklabels(ROWS, fontsize=7)
    ax.xaxis.tick_top()
    cbar = plt.colorbar(im, ax=ax, label=cbar_label)
    if args.metric == "coverage":
        cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: (f"{x/1e9:.3g}B" if x >= 1e9 else
                          f"{x/1e6:.3g}M" if x >= 1e6 else
                          f"{x/1e3:.3g}k" if x >= 1e3 else
                          f"{x:.3g}")
        ))
    title = args.title or cbar_label
    ax.set_title(f"{title} (min_cov={args.min_cov})")
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    plt.close()
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
