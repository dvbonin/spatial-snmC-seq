#!/usr/bin/env python3
"""
Metagene methylation profile per expression group, keratinocytes next to Merkel
carcinoma. One figure per region type (--mode genebody | promoter). --n-groups sets how
many equal-sized groups the genes are cut into (4 = quartiles, 5 = quintiles).

Genes are ranked by that cell type's own bulk RNA-seq level, so the groups are
within-panel and the two panels are NOT on a common expression scale. Worth remembering when reading Q1: the keratinocyte reference
has only 4 samples and 20.6% of its genes sit at exactly zero counts, which fills Q1
almost entirely with undetected genes; the Merkel reference has 8 samples and one such
gene, so its Q1 really is "lowest expressed".

Every gene counts once per bin: the CpGs a gene has in a bin are averaged into one value
for that gene, and the curve is the plain mean of those per-gene values across the
quartile. Nothing is weighted by coverage or by CpG count, so a long CpG-dense gene does
not outvote a sparse one -- at the cost of noise, since a gene whose only CpG in a bin
comes from a single cell contributes a full-weight 0 or 1.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import sys
sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import FONT_STACK, apply_style, save_figure
apply_style()
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd

SUBTYPE = {"HealthyKeratinocytes": "bK", "MerkelCarcinoma": "cMCC"}
DISPLAY = {"HealthyKeratinocytes": "Keratinocytes", "MerkelCarcinoma": "Merkel Cell Carcinoma"}
def labels_and_colors(n):
    lab = {q: f"Q{q}" for q in range(1, n + 1)}
    lab[1] += " (lowest expression)"
    lab[n] += " (highest expression)"
    return lab, {q: cm.viridis((q - 1) / (n - 1)) for q in range(1, n + 1)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kz-profile", required=True)
    p.add_argument("--mcc-profile", required=True)
    p.add_argument("--expression", required=True)
    p.add_argument("--mode", required=True, choices=["genebody", "promoter"])
    p.add_argument("--n-groups", type=int, default=4)
    p.add_argument("--flank", type=int, default=5000, help="promoter only, for the title")
    p.add_argument("--output", required=True)
    args = p.parse_args()
    n = args.n_groups
    LABELS, COLORS = labels_and_colors(n)
    groups = list(range(1, n + 1))

    expr = pd.read_csv(args.expression)
    cols = ["gene", "bin", "mean_beta"]      # n_cpg is written for reference, not used here
    profiles = {"HealthyKeratinocytes": pd.read_parquet(args.kz_profile, columns=cols),
                "MerkelCarcinoma": pd.read_parquet(args.mcc_profile, columns=cols)}

    binned_all = []
    for cell_type, df in profiles.items():
        e = expr[expr["subtype"] == SUBTYPE[cell_type]].set_index("gene")["mean_lognorm"]
        genes = [g for g in df["gene"].unique() if g in e.index]
        quartile = pd.qcut(e.loc[genes].rank(method="first"), n, labels=groups)

        df = df[df["gene"].isin(quartile.index)].copy()
        df["quartile"] = df["gene"].map(quartile)
        binned = df.groupby(["quartile", "bin"], observed=True)["mean_beta"].mean().reset_index()
        binned["cell_type"] = cell_type
        binned_all.append(binned)
        print(f"{cell_type}: {len(genes)} genes grouped "
              f"(of {df['gene'].nunique()} with any covered CpG in the region)", flush=True)

    binned_all = pd.concat(binned_all, ignore_index=True)
    binned_all.to_csv(args.output.replace(".png", "_binned.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), sharex=True, sharey=True)
    for ax, cell_type in zip(axes, ["HealthyKeratinocytes", "MerkelCarcinoma"]):
        for q in groups:
            s = binned_all[(binned_all["quartile"] == q) &
                           (binned_all["cell_type"] == cell_type)].sort_values("bin")
            ax.plot(s["bin"], s["mean_beta"], color=COLORS[q], linewidth=1.5, label=LABELS[q])
        if args.mode == "promoter":
            ax.axvline(0, color="grey", linestyle=":", linewidth=1)
        ax.set_title(DISPLAY[cell_type])
        ax.set_xlabel("Position Relative to TSS (bp)" if args.mode == "promoter"
                      else "Normalized Gene Body Position (0=TSS, 1=TES)")
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel("Mean CpG Methylation")
    axes[0].legend(frameon=False)
    grp = {4: "Quartile", 5: "Quintile"}.get(n, f"{n}-Group")
    fig.suptitle(f"Promoter Methylation by Expression {grp} (TSS +/-{args.flank}bp)"
                 if args.mode == "promoter" else f"Gene Body Methylation by Expression {grp}")
    fig.tight_layout()
    save_figure(fig, args.output, dpi=150)
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
