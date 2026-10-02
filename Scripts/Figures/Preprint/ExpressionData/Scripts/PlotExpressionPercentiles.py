#!/usr/bin/env python3
"""
The expression distribution the metagene groups are cut from: genes ranked by their
bulk RNA-seq level, binned into percentiles of 1% of genes each, mean level per
percentile, with the group boundaries marked (--n-groups, 4 = quartiles).

Genes come from whichever profile files are passed in, so the marked boundaries are the
ones actually used by that metagene figure -- pass the promoter profiles to describe the
promoter figure, the gene-body ones for that panel (they differ by a few hundred genes). The two panels have independent y axes:
the two RNA-seq series are separate experiments with different depths, and their levels
are only ever used to rank genes within a cell type, never compared across panels.
"""
import os
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.environ["SNMC_ROOT"], "Scripts/Figures/Preprint/Shared"))
from figure_style import apply_style, save_figure
apply_style()

SUBTYPE = {"HealthyKeratinocytes": "bK", "MerkelCarcinoma": "cMCC"}
DISPLAY = {"HealthyKeratinocytes": "Keratinocytes", "MerkelCarcinoma": "Merkel Cell Carcinoma"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kz-profile", required=True)
    p.add_argument("--mcc-profile", required=True)
    p.add_argument("--expression", required=True)
    p.add_argument("--n-groups", type=int, default=4)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    expr = pd.read_csv(args.expression)
    profiles = {"HealthyKeratinocytes": pd.read_parquet(args.kz_profile, columns=["gene"]),
                "MerkelCarcinoma": pd.read_parquet(args.mcc_profile, columns=["gene"])}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    rows = []
    for ax, cell_type in zip(axes, ["HealthyKeratinocytes", "MerkelCarcinoma"]):
        e = expr[expr["subtype"] == SUBTYPE[cell_type]].set_index("gene")["mean_lognorm"]
        e = e.loc[e.index.intersection(profiles[cell_type]["gene"].unique())].sort_values()
        pct = np.ceil(np.arange(1, len(e) + 1) / len(e) * 100).astype(int)
        m = pd.DataFrame({"pct": pct, "v": e.to_numpy()}).groupby("pct")["v"].mean()

        ax.plot(m.index, m.to_numpy(), color="#1B6CA8", linewidth=1.6)
        n = args.n_groups
        edges = [100 * i / n for i in range(1, n)]
        for q in edges:
            ax.axvline(q, color="grey", linestyle="--", linewidth=0.8)
        ymax = m.max()
        for i in range(n):
            ax.text(100 * (i + 0.5) / n, ymax * 0.96, f"Q{i + 1}",
                    ha="center", va="top", color="grey", fontsize=11)
        ax.set_title(f"{DISPLAY[cell_type]}  ({len(e)} genes)")
        ax.set_xlabel("Percentile of expression-ranked genes")
        ax.set_xlim(0, 100)
        ax.spines[["top", "right"]].set_visible(False)

        cuts = np.percentile(e.to_numpy(), edges)
        print(f"{cell_type}: {len(e)} genes, group cuts at "
              f"{' / '.join(f'{c:.3f}' for c in cuts)}, "
              f"range {e.min():.3f}-{e.max():.3f}", flush=True)
        rows.append(m.rename("mean_lognorm").reset_index().assign(cell_type=cell_type))

    axes[0].set_ylabel("Mean expression (mean log1p CPM)")
    fig.suptitle("Expression Distribution Underlying the "
                 + {4: "Quartiles", 5: "Quintiles"}.get(args.n_groups,
                                                        f"{args.n_groups} Groups"))
    fig.tight_layout()
    save_figure(fig, args.output, dpi=150)
    pd.concat(rows, ignore_index=True).to_csv(args.output.replace(".png", ".csv"), index=False)
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
