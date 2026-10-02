#!/usr/bin/env python3
import argparse
import gzip
import os
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_info_field(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
    return info


def summarize_vcf(vcf_path, cx_allow=None):
    total_cov   = 0
    weighted_bt = 0.0
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            if cx_allow is not None:
                info = parse_info_field(fields[7])
                if info.get("CX") not in cx_allow:
                    continue
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")
            if len(fmt_vals) != len(fmt_keys):
                continue
            fmt = dict(zip(fmt_keys, fmt_vals))
            if "CV" not in fmt or "BT" not in fmt:
                continue
            try:
                cv = int(fmt["CV"])
                bt = float(fmt["BT"])
            except ValueError:
                continue
            if cv > 0:
                total_cov   += cv
                weighted_bt += bt * cv
    if total_cov == 0:
        return None, 0
    return weighted_bt / total_cov, total_cov


def process_dir(vcf_dir, sample, plate, cx_allow=None):
    results = {}
    for fpath in glob.glob(os.path.join(vcf_dir, "*_VCF.gz")):
        well = os.path.basename(fpath).replace("_VCF.gz", "")
        bt, cov = summarize_vcf(fpath, cx_allow=cx_allow)
        results[(f"{sample}_{plate}", well)] = (bt, cov)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lambda-dir",   required=True)
    ap.add_argument("--puc19-dir",    required=True)
    ap.add_argument("--outdir",       required=True)
    ap.add_argument("--sample",       required=True)
    ap.add_argument("--plate",        required=True)
    ap.add_argument("--summary-out",  default=None,
                    help="Explicit path for spike_qc_summary.tsv "
                         "(default: <outdir>/spike_qc_summary.tsv)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Processing lambda (ALL contexts)...")
    lambda_all = process_dir(args.lambda_dir, args.sample, args.plate, cx_allow=None)

    print("Processing pUC19 (CpG only)...")
    puc19_cpg  = process_dir(args.puc19_dir,  args.sample, args.plate, cx_allow={"CG"})

    keys = sorted(set(lambda_all.keys()) | set(puc19_cpg.keys()))

    lam_vals           = []
    puc_vals           = []
    scatter_x          = []
    scatter_y          = []
    scatter_labels     = []

    if args.summary_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.summary_out)), exist_ok=True)
        tsv_path = args.summary_out
    else:
        tsv_path = os.path.join(args.outdir, "spike_qc_summary.tsv")
    with open(tsv_path, "w") as out:
        out.write("\t".join([
            "sample", "well",
            "lambda_all_bt", "lambda_all_cov", "lambda_conversion",
            "puc19_cpg_bt", "puc19_cpg_cov",
        ]) + "\n")

        for (sample, well) in keys:
            lam_bt,  lam_cov = lambda_all.get((sample, well), (None, 0))
            puc_bt,  puc_cov = puc19_cpg.get((sample, well),  (None, 0))
            lam_conv = (1 - lam_bt) if lam_bt is not None else None

            if lam_conv is not None:
                lam_vals.append(lam_conv)
            if puc_bt is not None:
                puc_vals.append(puc_bt)
            if lam_conv is not None and puc_bt is not None:
                scatter_x.append(lam_conv)
                scatter_y.append(puc_bt)
                scatter_labels.append(well)

            out.write("\t".join([
                sample, well,
                "NA" if lam_bt   is None else f"{lam_bt:.6g}",
                str(lam_cov),
                "NA" if lam_conv is None else f"{lam_conv:.6g}",
                "NA" if puc_bt   is None else f"{puc_bt:.6g}",
                str(puc_cov),
            ]) + "\n")

    plt.figure()
    plt.hist(lam_vals, bins=50)
    plt.xlabel("Lambda conversion (ALL contexts)")
    plt.ylabel("Wells")
    plt.title("Lambda Conversion (ALL)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "lambda_conversion_all.png"), dpi=200)
    plt.close()

    plt.figure()
    plt.hist(puc_vals, bins=50)
    plt.xlabel("pUC19 methylation (CpG only)")
    plt.ylabel("Wells")
    plt.title("pUC19 CpG Methylation")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "puc19_methylation_cpg.png"), dpi=200)
    plt.close()

    fig, ax = plt.subplots()
    for x, y, label in zip(scatter_x, scatter_y, scatter_labels):
        ax.text(x, y, label, fontsize=5, ha="center", va="center", alpha=0.8)
    ax.set_xlabel("Lambda conversion (ALL)")
    ax.set_ylabel("pUC19 methylation (CpG)")
    ax.set_title("Spike QC (Lambda ALL vs pUC19 CpG)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "lambdaALL_vs_puc19CPG.png"), dpi=200)
    plt.close()

    print(f"QC complete. Wrote: {tsv_path}")


if __name__ == "__main__":
    main()
