import argparse
import csv
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wells-tsv",  required=True)
    ap.add_argument("--beds-dir",   required=True, help="directory containing {well}_MethStats.tsv (from WellMethStats.py)")
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--sentinel",    required=True)
    args = ap.parse_args()

    wells = []
    with open(args.wells_tsv) as f:
        next(f)
        for line in f:
            w = line.split("\t")[0].strip()
            if w:
                wells.append(w)

    # mean_beta = each site's own M/(M+U) averaged over sites (site-weighted).
    # The read-pooled ratio sum(M)/sum(M+U) is deliberately not carried: nothing
    # downstream uses it, and mixing the two estimators across contexts was what
    # made the old table internally inconsistent. n_calls is kept so per-well
    # mean coverage (n_calls / n_sites) stays derivable.
    fieldnames = ["well"]
    for ctx in ("CG", "CHG", "CHH"):
        fieldnames += [f"{ctx}_n_sites", f"{ctx}_n_calls", f"{ctx}_n_methylated", f"{ctx}_mean_beta"]

    os.makedirs(os.path.dirname(args.summary_out), exist_ok=True)
    with open(args.summary_out, "w", newline="") as out_f:
        w_out = csv.DictWriter(out_f, fieldnames=fieldnames)
        w_out.writeheader()
        for well in wells:
            stats_path = f"{args.beds_dir}/{well}_MethStats.tsv"
            row = {"well": well}
            if os.path.exists(stats_path):
                with open(stats_path) as f:
                    for r in csv.DictReader(f, delimiter="\t"):
                        ctx = r["context"]
                        if ctx not in ("CG", "CHG", "CHH"):
                            continue
                        row[f"{ctx}_n_sites"] = r["n_sites"]
                        row[f"{ctx}_n_calls"] = r["n_calls"]
                        row[f"{ctx}_n_methylated"] = r["n_methylated"]
                        row[f"{ctx}_mean_beta"] = r["mean_beta"]
            w_out.writerow(row)

    os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
    with open(args.sentinel, "w") as f:
        f.write("done\n")
    print(f"Wrote: {args.summary_out}\nSentinel: {args.sentinel}")


if __name__ == "__main__":
    main()
