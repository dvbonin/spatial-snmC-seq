#!/usr/bin/env python3
"""
Raw read pairs per cell for sci-MET v3, as a bound rather than a measurement.

Combinatorial indexing pools every cell into one library, and neither GEO nor the
amethyst objects record reads per cell. What is public is the total: the read pairs of
every SRA run of the GEO sample (counts from ENA). Dividing by the number of cells that
made it into the object over-counts each cell's reads, since reads of barcodes that were
filtered out are in the total too -- so it is an upper bound on reads per cell, and
CpGs per read pair computed from it a lower bound on the assay's efficiency.

Writes a one-row table: total read pairs, cells, read pairs per cell, and the median
covered CpGs per cell from the metadata table.
"""
import argparse
import csv
import io
import re
import statistics
import time
import urllib.request

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=brief"
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run"
       "&fields=run_accession,read_count,library_layout,instrument_model&format=tsv")


def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except OSError:
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gsm", default="GSM8433386")   # sciMETv3 NovaSeq atlas
    p.add_argument("--metadata", required=True)     # ExtractAmethystMetadata.R output
    p.add_argument("--output", required=True)
    args = p.parse_args()

    srx = re.findall(r"SRX\d+", get(GEO.format(args.gsm)))
    runs = []
    for x in sorted(set(srx)):
        runs += list(csv.DictReader(io.StringIO(get(ENA.format(x))), delimiter="\t"))
    layouts = {r["library_layout"] for r in runs}
    if layouts != {"PAIRED"}:
        raise SystemExit(f"expected paired-end runs only, got {layouts}")
    total = sum(int(r["read_count"]) for r in runs)

    with open(args.metadata) as f:
        cg = [int(r["cg_cov"]) for r in csv.DictReader(f, delimiter="\t")]
    per_cell = total / len(cg)
    med = statistics.median(cg)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gsm", "srx", "runs", "total_read_pairs", "cells",
                    "read_pairs_per_cell_upper", "median_cg_cov", "cpg_per_raw_pair_lower"])
        w.writerow([args.gsm, len(set(srx)), len(runs), total, len(cg),
                    round(per_cell), med, round(med / per_cell, 4)])
    print(f"{args.gsm}: {len(runs)} runs, {total:,} read pairs, {len(cg):,} cells -> "
          f"<= {per_cell:,.0f} pairs/cell, >= {med / per_cell:.3f} CpGs/pair")


if __name__ == "__main__":
    main()
