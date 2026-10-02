#!/usr/bin/env python3
"""
One summary table per FirstPrep plate: well position, cell count, raw and
deduplicated read counts, CpG coverage, and the QC verdict.

Reads joined on WellPosition from three sources:
  well_stats.csv        Analysis.smk    cell count, ReadCounts_Main, CpGsCovered
  *_raw_reads.csv       CountRawBarcodes / MergeRawBarcodeCounts
  qc_pass.csv           Preprint/QC/ApplyQCToWellStats.R

The QC verdict is taken from qc_pass.csv rather than from well_stats' own
passes_cpg_qc column: a rerun of Analysis.smk's build_well_stats resets that
column to NA, while qc_pass.csv is the durable record of the decision.
"""
import argparse
import csv
import os

ap = argparse.ArgumentParser()
ap.add_argument("--well-stats-dir", required=True)
ap.add_argument("--raw-reads-dir", required=True)
ap.add_argument("--qc-pass", required=True)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--plates", default="P3,P4")
args = ap.parse_args()

qc = {}
with open(args.qc_pass) as f:
    for r in csv.DictReader(f):
        qc[(r["Plate"], r["WellPosition"])] = r["passes_cpg_qc"]

os.makedirs(args.out_dir, exist_ok=True)
for plate in args.plates.split(","):
    with open(os.path.join(args.raw_reads_dir, f"FirstPrep_{plate}_raw_reads.csv")) as f:
        raw = {r["WellPosition"]: r["RawReads"] for r in csv.DictReader(f)}

    out = os.path.join(args.out_dir, f"FirstPrep_{plate}_well_summary.csv")
    n = passed = 0
    with open(os.path.join(args.well_stats_dir, f"FirstPrep_{plate}_well_stats.csv")) as f, \
         open(out, "w", newline="") as g:
        w = csv.writer(g, lineterminator="\n")
        w.writerow(["WellPosition", "Cells", "RawReads", "Reads", "CpGs", "PassesQC"])
        for r in csv.DictReader(f):
            well = r["WellPosition"]
            verdict = qc.get((plate, well), "")
            w.writerow([well, r["CellCount"], raw[well],
                        r["ReadCounts_Main"], r["CpGsCovered"], verdict])
            n += 1
            passed += verdict == "TRUE"
    print(f"{plate}: {n} wells, {passed} pass QC -> {out}")
