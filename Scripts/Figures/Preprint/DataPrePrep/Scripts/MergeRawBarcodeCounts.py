#!/usr/bin/env python3
"""
Per-chunk raw barcode tallies -> raw reads per well, one CSV per plate.

Barcodes are corrected exactly as the pipeline does it (CorrectCBCTags.py with
--max-distance 1 --drop-ambiguous --drop-unmatched), so the result is comparable
to the ReadCounts_* columns rather than an exact-match undercount. All 1-mismatch
variants of the 384 whitelist barcodes are precomputed once (384 x 8 x 4 keys);
a variant reachable from two different wells is ambiguous and dropped.

Counts are read PAIRS (R1 only). RawReads = 2 x RawReadPairs, since every raw
pair has both mates before any filtering.
"""
import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--counts-dir", required=True)
ap.add_argument("--mapping", required=True)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()

wl = {}
with open(args.mapping) as f:
    for row in csv.DictReader(f, delimiter=";"):
        wl[row["Barcode"]] = row["WellPosition"]
L = len(next(iter(wl)))

# exact matches win; 1-mismatch variants resolve only when unambiguous
variants = {}
for bc, well in wl.items():
    for i in range(L):
        for alt in "ACGTN":
            if alt == bc[i]:
                continue
            v = bc[:i] + alt + bc[i + 1:]
            if v in wl:
                continue
            variants[v] = None if (v in variants and variants[v] != well) else well
lookup = dict(wl)
lookup.update({v: w for v, w in variants.items() if w is not None})
ambiguous = {v for v, w in variants.items() if w is None}

per_plate = defaultdict(Counter)
stats = defaultdict(Counter)
for path in sorted(glob.glob(os.path.join(args.counts_dir, "*.tsv"))):
    plate = re.search(r"_(P\d+)_C\d+$", os.path.basename(path)[:-4]).group(1)
    with open(path) as f:
        for line in f:
            bc, n = line.split("\t")
            n = int(n)
            stats[plate]["total"] += n
            if bc in wl:
                per_plate[plate][lookup[bc]] += n
                stats[plate]["exact"] += n
            elif bc in lookup:
                per_plate[plate][lookup[bc]] += n
                stats[plate]["corrected"] += n
            elif bc in ambiguous:
                stats[plate]["ambiguous"] += n
            else:
                stats[plate]["unmatched"] += n

os.makedirs(args.out_dir, exist_ok=True)
for plate, counts in sorted(per_plate.items()):
    out = os.path.join(args.out_dir, f"FirstPrep_{plate}_raw_reads.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["WellPosition", "RawReadPairs", "RawReads"])
        for well in sorted(set(wl.values())):
            p = counts.get(well, 0)
            w.writerow([well, p, 2 * p])
    s = stats[plate]
    print(f"{plate}: {s['total']/1e6:.1f}M pairs seen  "
          f"exact={s['exact']/1e6:.1f}M  corrected={s['corrected']/1e6:.1f}M  "
          f"ambiguous={s['ambiguous']/1e6:.1f}M  unmatched={s['unmatched']/1e6:.1f}M  "
          f"-> assigned={(s['exact']+s['corrected'])/1e6:.1f}M "
          f"({100*(s['exact']+s['corrected'])/s['total']:.1f}%)")
    print(f"   wrote {out}")
