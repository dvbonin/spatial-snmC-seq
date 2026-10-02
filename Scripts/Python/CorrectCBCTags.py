#!/usr/bin/env python3
import sys
from itertools import product

USAGE = """\
Usage:
  correct_cb_tags.py <mapping.csv> [--drop-ambiguous] [--drop-unmatched] [--max-distance N] [--stats-only] [--stats-out FILE]

The mapping CSV must be ';'-separated with header columns including
'Barcode' and 'WellPosition'. The 'Barcode' column is used as the
whitelist for CB-tag correction; the 'WellPosition' column is written
as an XP:Z:<well> tag for each read with a (corrected) match.

Options:
  --drop-ambiguous     Drop reads whose barcodes are within <=N mismatches of
                       multiple whitelist entries (ambiguous correction).
  --drop-unmatched     Drop reads whose barcodes exceed --max-distance from all
                       whitelist entries (no valid well assignment possible).
  --max-distance N     Maximum Hamming distance for correction (default 1).
  --stats-only         Do not output reads; only print summary stats.
  --stats-out FILE     Write statistics to FILE instead of stdout/stderr.
"""

DROP_AMB = False
DROP_UNMATCHED = False
MAX_DIST = 1
STATS_ONLY = False
STATS_PATH = None
WL_PATH = None

args = sys.argv[1:]
if not args:
    sys.stderr.write(USAGE)
    sys.exit(2)

if args[0].startswith("-"):
    sys.stderr.write(USAGE)
    sys.exit(2)
WL_PATH = args[0]

i = 1
while i < len(args):
    if args[i] == "--drop-ambiguous":
        DROP_AMB = True
    elif args[i] == "--drop-unmatched":
        DROP_UNMATCHED = True
    elif args[i] == "--stats-only":
        STATS_ONLY = True
    elif args[i] == "--max-distance":
        i += 1
        if i >= len(args):
            sys.stderr.write("ERROR: --max-distance requires a number\n")
            sys.exit(2)
        MAX_DIST = int(args[i])
    elif args[i] == "--stats-out":
        i += 1
        if i >= len(args):
            sys.stderr.write("ERROR: --stats-out requires a filename\n")
            sys.exit(2)
        STATS_PATH = args[i]
    else:
        sys.stderr.write(f"Unknown option: {args[i]}\n{USAGE}")
        sys.exit(2)
    i += 1

# ---------------- load mapping (Barcode + WellPosition) ----------------
ALPH = ("A", "C", "G", "T")
bc_to_well = {}
with open(WL_PATH) as fh:
    header = fh.readline().rstrip("\n").split(";")
    try:
        bc_col   = header.index("Barcode")
        well_col = header.index("WellPosition")
    except ValueError:
        sys.stderr.write("ERROR: mapping CSV must have 'Barcode' and 'WellPosition' columns.\n")
        sys.exit(2)
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split(";")
        bc   = f[bc_col].strip().upper()
        well = f[well_col].strip()
        if bc:
            bc_to_well[bc] = well
wlset = set(bc_to_well)
if not wlset:
    sys.stderr.write("ERROR: mapping CSV has no barcodes.\n")
    sys.exit(2)
L = len(next(iter(wlset)))
if any(len(x) != L for x in wlset):
    sys.stderr.write("ERROR: mapping barcodes must have equal length.\n")
    sys.exit(2)

# ---------------- build neighbor map ----------------
def generate_neighbors(seq, dist):
    if dist == 0:
        yield seq
        return
    if dist == 1:
        for i in range(len(seq)):
            for b in ALPH:
                if seq[i] != b:
                    yield seq[:i] + b + seq[i + 1 :]
    else:
        from itertools import combinations
        for idxs in combinations(range(len(seq)), dist):
            for repls in product(ALPH, repeat=dist):
                s = list(seq)
                valid = True
                for i, b in zip(idxs, repls):
                    if s[i] == b:
                        valid = False
                        break
                    s[i] = b
                if valid:
                    yield "".join(s)

neighbor = {}
for w in wlset:
    neighbor[w] = w  # exact self mapping
for w in wlset:
    for d in range(1, MAX_DIST + 1):
        for nb in generate_neighbors(w, d):
            if nb not in neighbor:
                neighbor[nb] = w
            else:
                if neighbor[nb] != w and neighbor[nb] != nb and neighbor[nb] is not None:
                    neighbor[nb] = None  # ambiguous

# ---------------- main loop ----------------
reads_total = 0
no_cb = 0
exact_hits = 0
one_mm_fixed = 0
ambiguous = 0
far = 0
dropped = 0

for line in sys.stdin:
    if line.startswith("@"):
        if not STATS_ONLY:
            sys.stdout.write(line)
        continue

    reads_total += 1
    f = line.rstrip("\n").split("\t")
    cb_idx = -1
    for i in range(11, len(f)):
        if f[i].startswith("CB:Z:"):
            cb_idx = i
            break
    if cb_idx < 0:
        no_cb += 1
        if not STATS_ONLY:
            sys.stdout.write(line)
        continue

    raw = f[cb_idx][5:].upper()
    corr = raw

    if raw not in neighbor:
        far += 1
        if DROP_UNMATCHED:
            dropped += 1
            continue
    else:
        hit = neighbor[raw]
        if hit is None:
            ambiguous += 1
            if DROP_AMB:
                dropped += 1
                continue  # skip
        else:
            if hit == raw:
                exact_hits += 1
            else:
                corr = hit
                one_mm_fixed += 1

    if not STATS_ONLY:
        if corr != raw:
            f[cb_idx] = "CB:Z:" + corr
            has_cr = any(x.startswith("CR:Z:") for x in f[11:])
            if not has_cr:
                f.append("CR:Z:" + raw)
        well = bc_to_well.get(corr)
        if well is not None:
            f.append("XP:Z:" + well)
        sys.stdout.write("\t".join(f) + "\n")

# ---------------- summary ----------------
def _pct(n, total):
    return f"{100*n/total:.1f}%" if total > 0 else "n/a"

assigned = exact_hits + one_mm_fixed
W = 74

def _row(stat, n, total, fate, desc):
    return (f"# {stat:<32} {n:>10,}   {_pct(n, total):>7}   {fate:<8}   {desc}")

summary = "\n".join([
    "# CorrectCBCTags summary",
    f"# {'─'*W}",
    f"# {'stat':<32} {'count':>10}   {'% total':>7}   {'fate':<8}   description",
    f"# {'─'*W}",
    f"# {'reads_total':<32} {reads_total:>10,}   {'100.0%':>7}   {'':8}   all reads entering the tagger",
    _row("no_cb_tag",                          no_cb,       reads_total, "dropped", "no CB:Z: tag found; passed through unmodified"),
    _row("exact_match",                        exact_hits,  reads_total, "kept",    "barcode matched a whitelist entry exactly"),
    _row("corrected (<="+str(MAX_DIST)+"mm)",  one_mm_fixed,reads_total, "kept",    "barcode corrected to nearest whitelist entry"),
    _row("ambiguous (<="+str(MAX_DIST)+"mm)",  ambiguous,   reads_total, "dropped", "equally close to >=2 whitelist entries"),
    _row("no_match (>"+str(MAX_DIST)+"mm)",    far,         reads_total, "dropped", "too far from all whitelist entries"),
    f"# {'─'*W}",
    _row("assigned XP tag",                    assigned,    reads_total, "kept",    "reads written with a valid well (XP:Z:) tag"),
    _row("dropped",                            dropped,     reads_total, "dropped", "reads removed from output stream"),
    f"# {'─'*W}",
]) + "\n"

if STATS_PATH:
    with open(STATS_PATH, "w") as stats:
        stats.write(summary)
else:
    sys.stderr.write(summary)
