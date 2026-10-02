#!/usr/bin/env python3
"""
Every CpG dyad in the reference, as the 0-based coordinate of its plus-strand C.

Needed to collapse the bulk plate's strands. Bismark .cov files carry no strand
column, and adjacency alone cannot settle it: in a CGCG stretch positions 2 and 3
are both CpG-context cytosines but belong to different dyads. With this set the
rule is exact -- a position p is the plus-strand C iff p is in the set, otherwise
it is the minus-strand C of the dyad starting at p-1.

Written once as an annotation rather than queried per file, so the bulk merge is
a join instead of tens of millions of FASTA lookups.

Output: chrom, start, end (= start + 2), matching the merged single-cell beds.
"""
import argparse
import gzip

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

opener = gzip.open if args.fasta.endswith(".gz") else open
n = 0
with opener(args.fasta, "rt") as fa, gzip.open(args.output, "wt") as out:
    chrom, offset, tail = None, 0, ""
    for line in fa:
        if line.startswith(">"):
            chrom, offset, tail = line[1:].split()[0], 0, ""
            continue
        # `tail` carries the last base of the previous line so a CG spanning the
        # line break is still found.
        seq = tail + line.strip().upper()
        i = seq.find("CG")
        while i != -1:
            out.write(f"{chrom}\t{offset - len(tail) + i}\t{offset - len(tail) + i + 2}\n")
            n += 1
            i = seq.find("CG", i + 1)
        offset += len(line.strip())
        tail = seq[-1:]
print(f"{n} CpG dyads written to {args.output}", flush=True)
