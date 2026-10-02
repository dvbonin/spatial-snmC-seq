#!/usr/bin/env python3
"""
GroupReadsToCells.py

Split one or more tagged BAMs into per-well BAMs by an aux tag (default XP =
well position from CorrectCBCTags.py). Combines reads from all input BAMs
into a single set of per-well outputs, named `{outdir}/{tag_value}.bam`
(e.g. `A01.bam`).

With --contig-split, reads are additionally routed by alignment contig into
`{outdir}/Main/{well}.bam` (chr1-22/X/Y) and `{outdir}/Alt/{well}.bam`
(everything else -- chrM, EBV, decoy, unplaced/unlocalized scaffolds), each
class getting its own wells.tsv / read_counts.csv. This is the only stage that
already decodes every read, so classifying by RNAME here costs no extra pass.

With --threads > 1, chunk BAMs are processed in parallel, each writing to a
temporary subdirectory. Per-well BAMs are then merged with `samtools cat`.

Outputs:
  --wells-out          TSV : <well> <records_written>      (one row per observed well)
  --counts-out         CSV : WellPosition;ReadCounts
  --contig-counts-out  TSV : <contig> <records_written>    (per-contig, all wells pooled)
  With --contig-split, wells.tsv / read_counts.csv are written per class into
  `{outdir}/{class}/` instead, and --wells-out/--counts-out are ignored.
"""

import argparse
import os
import resource
import shutil
import sys
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam

# Same chr1-22/X/Y convention as ProcessCells.py's PRIMARY_CHROMS and
# PlotContigAlignmentShare.py's PRIMARY_ORDER.
PRIMARY_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
MAIN_CLASS = "Main"
ALT_CLASS = "Alt"
CONTIG_CLASSES = (MAIN_CLASS, ALT_CLASS)


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--in", dest="inbams", required=True, nargs="+",
                   help="Input BAM(s) (tagged with XP/CB/etc.)")
    p.add_argument("--outdir", required=True, help="Output directory for per-well BAMs")
    p.add_argument("--tag", default="XP", help="Aux tag used to split (e.g. XP / CB)")
    p.add_argument("--max-open", type=int, default=400,
                   help="Max simultaneously open BAM writers (one per observed tag value). "
                        "A 384-well plate needs at least 384, or 768 with --contig-split.")
    p.add_argument("--threads", type=int, default=1,
                   help="Number of chunk BAMs to process in parallel.")
    p.add_argument("--index", action="store_true", help="Run `samtools index` on each output BAM")
    p.add_argument("--require-name-sorted", action="store_true",
                   help="Enforce name-sorted input and drop pairs with inconsistent tag values")
    p.add_argument("--contig-split", action="store_true",
                   help="Also route reads by contig class into Main/ and Alt/ subdirectories")
    p.add_argument("--wells-out", default=None,
                   help="Manifest TSV: <well> <records_written> (ignored with --contig-split)")
    p.add_argument("--counts-out", default=None,
                   help="Read counts CSV: WellPosition;ReadCounts (ignored with --contig-split)")
    p.add_argument("--contig-counts-out", default=None,
                   help="Per-contig record counts TSV: <contig> <records_written>")
    return p.parse_args()


def contig_class(aln):
    return MAIN_CLASS if aln.reference_name in PRIMARY_CHROMS else ALT_CLASS


class WriterCache:
    """One pysam writer per observed key, held open for the duration of the run.

    Keys may contain a "/" (e.g. "Alt/A01" under --contig-split), in which case
    the leading component becomes a subdirectory of outdir.
    """
    def __init__(self, outdir, header, max_open=400):
        self.outdir = outdir
        self.header = header
        self.max_open = max_open
        self.writers = {}

    def _path(self, key):
        return os.path.join(self.outdir, f"{key}.bam")

    def get(self, key):
        w = self.writers.get(key)
        if w is None:
            if len(self.writers) >= self.max_open:
                raise RuntimeError(
                    f"Too many open writers (>{self.max_open}); raise --max-open"
                )
            path = self._path(key)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            w = pysam.AlignmentFile(path, "wb", header=self.header)
            self.writers[key] = w
        return w

    def close_all(self):
        for w in self.writers.values():
            w.close()
        self.writers.clear()


def get_tag(aln, tag):
    if aln.has_tag(tag):
        return str(aln.get_tag(tag))
    return None


def process_one(in_path, writers, tag, require_name_sorted, stats, split_contigs=False):
    in_bam = pysam.AlignmentFile(in_path, "rb")
    prev_qname = None
    prev_rec = None
    prev_key = None

    for rec in in_bam.fetch(until_eof=True):
        stats["total"] += 1
        stats["per_contig"][rec.reference_name] += 1
        val = get_tag(rec, tag)
        if val is None:
            stats["no_tag"] += 1
            continue
        key = f"{contig_class(rec)}/{val}" if split_contigs else val

        if require_name_sorted:
            qn = rec.query_name
            if prev_qname is None:
                prev_qname, prev_rec, prev_key = qn, rec, key
                continue
            if qn == prev_qname:
                # Both the well tag and (under --contig-split) the contig class
                # must agree, so a pair straddling a primary and a non-primary
                # contig is dropped as a unit rather than orphaned into both.
                if key == prev_key:
                    w = writers.get(key)
                    w.write(prev_rec)
                    w.write(rec)
                    stats["written"] += 2
                    stats["per_key"][key] += 2
                else:
                    stats["dropped_inconsistent_pairs"] += 2
                prev_qname = prev_rec = prev_key = None
            else:
                stats["dropped_singletons"] += 1
                prev_qname, prev_rec, prev_key = qn, rec, key
        else:
            w = writers.get(key)
            w.write(rec)
            stats["written"] += 1
            stats["per_key"][key] += 1

    if require_name_sorted and prev_qname is not None:
        stats["dropped_singletons"] += 1

    in_bam.close()


def _empty_stats():
    return {
        "total": 0,
        "no_tag": 0,
        "written": 0,
        "dropped_inconsistent_pairs": 0,
        "dropped_singletons": 0,
        "per_key": defaultdict(int),
        "per_contig": defaultdict(int),
    }


def _merge_stats(total, chunk):
    for k, v in chunk.items():
        if k in ("per_key", "per_contig"):
            for name, count in v.items():
                total[k][name] += count
        else:
            total[k] += v


def process_chunk_isolated(in_path, chunk_dir, tag, require_name_sorted, max_open,
                           split_contigs):
    """Worker: process one chunk BAM into chunk_dir, return stats dict."""
    os.makedirs(chunk_dir, exist_ok=True)
    with pysam.AlignmentFile(in_path, "rb") as f:
        header = f.header
    writers = WriterCache(chunk_dir, header, max_open=max_open)
    stats = _empty_stats()
    process_one(in_path, writers, tag, require_name_sorted, stats, split_contigs)
    writers.close_all()
    stats["per_key"] = dict(stats["per_key"])
    stats["per_contig"] = dict(stats["per_contig"])
    return stats


def _write_well_tables(keys_counts, wells_out, counts_out):
    if wells_out:
        os.makedirs(os.path.dirname(wells_out), exist_ok=True)
        with open(wells_out, "w") as out:
            out.write("well\trecords_written\n")
            for k in sorted(keys_counts):
                out.write(f"{k}\t{keys_counts[k]}\n")
    if counts_out:
        os.makedirs(os.path.dirname(counts_out), exist_ok=True)
        with open(counts_out, "w") as out:
            out.write("WellPosition;ReadCounts\n")
            for k in sorted(keys_counts):
                out.write(f"{k};{keys_counts[k]}\n")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    soft_fds = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    if args.max_open >= soft_fds - 64:
        sys.stderr.write(
            f"# WARNING: --max-open {args.max_open} is close to the open-file limit "
            f"({soft_fds}); raise it with `ulimit -n` if writers fail to open\n"
        )

    stats = _empty_stats()

    if args.threads > 1 and len(args.inbams) > 1:
        chunk_dirs = [
            os.path.join(args.outdir, f"_tmp_chunk_{i:03d}")
            for i in range(len(args.inbams))
        ]

        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(
                    process_chunk_isolated, in_path, chunk_dir,
                    args.tag, args.require_name_sorted, args.max_open,
                    args.contig_split
                ): in_path
                for in_path, chunk_dir in zip(args.inbams, chunk_dirs)
            }
            for f in as_completed(futures):
                in_path = futures[f]
                sys.stderr.write(f"# finished {in_path}\n")
                _merge_stats(stats, f.result())

        # Collect all per-well BAMs observed across chunk dirs, keyed by the
        # path relative to the chunk dir so Main/ and Alt/ stay separate.
        all_keys = set()
        for chunk_dir in chunk_dirs:
            for root, _dirs, fnames in os.walk(chunk_dir):
                for fname in fnames:
                    if fname.endswith(".bam"):
                        rel = os.path.relpath(os.path.join(root, fname), chunk_dir)
                        all_keys.add(rel[:-4])

        # Merge per-well BAMs in chunk order (preserves within-chunk name sort)
        for key in sorted(all_keys):
            parts = [
                os.path.join(d, f"{key}.bam")
                for d in chunk_dirs
                if os.path.exists(os.path.join(d, f"{key}.bam"))
            ]
            out_bam = os.path.join(args.outdir, f"{key}.bam")
            os.makedirs(os.path.dirname(out_bam), exist_ok=True)
            subprocess.run(
                ["samtools", "cat", "-o", out_bam] + parts, check=True,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )

        for chunk_dir in chunk_dirs:
            shutil.rmtree(chunk_dir, ignore_errors=True)

    else:
        with pysam.AlignmentFile(args.inbams[0], "rb") as first_in:
            header = first_in.header
        writers = WriterCache(args.outdir, header, max_open=args.max_open)
        for in_path in args.inbams:
            sys.stderr.write(f"# processing {in_path}\n")
            process_one(in_path, writers, args.tag, args.require_name_sorted, stats,
                        args.contig_split)
        writers.close_all()

    if args.contig_split:
        # One wells.tsv / read_counts.csv per contig class, so merge_wells can
        # treat Main and Alt as ordinary per-type chunk directories.
        for cls in CONTIG_CLASSES:
            prefix = f"{cls}/"
            per_well = {
                k[len(prefix):]: v
                for k, v in stats["per_key"].items() if k.startswith(prefix)
            }
            _write_well_tables(
                per_well,
                os.path.join(args.outdir, cls, "wells.tsv"),
                os.path.join(args.outdir, cls, "read_counts.csv"),
            )
    else:
        _write_well_tables(stats["per_key"], args.wells_out, args.counts_out)

    if args.contig_counts_out:
        os.makedirs(os.path.dirname(args.contig_counts_out), exist_ok=True)
        with open(args.contig_counts_out, "w") as out:
            out.write("contig\trecords\n")
            for c in sorted(stats["per_contig"], key=lambda x: (x is None, x)):
                out.write(f"{c}\t{stats['per_contig'][c]}\n")

    if args.index:
        for root, _dirs, fnames in os.walk(args.outdir):
            for fname in sorted(fnames):
                if not fname.endswith(".bam"):
                    continue
                path = os.path.join(root, fname)
                try:
                    if os.path.getsize(path) == 0:
                        continue
                    subprocess.run(["samtools", "index", path], check=False,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                except Exception:
                    pass

    sys.stderr.write(
        "# group_reads_to_wells summary\n"
        f"inputs\t{len(args.inbams)}\n"
        f"tag\t{args.tag}\n"
        f"contig_split\t{args.contig_split}\n"
        f"total_records\t{stats['total']}\n"
        f"no_tag\t{stats['no_tag']}\n"
        f"written\t{stats['written']}\n"
        f"dropped_inconsistent_pairs\t{stats['dropped_inconsistent_pairs']}\n"
        f"dropped_singletons\t{stats['dropped_singletons']}\n"
        f"keys_observed\t{len(stats['per_key'])}\n"
    )


if __name__ == "__main__":
    main()
