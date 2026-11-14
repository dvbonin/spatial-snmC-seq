#!/usr/bin/env python3
"""
Covered CpGs and raw read pairs per published cell, from a manifest of BuildManifests.py.

CpGs are counted on the same convention as our own assay: strand-merged CpG
dinucleotides, a dinucleotide counting once if either of its cytosines is covered by at
least one read. A plus-strand C at position p and a minus-strand C at p+1 (1-based, as
both Bismark reports and allc files are) are the same dinucleotide. Only the nuclear
chromosomes are counted -- lambda, chrM and unplaced contigs are dropped, as in our own
pipeline. Strand-separate covered CpG cytosines are reported alongside.

Raw read pairs are taken from the SRA (spots summed over all runs of the GSM) or, for
NeMO cells, counted as the records of the cell's R1 FASTQ.

Files are streamed and discarded once counted; nothing is kept on disk but the table.
"""
import argparse
import csv
import gzip
import http.client
import io
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from functools import partial

NUCLEAR = re.compile(r"^(chr)?([0-9]+|X|Y)$")


def open_url(url):
    return urllib.request.urlopen(url, timeout=600)


def retry(fn, tries=8):
    """NCBI and NeMO drop connections mid-transfer often enough to matter, and an
    empty-body episode can last minutes -- hence the long, growing back-off."""
    for i in range(tries):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404) or i == tries - 1:   # permanent: no point waiting
                raise
            time.sleep(30 * (i + 1))
        except (OSError, http.client.HTTPException):
            if i == tries - 1:
                raise
            time.sleep(30 * (i + 1))


def spooled(url, tmp):
    """Download to a temporary file (gzip/tar need a seekable or fully streamed source)."""
    def get():
        f = tempfile.TemporaryFile(dir=tmp)
        with open_url(url) as r:
            expected = r.headers.get("Content-Length")
            shutil.copyfileobj(r, f, 1 << 22)
        # under load NCBI can answer 200 with a short or empty body; a short file would
        # otherwise be counted as a cell with low (or no) coverage
        got = f.tell()
        if got == 0 or (expected is not None and got != int(expected)):
            f.close()
            raise OSError(f"{url}: got {got} bytes, expected {expected}")
        f.seek(0)
        return f
    return retry(get)


def count_calls(lines, fmt):
    """Covered CpG cytosines and strand-merged dinucleotides from per-cytosine calls."""
    cyt = 0
    dinuc = set()
    n = 0
    for line in lines:
        n += 1
        x = line.split("\t", 6)
        chrom, pos, strand = x[0], int(x[1]), x[2]
        if not NUCLEAR.match(chrom):
            continue
        if fmt == "bismark_cpg":     # chr pos strand meth unmeth context trinuc
            if x[5] != "CG" or int(x[3]) + int(x[4]) < 1:
                continue
        else:                         # allc: chr pos strand context mc cov methylated
            if not x[3].startswith("CG") or int(x[5]) < 1:
                continue
        cyt += 1
        dinuc.add((chrom, pos if strand == "+" else pos - 1))
    if n == 0 or not dinuc:
        raise ValueError("no covered CpGs in the calls file")
    return cyt, len(dinuc)


def text_lines(fileobj):
    return io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj), encoding="ascii")


def cpgs(row, tmp):
    f = spooled(row["calls_url"], tmp)
    try:
        if row["calls_format"] == "allc_tar":
            with tarfile.open(fileobj=f, mode="r:") as t:
                m = next(m for m in t if m.name.endswith(".allc.tsv.gz"))
                return count_calls(text_lines(t.extractfile(m)), "allc")
        return count_calls(text_lines(f), row["calls_format"])
    finally:
        f.close()


def fastq_r1_records(url, tmp):
    f = spooled(url, tmp)
    try:
        with tarfile.open(fileobj=f, mode="r:") as t:
            m = next(m for m in t if m.name.endswith(".R1.fq.gz"))
            n = sum(1 for _ in gzip.GzipFile(fileobj=t.extractfile(m)))
        if n % 4:
            raise ValueError(f"{url}: R1 line count {n} is not a multiple of 4")
        return n // 4
    finally:
        f.close()


def sra_spots(gsm):
    """Spots summed over every SRA run of a GSM, with the library layout of those runs."""
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={gsm}"
    def get():
        with open_url(url) as r:
            return list(csv.DictReader(io.TextIOWrapper(r, encoding="utf-8")))
    rows = retry(get)
    rows = [x for x in rows if x.get("Run")]
    if not rows:
        raise ValueError(f"{gsm}: no SRA runs")
    layouts = {x["LibraryLayout"] for x in rows}
    return sum(int(x["spots"]) for x in rows), "/".join(sorted(layouts)), len(rows)


def process(row, tmp, cache):
    """One cell, cached as JSON so an interrupted run resumes where it stopped.
    Errors come back as plain strings: network exceptions hold open sockets, which
    cannot cross the process boundary."""
    path = os.path.join(cache, f"{row['cell']}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    try:
        out = count_cell(row, tmp)
    except urllib.error.HTTPError as e:
        if e.code not in (403, 404):
            return {"cell": row["cell"], "error": f"HTTPError: {e}"}
        # access-restricted or withdrawn upstream: not a transient failure, and not
        # something a rerun fixes -- the cell is left out and listed
        out = {"cell": row["cell"], "unavailable": f"HTTP {e.code}"}
    except Exception as e:
        return {"cell": row["cell"], "error": f"{type(e).__name__}: {e}"}
    with open(path + ".part", "w") as f:
        json.dump(out, f)
    os.replace(path + ".part", path)
    return out


def count_cell(row, tmp):
    cyt, din = cpgs(row, tmp)
    if row["reads_source"] == "sra":
        pairs, layout, runs = sra_spots(row["cell"])
    else:
        pairs, layout, runs = fastq_r1_records(row["reads_url"], tmp), "PAIRED", 1
    print(f"  {row['method']} {row['cell']}: {din:,} CpGs, {pairs:,} read pairs", flush=True)
    return dict(method=row["method"], cell=row["cell"], title=row["title"],
                covered_cytosines=cyt, covered_dinucleotides=din,
                raw_read_pairs=pairs, layout=layout, sra_runs=runs,
                cpg_per_raw_pair=round(din / pairs, 6))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--tmp", default=None)
    p.add_argument("--cache", default=None, help="per-cell results (default: <output>.cells)")
    args = p.parse_args()

    with open(args.manifest) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cache = args.cache or args.output + ".cells"
    os.makedirs(cache, exist_ok=True)
    # processes, not threads: parsing the calls is CPU-bound and would share one core
    with ProcessPoolExecutor(args.jobs) as ex:
        out = list(ex.map(partial(process, tmp=args.tmp, cache=cache), rows))

    failed = [r for r in out if "error" in r]
    if failed:
        for r in failed:
            print(f"FAILED {r['cell']}: {r['error']}", file=sys.stderr)
        sys.exit(f"{len(failed)} of {len(out)} cells failed; rerun to retry them "
                 f"(finished cells are cached in {cache})")

    gone = [r for r in out if "unavailable" in r]
    out = [r for r in out if "unavailable" not in r]
    for r in gone:
        print(f"UNAVAILABLE {r['cell']}: {r['unavailable']} -- left out", file=sys.stderr)
    single = [r for r in out if r["layout"] != "PAIRED"]
    if single:
        print(f"WARNING: {len(single)} cells are not paired-end; spots are reads, not pairs",
              file=sys.stderr)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"Wrote {args.output} ({len(out)} cells)")


if __name__ == "__main__":
    main()
