#!/usr/bin/env python3
"""
Which published cells enter the assay comparison, and where their files live.

One manifest per method, one row per cell: the cell identifier, its title, the URL of
its per-cell methylation calls and the format of that file, and where its raw read
count comes from (an SRA lookup by GSM, or the cell's own R1 FASTQ).

  scBS-seq   Smallwood 2014, GSE56879. Every single cell with a CpG report. Excluded:
             the CT_* negative controls, the two *_bulk samples, and the two
             "(deeper)" libraries -- resequencing of MII_2 and MII_5, which would
             otherwise enter the comparison twice.
  snmC-seq2  Luo 2018, GSE112471 (both platform series). Every --snmc2-stride-th
             sample in GSM order, starting with the first: a systematic sample that
             spreads evenly over the pools and species of the series.
  snmC-seq3  Tian 2023 human brain atlas on NeMO. A seeded random sample of
             --snmc3-n cells among those with both an allc and a FASTQ archive,
             taken from the fetch.txt lists of the published BDBags.

sci-MET v3 has no per-cell files to fetch; its cells come from the amethyst object
(ExtractAmethystMetadata.R).
"""
import argparse
import csv
import gzip
import http.client
import io
import os
import random
import re
import tarfile
import time
import urllib.request

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
NEMO_BAGS = ("https://data.nemoarchive.org/publication_release/"
             "Ecker_Human_Methylome_2023-bag-2025-03-12")
NEMO_COUNTS = "Analysis_bag_2_Ecker_Human_Methylome_2023_Analysis_Counts"
NEMO_RAW = "Raw_data_bag_1_Ecker_Human_Methylome_2023_Raw"

FIELDS = ["method", "cell", "title", "calls_url", "calls_format", "reads_source", "reads_url"]


def fetch(url, tries=5):
    """GET with retries: NCBI drops connections mid-transfer often enough to matter."""
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.read()
        except (OSError, http.client.HTTPException):
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def series_matrix(gse, platform=None):
    """Sample columns of a GEO series matrix as a list of dicts (one per sample)."""
    stem = f"{gse}-{platform}" if platform else gse
    url = f"{GEO}/{gse[:-3]}nnn/{gse}/matrix/{stem}_series_matrix.txt.gz"
    text = gzip.decompress(fetch(url)).decode()
    cols = {}
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        key, *vals = line.split("\t")
        vals = [v.strip('"') for v in vals]
        # repeated keys (supplementary_file_1, _2 ...) are distinct; characteristics repeat
        # under the same key and are not needed here
        cols.setdefault(key[len("!Sample_"):], vals)
    n = len(cols["geo_accession"])
    return [{k: v[i] for k, v in cols.items()} for i in range(n)]


def https(url):
    return re.sub(r"^ftp://", "https://", url)


def scbs():
    rows = []
    for s in series_matrix("GSE56879"):
        t = s["title"]
        if t.startswith("CT_") or "bulk" in t or "deeper" in t:
            continue
        url = next((https(v) for k, v in s.items()
                    if k.startswith("supplementary_file") and v.endswith(".CpG.txt.gz")), None)
        if url:
            rows.append(dict(method="scBS-seq", cell=s["geo_accession"], title=t,
                             calls_url=url, calls_format="bismark_cpg",
                             reads_source="sra", reads_url=""))
    return rows


def snmc2(stride):
    samples = []
    for platform in ("GPL20301", "GPL24247"):
        samples += series_matrix("GSE112471", platform)
    samples.sort(key=lambda s: int(s["geo_accession"][3:]))
    rows = []
    for s in samples[::stride]:
        url = next(https(v) for k, v in s.items()
                   if k.startswith("supplementary_file") and "allc" in v)
        rows.append(dict(method="snmC-seq2", cell=s["geo_accession"], title=s["title"],
                         calls_url=url, calls_format="allc",
                         reads_source="sra", reads_url=""))
    return rows


def bag_fetch_list(bag):
    """URL per cell from a BDBag's fetch.txt (the bag archive holds only manifests)."""
    data = fetch(f"{NEMO_BAGS}/{bag}.tgz")
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as t:
        text = t.extractfile(f"{bag}/fetch.txt").read().decode()
    out = {}
    for line in text.splitlines():
        url = line.split("\t")[0]
        cell = os.path.basename(url).split(".")[0]
        out[cell] = url
    return out


def snmc3(n, seed):
    allc = {c: u for c, u in bag_fetch_list(NEMO_COUNTS).items() if u.endswith(".allc.tsv.tar")}
    fastq = {c: u for c, u in bag_fetch_list(NEMO_RAW).items() if u.endswith(".fastq.tar")}
    cells = sorted(set(allc) & set(fastq))
    pick = sorted(random.Random(seed).sample(cells, n))
    print(f"snmC-seq3: {len(cells):,} cells with allc and FASTQ, sampled {n} (seed {seed})")
    return [dict(method="snmC-seq3", cell=c, title=c, calls_url=allc[c],
                 calls_format="allc_tar", reads_source="fastq_tar", reads_url=fastq[c])
            for c in pick]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", required=True)
    p.add_argument("--snmc2-stride", type=int, default=4)
    p.add_argument("--snmc3-n", type=int, default=60)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for name, rows in [("scBS-seq", scbs()),
                       ("snmC-seq2", snmc2(args.snmc2_stride)),
                       ("snmC-seq3", snmc3(args.snmc3_n, args.seed))]:
        out = os.path.join(args.out_dir, f"manifest_{name}.tsv")
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, FIELDS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"{name}: {len(rows)} cells -> {out}")


if __name__ == "__main__":
    main()
