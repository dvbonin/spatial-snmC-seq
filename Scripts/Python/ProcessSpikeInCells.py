import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


def run(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\nSTDERR:\n{result.stderr}")


def process_spikein_well(well, split_dir, vcf_dir, spike_ref, tmp_base, t):
    tmp    = f"{tmp_base}/{well}"
    bam_in = f"{split_dir}/{well}.bam"
    os.makedirs(tmp, exist_ok=True)

    run(f"samtools sort -@ {t} -o {tmp}/coord.bam {bam_in}")
    run(f"samtools index -@ {t} {tmp}/coord.bam")
    run(f"biscuit pileup -@ {t} -o {tmp}/pileup.vcf {spike_ref} {tmp}/coord.bam")

    vcf_out = f"{vcf_dir}/{well}_VCF.gz"
    run(f"bgzip -c {tmp}/pileup.vcf > {vcf_out}")
    run(f"tabix -p vcf {vcf_out}")
    run(f"rm -rf {tmp}")
    return well


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wells-tsv", required=True)
    parser.add_argument("--split-dir", required=True)
    parser.add_argument("--vcf-dir",   required=True)
    parser.add_argument("--spike-ref", required=True)
    parser.add_argument("--tmp-dir",   required=True)
    parser.add_argument("--threads",   type=int, default=8)
    parser.add_argument("--sentinel",  required=True)
    args = parser.parse_args()

    wells = []
    with open(args.wells_tsv) as f:
        next(f)
        for line in f:
            w = line.split("\t")[0].strip()
            if w:
                wells.append(w)

    for d in [args.vcf_dir, args.tmp_dir]:
        os.makedirs(d, exist_ok=True)

    threads_per_well = 2
    n_parallel = max(1, args.threads // threads_per_well)
    print(f"Processing {len(wells)} spike-in wells "
          f"({n_parallel} parallel, {threads_per_well} threads each)...")

    with ThreadPoolExecutor(max_workers=n_parallel) as pool:
        futures = {
            pool.submit(
                process_spikein_well, well,
                args.split_dir, args.vcf_dir, args.spike_ref,
                args.tmp_dir, threads_per_well,
            ): well
            for well in wells
        }
        for i, fut in enumerate(as_completed(futures), 1):
            fut.result()
            if i % 10 == 0:
                print(f"  {i}/{len(wells)} done")

    os.makedirs(os.path.dirname(args.sentinel), exist_ok=True)
    with open(args.sentinel, "w") as f:
        f.write("done\n")
    print(f"Done. Sentinel: {args.sentinel}")


if __name__ == "__main__":
    main()
