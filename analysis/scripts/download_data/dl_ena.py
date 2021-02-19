#!/usr/bin/env python3

"""
- Downloads fastqs from an ENA run or sample ID
- If from sample ID, concatenates with same pair (all *_1.fastq together, all *_2.fastq together)
- Throws if there are not exactly two read files (assumed paired-end)
"""

from glob import glob
from pathlib import Path
import shutil
import subprocess
import sys


def ena_download(outdir, ena_ID):
    if ena_ID.startswith("ERS"):
        command = "enaGroupGet"
    elif ena_ID.startswith("ERR"):
        command = "enaDataGet"
    else:
        raise ValueError(f"{ena_ID} must be of form ERR.*|ERS.*")
    command += f" --format fastq {ena_ID}"
    print("Start:", command)
    subprocess.check_output(command, shell=True, cwd=outdir)
    print("Finish:", command)


def make_two_read_files(outdir, ena_IDs) -> None:
    if ena_IDs.startswith("ERS"):
        run_ids = []
        for ena_ID in ena_IDs.split(","):
            run_ids += glob(f"{outdir}/{ena_ID}/*")
    else:
        run_ids = glob(f"{outdir}/*")

    for run_id in run_ids:
        assert_two_paired_fastq_files(run_id)

    for i in [1, 2]:
        reads = []
        for run_id in run_ids:
            reads += glob(f"{run_id}/*_{i}.fastq.gz")
        reads = " ".join(reads)
        if len(run_ids) > 1:
            command = f"cat {reads} > reads_{i}.fastq.gz"
        else:
            command = f"mv {reads} reads_{i}.fastq.gz"

        print("Start:", command)
        subprocess.check_output(command, shell=True, cwd=outdir)
        print("Finish:", command)

    for ena_ID in ena_IDs.split(","):
        shutil.rmtree(f"{outdir}/{ena_ID}")


def assert_two_paired_fastq_files(dirname: Path) -> None:
    files = glob(f"{dirname}/*")
    if not len(files) == 2:
        raise Exception(f"Expected two files in {dirname}")
    found_1 = False
    found_2 = False
    for fname in files:
        if not fname.endswith(".fastq.gz"):
            raise Exception(f"Expected a .fastq.gz file, got {fname}")
        if "_1" in fname:
            found_1 = fname
        elif "_2" in fname:
            found_2 = fname
    if not found_1 or not found_2:
        raise Exception(f"Expected _1 and _2 in the two files in {dirname}")


def main():
    try:
        ena_IDs = sys.argv[1]
        outdir = Path(sys.argv[2]).resolve()
    except:
        print(f"Usage: {sys.argv[0]} ena_ID output_dir \n"
                "Each ena_ID looks like ERR*/ERS*; multiple ena_IDs can be passed as comma-separated")
        exit(1)

    outdir.mkdir(parents=True, exist_ok=True)
    exist = [Path(f"{outdir}/reads_{i}.fastq.gz").exists() for i in [1, 2]]
    if sum(exist) == 2:
        print(f"{outdir} already has reads, nothing to do", file=sys.stderr)
        return

    for ena_ID in ena_IDs.split(","):
        ena_download(outdir, ena_ID)
    make_two_read_files(outdir, ena_IDs)


if __name__ == "__main__":
    main()
