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


def ena_download_enabrowsertools(outdir, ena_ID):
    if ena_ID.startswith("ERS"):
        command = "enaGroupGet"
    elif ena_ID.startswith("ERR"):
        command = "enaDataGet"
    command += f" --format fastq {ena_ID}"
    print("Start:", command)
    subprocess.check_output(command, shell=True, cwd=outdir)
    print("Finish:", command)


def ena_download_fastq_dl(outdir, ena_ID):
    command = f"fastq-dl -a {ena_ID} --provider ENA --max-attempts 2"
    print("Start:", command)
    subprocess.check_output(command, shell=True, cwd=outdir)
    print("Finish:", command)


def make_two_read_files(outdir) -> None:
    """
    This function assumes we downloaded with `fastq-dl`, not `enaBrowserTools`.
    For how it was done using the latter, see commit c5e6467fd762d17311c55db425d0f3e4bf881f99.
    (enaBrowserTools makes one directory per sample/run accession, whereas fastq-dl puts all reads in `outdir`)
    """
    globbed_runs = Path(outdir).glob("ERR*")
    to_remove = list()
    run_ids = set()
    for path in globbed_runs:
        run_id = path.name.split("_")[0]
        run_ids.add(run_id)
        to_remove.append(path)

    for run_id in run_ids:
        assert_two_paired_fastq_files(outdir,run_id)

    for i in [1, 2]:
        reads = []
        for run_id in run_ids:
            reads.append(str(outdir / f"{run_id}_{i}.fastq.gz"))
        reads = " ".join(reads)
        if len(run_ids) > 1:
            command = f"cat {reads} > {outdir}/reads_{i}.fastq.gz"
        else:
            command = f"mv {reads} {outdir}/reads_{i}.fastq.gz"

        print("Start:", command)
        subprocess.check_output(command, shell=True)
        print("Finish:", command)
    for path in to_remove:
        path.unlink()

def assert_two_paired_fastq_files(outdir: Path, run_id: str) -> None:
    files = list(map(str,outdir.glob(f"{run_id}*")))
    if not len(files) == 2:
        raise Exception(f"Expected two files in {outdir}")
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
        raise Exception(f"Expected _1 and _2 in the two files in {outdir}")


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
        if not ena_ID.startswith("ERS") and not ena_ID.startswith("ERR"):
            ValueError(f"{ena_ID} must be of form ERR.*|ERS.*")
        ena_download_fastq_dl(outdir, ena_ID)
    make_two_read_files(outdir)


if __name__ == "__main__":
    main()
