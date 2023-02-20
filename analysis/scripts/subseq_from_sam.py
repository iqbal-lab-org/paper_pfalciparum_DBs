from pysam import AlignmentFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import sys
import gzip
from pathlib import Path
from itertools import chain

def usage():
    print(f"""
           usage: {sys.argv[0]} alignment_dir fasta_dir gene1,gene2,...,

           For each [sam/bam] file in `alignment_dir`, if a fasta genome with the same filename
           prefix exists, will extract gene sequences you asked for using their mapped
           coordinates in the sam file.
           Outputs to stdout.
           """)
    exit(0)

def total_stem(fname: Path):
    result = fname
    while "." in result.name:
        result = Path(result.stem)
    return result.name

def load_fasta(fname: Path):
    if fname.suffix.endswith("gz"):
        fhandle = gzip.open(fname,"rt")
    else:
        fhandle = open(fname,"r")
    result = dict()
    records = SeqIO.parse(fhandle, "fasta")
    for record in records:
        result[record.id] = str(record.seq)
    fhandle.close()
    return result

if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()
    alignment_dir = Path(sys.argv[1])
    fasta_dir = Path(sys.argv[2])
    gene_set = set(sys.argv[3].split(","))
    if not alignment_dir.exists() or not fasta_dir.exists():
        raise ValueError(f"{alignment_dir} or {fasta_dir} not found")

    align_files = alignment_dir.glob("*.[sb]am")
    fasta_files = chain.from_iterable([fasta_dir.glob("*fasta"),
            fasta_dir.glob("*fasta.gz")])
    available_fastas = {total_stem(fname): fname for fname in fasta_files}

    for afname in align_files:
        sample_name = afname.stem
        if sample_name not in available_fastas:
            print(f"{afname} has no sample match in {alignment_dir}",file=sys.stderr)
            continue
        chroms = load_fasta(available_fastas[sample_name])
        mode = "r"
        if afname.suffix == "bam":
            mode += "b"
        af = AlignmentFile(afname, mode)
        for record in af.fetch():
            if record.query_name in gene_set:
                start = record.reference_start 
                if start < 0:
                    print(f"{record.query_name} in {afname} is unaligned",file=sys.stderr)
                    continue
                end = start + len(record.query_sequence)
                chrom = record.reference_name
                matching_string = chroms[chrom][start:end]
                output_id = f"{record.query_name}@{chrom}:{start}-{end}"
                print(f">{output_id}\n{matching_string}",file=sys.stdout)
