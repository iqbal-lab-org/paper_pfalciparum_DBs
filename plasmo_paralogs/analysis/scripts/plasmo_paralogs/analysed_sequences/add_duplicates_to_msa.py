import sys

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

from plasmo_paralogs.common_utils.msa import get_sample_id

def usage():
    print(f"Usage: {sys.argv[0]} msa_fname dedup_tsv_fname out_fname", file = sys.stderr)
    exit(0)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    msa_fname = sys.argv[1]
    dedup_tsv_fname = sys.argv[2]
    out_fname = sys.argv[3]
    dedup_map = dict()
    with open(dedup_tsv_fname) as fhandle_in:
        for i,line in enumerate(fhandle_in):
            if i == 0:
                continue
            elems = line.split("\t")
            dedup_map[elems[0]] = elems[1].strip().split(",")
    output_records = list()
    for record in AlignIO.read(open(msa_fname),"fasta"):
        sample_id = get_sample_id(record)
        output_records.append(record)
        for duplicate_id in dedup_map[record.id]:
            if duplicate_id != "":
                output_records.append(SeqRecord(record.seq, duplicate_id, description=""))

    with open(out_fname, "w") as out_fhandle:
        SeqIO.write(output_records, out_fhandle, "fasta")
