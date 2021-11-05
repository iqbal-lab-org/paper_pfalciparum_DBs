import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fin = sys.argv[1]
fout = sys.argv[2]

with open(fin) as fhandle_in:
    total = 0
    tot_len = 0
    records = SeqIO.parse(fin,"fasta")
    dedup_records = defaultdict(list)
    for record in records:
        total += 1
        rec_string = str(record.seq).upper()
        tot_len += len(rec_string)
        dedup_records[rec_string].append(SeqRecord(record.seq,id=record.id,description=""))

avg_len = round(tot_len / total,1)
print(f"Num inputted: {total} with avg length {avg_len}, num deduplicated: {len(dedup_records)}")


with open(fout,"w") as fhandle_out:
    for rec_list in dedup_records.values():
        SeqIO.write(rec_list[0], fhandle_out, "fasta")

fout_tsv = Path(fout).with_suffix(".tsv")

with open(fout_tsv,"w") as fhandle_out:
    fhandle_out.write("Representative\tIdenticals\n")
    for rec_list in dedup_records.values():
        ids = ",".join(map(lambda rec: rec.id, rec_list[1:]))
        fhandle_out.write(f"{rec_list[0].id}\t{ids}\n")

