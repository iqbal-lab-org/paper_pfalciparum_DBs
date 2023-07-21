from typing import Iterable, Union

from Bio.Align import MultipleSeqAlignment as MSA
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from collections import defaultdict

from plasmo_paralogs.common_utils import ID_DELIM

def split_alignment_by_gene(alignment: MSA):
    result = defaultdict(lambda : MSA([]))
    for row in alignment:
        gene_name = get_gene_name(row)
        result[gene_name].append(row)
    return result

def get_gene_name(record: Union[SeqRecord,str], already_ID: bool = False):
    if already_ID:
        str_to_process = record
    else:
        str_to_process = record.id
    return str_to_process.split(ID_DELIM)[0].strip()

def get_sample_id(record: SeqRecord, already_ID: bool = False):
    if already_ID:
        str_to_process = record
    else:
        str_to_process = record.id
    return str_to_process.split(ID_DELIM)[1].strip()

def ensure_valid_alignment_slice(alignment: MSA, start: int, stop: int) -> None:
    if start < 0 or stop > len(alignment[0]):
        raise ValueError(f"start and stop must be >=0 and <= {len(alignment[0])}")

def extract_alignment_columns(alignment: MSA, start: int, stop: int) -> Iterable[str]:
    ensure_valid_alignment_slice(alignment, start, stop)
    sub_alignment = alignment[:,start: stop]
    return map(lambda rec: str(rec.seq), sub_alignment)

def ungap(seqs: Iterable[str]) -> Iterable[str]:
    for seq in seqs:
        if seq == "-":
            continue
        else:
            yield seq.replace("-","")
