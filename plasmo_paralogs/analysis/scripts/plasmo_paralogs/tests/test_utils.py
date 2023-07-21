from typing import List
from plasmo_paralogs.common_utils.msa import MSA

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

TEST_GENE_NAME = "gene1"

def make_alignment(seqs: List[str], ids: List[str] = None) -> MSA:
    seq_lengths = set(map(len, seqs))
    assert (
        len(seq_lengths) == 1
    ), "Sequences are not the same length, does not represent an alignment"
    if ids is None:
        seqrecords = [SeqRecord(Seq(seq), id=f"{TEST_GENE_NAME}_s{i}") for i, seq in enumerate(seqs)]
    else:
        seqrecords = [SeqRecord(Seq(seq), id=ID) for seq, ID in zip(seqs, ids)]
    return MSA(seqrecords)
