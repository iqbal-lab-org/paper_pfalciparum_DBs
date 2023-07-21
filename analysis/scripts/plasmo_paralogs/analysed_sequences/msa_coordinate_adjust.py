from Bio import SeqIO
import click

import sys

GAP_CHAR = "-"
FIELD_DELIM = ":"

@click.command()
@click.option('--seq_id', help=
        """Fasta ID to get coordinates from""", required=True
        )
@click.option('--coords', help=
    f"""Coordinates of feature to translate. Format: 'start{FIELD_DELIM}stop', 1-based""", 
    required=True
    )
@click.option('--reverse',is_flag=True)
@click.argument('msa_in')
def main(seq_id, coords,reverse,msa_in):
    """
    In default mode, converts coordinates on the sequence with ID `seq_id` into coordinates in 
    the alignment `msa_in`.

    If flag `--reverse` is used, converts coordinates in the alignment space of `msa_in`
    into coordinates in sequence with ID `seq_id`.
    """
    used_seq = None
    for rec in SeqIO.parse(msa_in,"fasta"):
        if rec.id == seq_id:
            used_seq = str(rec.seq)
            break
    if used_seq is None:
        raise ValueError(f"{seq_id} not found in {msa_in}")
    input_start, input_end = map(int,coords.split(FIELD_DELIM))
    result_start, result_end = coordinate_adjust(input_start - 1, input_end - 1, reverse)
    print(FIELD_DELIM.join(map(str,[result_start + 1,result_end + 1])),file=sys.stdout)

def coordinate_adjust(sequence, start: int, end: int, reverse: bool = False):
    """
    Converts coordinates between (ungapped) sequence <-> (possibly gapped) msa
    Positions are 0-based.
    """
    ungapped_idx = -1
    for i, char in enumerate(sequence):
        if char != GAP_CHAR:
            ungapped_idx += 1
        if not reverse and ungapped_idx == start:
            result_start = i
        if not reverse and ungapped_idx == end:
            result_end = i
        if reverse and i == start:
            result_start = ungapped_idx
        if reverse and i == end:
            result_end = ungapped_idx
    return (result_start, result_end)

if __name__ == "__main__":
    main()
