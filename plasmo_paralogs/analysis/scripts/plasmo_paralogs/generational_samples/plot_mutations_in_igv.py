from typing import Tuple
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass

import click
import pandas as pd
from pyfaidx import Fasta
import edlib
from pysam import AlignmentFile, AlignedSegment

from plasmo_paralogs.common_utils import ID_DELIM
from plasmo_paralogs.common_utils.msa import get_gene_name, get_sample_id

DBs = {"DBLMSP", "DBLMSP2"}


def find_closest_parent(mutation_row) -> Tuple[str, str]:
    p1 = mutation_row.Parent_1_edit_distance
    p2 = mutation_row.Parent_2_edit_distance
    if p1 < p2:
        return mutation_row.Parent_1_Sample_ID, mutation_row.Parent_2_Sample_ID
    else:
        return mutation_row.Parent_2_Sample_ID, mutation_row.Parent_1_Sample_ID


def get_paralog(gene_name):
    assert gene_name in DBs
    return DBs.difference({gene_name}).pop()

class OutputSeq:
    def __init__(self, fasta_id, description):
        self.fasta_id = fasta_id
        self.description = description

    def add_seq(self,seq: str):
        self.seq = seq



@click.command()
@click.argument("seq_ifname", type=click.Path(exists=True))
@click.argument("generational_mutations_tsv", type=click.Path(exists=True))
@click.argument("out_dirname")
def main(seq_ifname, generational_mutations_tsv, out_dirname):
    seqs = Fasta(seq_ifname)
    df = pd.read_csv(generational_mutations_tsv, sep="\t")
    df = df[df["passes_all_filters"]]
    df_DBs = df[df["Gene"].isin(DBs)]

    sams = defaultdict(list)
    Path(out_dirname).mkdir(exist_ok=True, parents=True)
    for mutation_row in df_DBs.itertuples():
        target_gene = mutation_row.Gene
        closest_parent, other_parent = find_closest_parent(mutation_row)
        parent_seq_name = f"{target_gene}{ID_DELIM}{closest_parent}"
        new_dirname = Path(out_dirname) / parent_seq_name
        new_dirname.mkdir(exist_ok=True, parents=True)
        parent_seq = seqs[parent_seq_name][:].seq
        # Make fasta
        with open(f"{out_dirname}/{parent_seq_name}.fasta","w") as fout:
            fout.write(f">{parent_seq_name}\n")
            fout.write(parent_seq + "\n")

        # Add seqs to align
        seqs_to_align = list()

        seqs_to_align.append(
                OutputSeq(
                    f"{target_gene}{ID_DELIM}{mutation_row.Child_Sample_ID}",
                    "child_seq"
                    ))

        seqs_to_align.append(
            OutputSeq(
            f"{target_gene}{ID_DELIM}{other_parent}",
            "cross_partner_ortholog"
            ))

        paralog_gene = get_paralog(target_gene)
        seqs_to_align.append(OutputSeq(f"{paralog_gene}{ID_DELIM}{other_parent}","cross_partner_paralog"))

        seqs_to_align.append(OutputSeq(f"{paralog_gene}{ID_DELIM}{closest_parent}","parent_paralog"))

        # Make sams
        for elem in seqs_to_align:
            elem.add_seq(seqs[elem.fasta_id][:].seq)
            with AlignmentFile(
                f"{new_dirname}/{elem.description}.sam",
                "w",
                reference_names=[parent_seq_name],
                reference_lengths=[len(parent_seq)],
            ) as outf:
                seg = AlignedSegment(header=outf.header)
                seg.query_name = elem.fasta_id
                seg.query_sequence = elem.seq
                alignment = edlib.align(elem.seq, parent_seq, task="path")
                seg.cigarstring = alignment["cigar"]
                seg.reference_name = parent_seq_name
                seg.reference_start = 0
                seg.mapping_quality = 40
                outf.write(seg)


if __name__ == "__main__":
    main()
