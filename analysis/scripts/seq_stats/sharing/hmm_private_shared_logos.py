"""
Separates MSA into private and shared sequences
TODO: produce logos from Skylign, using its rest API (http://skylign.org/help/#api_use)
"""
from pathlib import Path
from typing import List, Optional, Dict
from collections import defaultdict
import subprocess

import click
import sqlite3
#import requests
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from common_utils import (
    REGION_CLICK_HELP,
    VALID_SQLITE_DATA,
    VALID_SEQTYPES,
    ID_DELIM,
    REGION_DELIM,
    SHARED_FEATURE_DELIM,
)
from common_utils.sqlite import (
    get_sqlite_table_name,
    drop_if_exists,
)
from common_utils.sequences import get_pos_tuple, extract_kmer
from seq_stats.sharing import get_kmer_dict

SeqRecords = List[SeqRecord]
SplitMSA = Dict[str, SeqRecords]

GAP = "-"
UNSPEC_AA = "X"


def remove_shared_chars_from_private_seqs(split_msa: SplitMSA) -> None:
    shared_features = {
        feature for feature in split_msa if SHARED_FEATURE_DELIM in feature
    }
    shared_chars = defaultdict(set)
    for feature in shared_features:
        for rec in split_msa[feature]:
            for i, char in enumerate(rec):
                shared_chars[i].add(char)
    for feature in set(split_msa.keys()).difference(shared_features):
        for j in range(len(split_msa[feature])):
            target_rec = split_msa[feature][j]
            target_seq = str(target_rec.seq)
            new_target_seq = ""
            for pos, char in enumerate(target_seq):
                if char != GAP and char in shared_chars[pos]:
                    new_target_seq += GAP
                else:
                    new_target_seq += char
            split_msa[feature][j] = SeqRecord(
                Seq(new_target_seq), id=target_rec.id, description=""
            )


def process_all_gap_columns(split_msa: SplitMSA) -> None:
    """
    Replaces all-gap columns with unspecified amino acid ('X'), so that they appear
    in resulting hmm/logo
    """
    for feature in split_msa.keys():
        for i, rec in enumerate(split_msa[feature]):
            if set(rec) != {GAP}:
                break
        target_rec = split_msa[feature][i]
        target_seq = str(target_rec.seq)
        new_target_seq = ""
        for pos, char in enumerate(target_seq):
            chars = set([elem[pos] for elem in split_msa[feature]])
            if chars == {GAP}:
                new_target_seq += UNSPEC_AA
            else:
                new_target_seq += char
        split_msa[feature][i] = SeqRecord(
            Seq(new_target_seq), id=target_rec.id, description=""
        )


def split_msas(msa_fname, kmer_dict, out_dirname, seq_region: Optional[str]):
    all_feature_IDs = set()
    for elem in kmer_dict.values():
        for kmer in elem:
            all_feature_IDs.add(kmer.get_feature_class_ID())
    kmer_size = len(kmer.kmer)
    result: SplitMSA = {feature: list() for feature in all_feature_IDs}
    kmer_mapping = dict()
    for pos, kmers in kmer_dict.items():
        for kmer in kmers:
            key = f"{pos}{ID_DELIM}{kmer.kmer}"
            kmer_mapping[key] = kmer
    for record in AlignIO.read(open(msa_fname), "fasta"):
        sequence = str(record.seq)
        new_seqs = {feature: "" for feature in all_feature_IDs}
        first = True
        for start, center_pos, stop in get_pos_tuple(
            kmer_size, [sequence], region=seq_region
        ):
            kmer_seq = extract_kmer(sequence, start, stop)
            query = f"{center_pos}{ID_DELIM}{kmer_seq}"
            assigned_class = None
            if query in kmer_mapping:
                assigned_class = kmer_mapping[query].get_feature_class_ID()
            for feature in new_seqs:
                if feature == assigned_class:
                    to_add = kmer_mapping[query].kmer
                else:
                    to_add = GAP * kmer_size
                if not first:
                    to_add = to_add[-1]
                new_seqs[feature] += to_add
            if first:
                first = False
        assert len(set([len(elem) for elem in new_seqs.values()])) == 1
        for feature, seq in new_seqs.items():
            rec_id = f"{feature}{ID_DELIM}{record.id}"
            result[feature].append(SeqRecord(Seq(seq), id=rec_id, description=""))
    remove_shared_chars_from_private_seqs(result)
    process_all_gap_columns(result)
    return result


@click.command()
@click.argument("msa_fname")
@click.argument("out_dirname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option(
    "--kmer_sizes",
    type=str,
    help="kmer sizes to compute; comma-delimited",
    default="10",
    show_default=True,
)
@click.option("--gene_name", required=True)
@click.option("--seqtype", type=click.Choice(VALID_SEQTYPES), required=True)
@click.option("--seq_region", type=str, help=REGION_CLICK_HELP)
def main(
    msa_fname,
    out_dirname,
    sqlite_db_fpath,
    kmer_sizes,
    gene_name,
    seqtype,
    seq_region,
):
    """
    `out_dirname`: where to place the input MSA split into private/shared aligned
    sequences
    """
    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()
    input_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[1]
    )
    Path(out_dirname).mkdir(exist_ok=True, parents=True)
    for kmer_size in map(int, kmer_sizes.split(",")):
        kmer_dict = get_kmer_dict(sqlite_db_fpath, input_table_name, kmer_size)
        result = split_msas(msa_fname, kmer_dict, out_dirname, seq_region)
        if seq_region is None:
            to_add = ""
        else:
            start, end = seq_region.split(REGION_DELIM)
            to_add = f"{ID_DELIM}{start}{ID_DELIM}{end}"
        for feature, records in result.items():
            of_prefix = f"{out_dirname}/{feature}_k{kmer_size}{to_add}"
            SeqIO.write(
                records, f"{of_prefix}.msa", "fasta"
            )
            subprocess.run(
                f"hmmbuild --fragthresh 0 --enone --symfrac 0 --pnone {of_prefix}.hmm {of_prefix}.msa",
                shell=True,
                check=True,
            )


if __name__ == "__main__":
    main()
