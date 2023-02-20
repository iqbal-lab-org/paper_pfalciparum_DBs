import sqlite3
from itertools import combinations_with_replacement
from collections import Counter, defaultdict
from typing import Optional

import click
from Bio import SeqIO
import numpy as np

from common_utils import ID_DELIM, VALID_SEQTYPES, VALID_TOOLNAMES, VALID_SQLITE_DATA
from common_utils.sqlite import get_sqlite_table_name, drop_if_exists
from common_utils.msa import get_gene_name
from common_utils.metrics import MetricsRecorder


class CustomMetric(MetricsRecorder):
    _headers = ["pos", "gene_ID", "tool_name", "homozygosity"]
    _required_headers = _headers


@click.command()
@click.argument("msa_fname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option("--seqtype", required=True, type=click.Choice(VALID_SEQTYPES))
@click.option("--gene_name", required=True)
@click.option("--tool_name", required=True, type=click.Choice(VALID_TOOLNAMES))
def main(msa_fname, sqlite_db_fpath, seqtype, gene_name, tool_name):
    sub_alignments = compute_column_counts(msa_fname)

    gene_ids = list(sub_alignments.keys())
    if len(gene_ids) > 1:
        valid_gene_names = {"DBs"}
    else:
        valid_gene_names = set(gene_ids)
    if gene_name not in valid_gene_names:
        raise ValueError(
            f"{gene_name} is not a valid gene name for msa file {msa_fname}"
            f"Valid gene names: {valid_gene_names}"
        )
    gene_id_combinations = list(combinations_with_replacement(gene_ids, 2))
    result = list()
    for pos in range(len(sub_alignments[gene_ids[0]])):
        for combination in gene_id_combinations:
            if combination[0] == combination[1]:
                gene_ID = combination[0]
            else:
                gene_ID = ID_DELIM.join(combination)
            counter1 = sub_alignments[combination[0]][pos]
            counter2 = sub_alignments[combination[1]][pos]
            prob_of_same = compute_prob_of_same(counter1, counter2)
            metric_value = prob_of_same
            result.append(
                CustomMetric(
                    pos=pos,
                    gene_ID=gene_ID,
                    tool_name=tool_name,
                    homozygosity=metric_value,
                )
            )

    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()
    table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[0], tool_name=tool_name
    )
    drop_if_exists(sqlite_db_fpath, db_cursor, table_name)
    db_cursor.execute(
        f"create table {table_name} (pos,gene_ID,tool_name,homozygosity)",
    )
    for elem in result:
        db_cursor.execute(
            f"insert into {table_name} values (?,?,?,?)",
            (elem["pos"], elem["gene_ID"], elem["tool_name"], elem["homozygosity"]),
        )
    con.commit()
    con.close()


def compute_column_counts(msa_fname):
    sub_alignments = defaultdict(list)
    for record in SeqIO.parse(msa_fname, "fasta"):
        gene_name = get_gene_name(record)
        sub_alignments[gene_name].append([char for char in str(record.seq)])
    for gene_name in sub_alignments:
        array = np.array(sub_alignments[gene_name]).T
        counters = [False for _ in range(array.shape[0])]
        for idx in range(array.shape[0]):
            counters[idx] = to_relative(Counter("".join(array[idx])))
        sub_alignments[gene_name] = counters
    return sub_alignments


def to_relative(counter: Counter) -> None:
    total = sum(counter.values())
    for key in counter:
        counter[key] /= total
    return counter


def zero_one_bounded(number: float) -> bool:
    return 0 <= number <= 1


def compute_prob_of_same(collection1: Counter, collection2: Counter) -> Optional[float]:
    all_elements = set(collection1.keys()).union(set(collection2.keys()))
    # assert all(map(zero_one_bounded,
    #    chain(collection1.values(),collection2.values())
    #    )),(
    #        "Count arguments must be relative frequencies (0<=x<=1)")
    if any(map(lambda col: sum(col.values()) == 0, [collection1, collection2])):
        return None
    prob_of_same = 0
    for element in all_elements:
        prob_of_same += collection1[element] * collection2[element]
    return prob_of_same


if __name__ == "__main__":
    main()
