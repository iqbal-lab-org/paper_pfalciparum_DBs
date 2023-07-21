"""
Computes whether each kmer in each sample of an MSA is shared or not, using 
`shared_definition_by_geo` module. Also stores sample geospatial info (
country of origin, sampling date)
"""
import sqlite3

import click
from Bio import AlignIO

from plasmo_paralogs.common_utils import REGION_CLICK_HELP, VALID_SQLITE_DATA, VALID_SEQTYPES, ID_DELIM
from plasmo_paralogs.common_utils.sqlite import (
    get_sqlite_table_name,
    get_sample_geospatial,
    drop_if_exists,
)
from plasmo_paralogs.common_utils.msa import get_gene_name, get_sample_id
from plasmo_paralogs.common_utils.sequences import get_pos_tuple, extract_kmer
from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.seq_stats.sharing import get_kmer_dict, kmer_sharing_mapping


class AssignmentMetric(MetricsRecorder):
    _headers = [
        "pos",
        "sample_ID",
        "gene",
        "kmer_seq",
        "is_shared",
        "country",
        "year",
    ]
    _required_headers = _headers
    TABLE_HEADERS = "(" + ",".join(_headers) + ")"



@click.command()
@click.argument("msa_fname")
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
    sqlite_db_fpath,
    kmer_sizes,
    gene_name,
    seqtype,
    seq_region,
):
    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()

    output_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[2]
    )
    drop_if_exists(sqlite_db_fpath, db_cursor, output_table_name)
    db_cursor.execute(
        f"create table {output_table_name} {AssignmentMetric.TABLE_HEADERS}"
    )

    input_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[1]
    )
    sample_geospatial = get_sample_geospatial(sqlite_db_fpath)
    for kmer_size in map(int, kmer_sizes.split(",")):
        kmer_dict = get_kmer_dict(sqlite_db_fpath, input_table_name, kmer_size)
        kmer_mapping = kmer_sharing_mapping(kmer_dict)
        for record in AlignIO.read(open(msa_fname), "fasta"):
            sample_ID = get_sample_id(record)
            if sample_ID == "ref":
                continue
            gene_name = get_gene_name(record)
            geospatial = sample_geospatial[sample_ID]
            sequence = str(record.seq)
            for start, center_pos, stop in get_pos_tuple(
                kmer_size, [sequence], seq_region
            ):
                kmer_seq = extract_kmer(sequence, start, stop)
                kmer_mapping_key = f"{center_pos}{ID_DELIM}{kmer_seq}"
                is_shared = kmer_mapping[kmer_mapping_key]
                db_cursor.execute(
                    f"insert into {output_table_name} values (?,?,?,?,?,?,?)",
                    (
                        center_pos,
                        sample_ID,
                        gene_name,
                        kmer_seq,
                        is_shared,
                        geospatial.country,
                        geospatial.year,
                    ),
                )

    con.commit()
    con.close()


if __name__ == "__main__":
    main()
