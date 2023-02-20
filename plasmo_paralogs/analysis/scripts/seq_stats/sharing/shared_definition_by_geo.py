"""
Computes, for each kmer in an alignment of multiple genes, whether they are shared or private,
across different geographical groups 
Used to define shared sequences for the DBs.
"""
import sqlite3
import click
from collections import defaultdict
from itertools import chain

from common_utils import (
    VALID_SEQTYPES,
    ID_DELIM,
    VALID_SQLITE_DATA,
)
from seq_stats.sharing import MIN_NUM_SAMPLES_PER_COUNTRY, GEO_DEF_NAME
from common_utils.metrics import MetricsRecorder
from common_utils.sequences import _all_kmer_positions, extract_kmers
from common_utils.sqlite import (
    get_sqlite_table_name,
    drop_if_exists,
    get_sample_geospatial,
)
from common_utils.msa import get_sample_id, split_alignment_by_gene, AlignIO


def load_sample_countries_in_alignment(alignment, sqlite_db_fpath):
    sample_geospatial = get_sample_geospatial(sqlite_db_fpath)
    samples_by_country = defaultdict(set)
    for record in alignment:
        sample_ID = get_sample_id(record)
        geospatial = sample_geospatial[sample_ID].country
        samples_by_country[geospatial].add(sample_ID)

    all_samples = set(chain.from_iterable(samples_by_country.values()))
    # Restrict to well-sampled countries
    samples_by_country = {
        key: val
        for key, val in samples_by_country.items()
        if len(val) >= MIN_NUM_SAMPLES_PER_COUNTRY
    }
    all_well_sampled_samples = set(chain.from_iterable(samples_by_country.values()))
    samples_by_country[GEO_DEF_NAME] = all_samples
    samples_by_country[f"all_countries_min{MIN_NUM_SAMPLES_PER_COUNTRY}_samples"] = all_well_sampled_samples
    return samples_by_country


class GeoSharing(MetricsRecorder):
    _headers = ["pos", "kmer_seq", "geo_definition", "seqtype", "sharing_ID"]
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
def main(msa_fname, sqlite_db_fpath, kmer_sizes, gene_name, seqtype):
    with open(msa_fname) as fhandle_in:
        alignment = AlignIO.read(fhandle_in, "fasta")

    samples_by_country = load_sample_countries_in_alignment(alignment, sqlite_db_fpath)

    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()

    sharing_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[1]
    )
    drop_if_exists(sqlite_db_fpath, db_cursor, sharing_table_name)
    db_cursor.execute(f"create table {sharing_table_name} {GeoSharing.TABLE_HEADERS}")

    shared_positional_kmer_freqs = defaultdict(int)
    shared_id = None

    for region, kept_samples in samples_by_country.items():
        filtered_alignment = [
            record for record in alignment if get_sample_id(record) in kept_samples
        ]

        sub_alignments = split_alignment_by_gene(filtered_alignment)
        for kmer_size in map(int, kmer_sizes.split(",")):
            for start, center_pos, stop in _all_kmer_positions(
                kmer_size, len(alignment[0])
            ):
                kmer_store = dict()
                for gene_id, sub_alignment in sub_alignments.items():
                    for extracted_kmer in extract_kmers(
                        sub_alignment, start, stop, gene_id
                    ):
                        if extracted_kmer.kmer in kmer_store:
                            kmer_store[extracted_kmer.kmer].merge_with(extracted_kmer)
                        else:
                            kmer_store[extracted_kmer.kmer] = extracted_kmer
                for unique_kmer in kmer_store.values():
                    db_cursor.execute(
                        f"insert into {sharing_table_name} values (?,?,?,?,?)",
                        (
                            center_pos,
                            unique_kmer.kmer,
                            region,
                            seqtype,
                            unique_kmer.get_feature_class_ID(),
                        ),
                    )
                    if unique_kmer.is_shared and not region.startswith("all"):
                        shared_positional_kmer_freqs[
                            f"{center_pos}{ID_DELIM}{unique_kmer.kmer}"
                        ] += 1
                        if shared_id is None:
                            shared_id = unique_kmer.get_feature_class_ID()
    for elem, count in shared_positional_kmer_freqs.items():
        if count >= 2:
            pos, kmer_seq = elem.split(ID_DELIM)
            db_cursor.execute(
                f"insert into {sharing_table_name} values (?,?,?,?,?)",
                (
                    int(pos),
                    kmer_seq,
                    "shared_in_min_two_countries",
                    seqtype,
                    shared_id,
                ),
            )

    country_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[3]
    )
    drop_if_exists(sqlite_db_fpath, db_cursor, country_table_name)
    db_cursor.execute(f"create table {country_table_name} (country,num_samples)")
    for region, kept_samples in samples_by_country.items():
        if region.startswith("all"):
            continue
        db_cursor.execute(
            f"insert into {country_table_name} values (?,?)",
            (region, len(kept_samples)),
        )
    con.commit()
    con.close()


if __name__ == "__main__":
    main()
