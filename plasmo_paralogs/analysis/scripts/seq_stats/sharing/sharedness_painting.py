import sqlite3
from pathlib import Path
from typing import Dict, List, Optional

from loguru import logger
import click
from numpy import array as np_array

from common_utils import REGION_CLICK_HELP
from common_utils.msa import AlignIO
from common_utils.metrics import PositionalMetric
from common_utils.sequences import (
    get_pos_tuple,
    extract_kmers,
)
from seq_stats.sharing import UniqueKmer, KmerDict, get_kmer_dict
from common_utils.sqlite import get_sqlite_table_name


class CustomPositionalMetricKmer(PositionalMetric):
    _headers = PositionalMetric._headers + ["kmer_size"]


class CustomPositionalMetricSampleID(CustomPositionalMetricKmer):
    _headers = CustomPositionalMetricKmer._headers + ["sample_id", "colour_mapping"]


SpreadValueDict = Dict[int, Dict[str, float]]


def either_side_yielder(start: int, max_vals: int):
    used_vals = 0
    prev_result = start
    while used_vals < max_vals:
        if used_vals % 2 == 1:
            result = prev_result + used_vals
        else:
            result = prev_result - used_vals
        yield result
        used_vals += 1
        prev_result = result


class CategoricallySpreadValueMapper:
    """
    Assigns features (here, kmers) a value between 0 and 1, where:
    - Classes get evenly spread out in the space of values
    - Classes get assigned non-overlapping value ranges
    The amount of non-overlap is determined by `CLASS_SPREAD`
    """

    CLASS_SPREAD = 1.4

    def __init__(self, max_k: int, class_names: List[str]):
        num_classes = len(class_names)
        spread_density = int(num_classes * max_k * self.CLASS_SPREAD)
        self._value_range = np_array(range(0, spread_density)) / spread_density
        self._feature_store = {class_name: list() for class_name in class_names}
        self._feature_value_positions = {}
        step = int(spread_density / num_classes)
        for i, class_name in enumerate(class_names):
            interval_start = i * step
            interval_stop = (i + 1) * step
            midpoint = int((interval_start + interval_stop) / 2)
            self._feature_value_positions[class_name] = list(
                either_side_yielder(midpoint, max_k)
            )

    def add_feature(self, class_name: str, feature):
        if class_name not in self._feature_store:
            raise ValueError(f"{class_name} not tracked")
        self._feature_store[class_name].append(feature)

    def clear_all_features(self):
        new_feature_store = {class_name: list() for class_name in self._feature_store}
        self._feature_store = new_feature_store

    def get_feature_mapping(self):
        result = dict()
        for class_name, features in self._feature_store.items():
            sorted_features = sorted(features)
            for i, feature in enumerate(sorted_features):
                value_pos = self._feature_value_positions[class_name][i]
                result[feature] = self._value_range[value_pos]
        return result


def sort_class_names(class_names: List[str]) -> List[str]:
    """
    Makes sure the combinatorial classes get middle values
    """
    sorted_class_names = list()
    to_splice = list()
    for cl in class_names:
        if UniqueKmer.FEATURE_DELIM not in cl:
            sorted_class_names.append(cl)
        else:
            to_splice.append(cl)
    midpoint = int(len(sorted_class_names) / 2)
    sorted_class_names = (
        sorted_class_names[:midpoint] + to_splice + sorted_class_names[midpoint:]
    )
    return sorted_class_names


def map_features_to_spread_values(
    kmer_dict: KmerDict, class_names: List[str]
) -> SpreadValueDict:
    result = {}
    max_k = max(map(lambda elem: len(elem), kmer_dict.values()))
    sorted_class_names = sort_class_names(class_names)
    value_mapper = CategoricallySpreadValueMapper(max_k, sorted_class_names)
    for pos, unique_kmers in kmer_dict.items():
        for kmer in unique_kmers:
            value_mapper.add_feature(kmer.get_feature_class_ID(), kmer.kmer)
        result[pos] = value_mapper.get_feature_mapping()
        value_mapper.clear_all_features()
    return result


@click.command()
@click.argument("msa_fname")
@click.argument("fout_prefix")
@click.argument("sqlite_db_fpath")
@click.option(
    "--kmer_sizes",
    type=str,
    help="Comma-delimited set of kmer sizes",
    default="5,10,15",
    show_default=True,
)
@click.option(
    "--region",
    type=str,
    help=REGION_CLICK_HELP + "By default, processes whole alignment.",
)
def main(msa_fname, fout_prefix, sqlite_db_fpath, kmer_sizes, region):
    if not Path(sqlite_db_fpath).exists():
        logger.error(f"{sqlite_db_fpath} not found")
    with open(msa_fname) as fhandle_in:
        alignment = AlignIO.read(fhandle_in, "fasta")

    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()

    reg_name = "full" if region is None else region
    fout_prefix = f"{fout_prefix}_reg{reg_name}"
    sharedness_fname = f"{fout_prefix}_sharedness.tsv"
    counts_fname = f"{fout_prefix}_counts.tsv"
    # Write file headers
    with open(counts_fname, "w") as fhandle_out:
        fhandle_out.write(CustomPositionalMetricKmer.get_header())
    with open(sharedness_fname, "w") as fhandle_out:
        fhandle_out.write(CustomPositionalMetricSampleID.get_header())
    # Populate tables
    for kmer_size in map(int, kmer_sizes.split(",")):
        table_name = get_sqlite_table_name(msa_fname, kmer_size)
        kmer_dict: KmerDict = get_kmer_dict(sqlite_db_fpath, table_name)
        all_class_names = sorted(
            [
                elem[0]
                for elem in db_cursor.execute(
                    f"select distinct feature_ID from {table_name}"
                )
            ]
        )
        value_mapper = map_features_to_spread_values(kmer_dict, all_class_names)
        compute_sharedness_at_given_k(
            kmer_dict, alignment, kmer_size, sharedness_fname, value_mapper, region
        )
        count_kmers_at_given_k(
            db_cursor, kmer_dict, alignment, kmer_size, counts_fname, region, table_name
        )
    con.close()


def compute_sharedness_at_given_k(
    kmer_dict,
    alignment,
    kmer_size: int,
    out_fname: str,
    value_mapper: SpreadValueDict,
    region: Optional[str],
) -> None:
    fhandle_out = open(out_fname, "a")
    sample_sharedness = {rec.id: dict() for rec in alignment}
    pos_tuple = get_pos_tuple(kmer_size, alignment, region)

    for start, center_pos, stop in pos_tuple:
        # Get sharedness
        for extracted_kmer in extract_kmers(alignment, start, stop, gene_id="_"):
            sample_id = list(extracted_kmer.instances)[0]
            matching_kmers = [
                elem
                for elem in kmer_dict[center_pos]
                if elem.kmer == extracted_kmer.kmer
            ]
            assert (
                len(matching_kmers) == 1
            ), f"Error: more than one sequence at given pos: {matching_kmers}"
            colour_mapping = value_mapper[center_pos][matching_kmers[0].kmer]
            SampleSharedness = CustomPositionalMetricSampleID(
                pos=center_pos,
                feature_ID=None,
                metric_ID="sharedness_class",
                metric_value=matching_kmers[0].get_feature_class_ID(),
                kmer_size=kmer_size,
                sample_id=sample_id,
                colour_mapping=colour_mapping,
            )
            fhandle_out.write(str(SampleSharedness))
    fhandle_out.close()


def count_kmers_at_given_k(
    db_cursor, kmer_dict, alignment, kmer_size: int, out_fname: str, region, table_name
) -> None:
    positional_sharedness: List[KmerSizePositionalMetric] = list()
    pos_tuple = get_pos_tuple(kmer_size, alignment, region)

    # Load feature names
    query_result = db_cursor.execute(
        f"select distinct feature_ID from {table_name}",
    ).fetchall()
    all_feature_combinations = [elem[0] for elem in query_result]

    for start, center_pos, stop in pos_tuple:
        num_kmers_per_feature = dict()
        for feature_ID in all_feature_combinations:
            num_kmers = len(
                [
                    elem
                    for elem in kmer_dict[center_pos]
                    if elem.get_feature_class_ID() == feature_ID
                ]
            )
            positional_sharedness.append(
                CustomPositionalMetricKmer(
                    pos=center_pos,
                    feature_ID=feature_ID,
                    metric_ID="unique_kmer_count",
                    metric_value=num_kmers,
                    kmer_size=kmer_size,
                )
            )

    with open(out_fname, "a") as fhandle_out:
        for elem in positional_sharedness:
            fhandle_out.write(str(elem))


if __name__ == "__main__":
    main()
