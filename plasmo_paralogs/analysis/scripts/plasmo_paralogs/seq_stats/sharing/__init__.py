"""
Definition: a kmer is called `shared` if it seen across >1 gene in a multi-gene MSA (in
practice, of DBLMSP and DBLMSP2)
"""
import sqlite3
from collections import defaultdict
from typing import Set, Iterable, Dict, List

from plasmo_paralogs.common_utils import ID_DELIM, SHARED_FEATURE_DELIM

StringSet = Set[str]

# Definition of 'shared'
MIN_NUM_SAMPLES_PER_COUNTRY = 50
GEO_DEF_NAME = f"all_countries_min1_samples"


class UniqueKmer:

    def __init__(self, kmer: str, feature_classes: StringSet, instances: StringSet):
        self.kmer = kmer
        self.feature_classes = feature_classes
        self.instances = instances

    def __hash__(self):
        return hash(self.kmer)

    def __eq__(self, other):
        return self.kmer == other.kmer

    def __repr__(self):
        return f"kmer: {self.kmer}, feature_classes: {self.feature_classes}, instances: {self.instances}"

    def _add_features(self, features: StringSet):
        self.feature_classes.update(features)

    def _add_instances(self, instances: StringSet):
        for instance in instances:
            if instance in self.instances:
                raise ValueError(f"{instance} already in set of instances")
            self.instances.add(instance)

    def get_feature_class_ID(self) -> str:
        return SHARED_FEATURE_DELIM.join(sorted(list(self.feature_classes)))

    def merge_with(self, other):
        assert self.kmer == other.kmer, "Error: kmers must be identical to merge"
        self._add_features(other.feature_classes)
        self._add_instances(other.instances)

    @property
    def is_shared(self) -> bool:
        return len(self.feature_classes) > 1

    def load_features(self, feature_str: str) -> None:
        self._add_features(feature_str.split(SHARED_FEATURE_DELIM))


UniqueKmers = Iterable[UniqueKmer]
KmerDict = Dict[int, List[UniqueKmer]]


def get_kmer_dict(sqlite_db_fpath, table_name: str, kmer_size: int):
    """
    Loads a map {pos : UniqueKmer} from sqlite table
    """
    con = sqlite3.connect(sqlite_db_fpath)
    con.row_factory = sqlite3.Row
    db_cursor = con.cursor()
    result: KmerDict = defaultdict(list)
    for row in db_cursor.execute(
        f'select * from {table_name} where geo_definition = "{GEO_DEF_NAME}" and length(kmer_seq) = {kmer_size}'
    ):
        position = row["pos"]
        kmer = row["kmer_seq"]
        features = row["sharing_ID"]
        new_unique_kmer = UniqueKmer(kmer, set(), set())
        new_unique_kmer.load_features(features)
        result[position].append(new_unique_kmer)
    con.close()
    return result


def kmer_sharing_mapping(kmer_dict):
    """
    Produces map {kmer: bool} from `KmerDict`
    Bool is whether the kmer is shared in the MSA the table was built from
    """
    kmer_mapping = dict()
    for pos, kmers in kmer_dict.items():
        for kmer in kmers:
            key = f"{pos}{ID_DELIM}{kmer.kmer}"
            kmer_mapping[key] = kmer.is_shared
    return kmer_mapping

def load_kmer_sharing_mapping(table_name,con):
    result = dict()
    con.row_factory = sqlite3.Row
    db_cursor = con.cursor()
    for row in db_cursor.execute(
        f'select * from {table_name}'
    ):
        pos = row["pos"]
        kmer_seq = row["kmer_seq"]
        key = f"{pos}{ID_DELIM}{kmer_seq}"
        result[key] = bool(row["is_shared"])
    return result
