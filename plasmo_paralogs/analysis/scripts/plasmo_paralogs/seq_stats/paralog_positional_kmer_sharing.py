"""
Computes how often aligned kmers are identical/different across paralogs in a 
given population of sequences. Computes this both empirically- looking at 
actual aligned paralog pairs on a genome- and randomly, by randomly matching paralog
pairs.

Supports looking at specific sequence regions/geographical regions only.

FPMK: the fraction of paralogs (on same genome) with identical aligned kmers, at each 
position, across the population. This can be computed directly from the alignment.

FPPK: same as FPMK but for private kmers. This can only be computed from a prior 
definition of shared/private kmers in the population, as a mismatch between paralogs
in same genome does not imply the kmers are not shared globally.
"""
import sqlite3
from typing import Set, List, Iterable, Tuple
import re
from itertools import repeat, starmap
from collections import defaultdict
from random import choice as rand_choice

from loguru import logger
from numpy import zeros as np_zeros, array as np_array
import click
from Bio import AlignIO

from plasmo_paralogs.seq_stats.sharing import get_kmer_dict, kmer_sharing_mapping
from plasmo_paralogs.common_utils import REGION_CLICK_HELP, ID_DELIM, VALID_SEQTYPES, VALID_SQLITE_DATA
from plasmo_paralogs.common_utils.msa import get_gene_name, get_sample_id
from plasmo_paralogs.common_utils.metrics import MetricsRecorder
from plasmo_paralogs.common_utils.sequences import get_pos_tuple, extract_kmer
from plasmo_paralogs.common_utils.sqlite import (
    get_sqlite_table_name,
    get_sample_ids_matching_geo_region,
    drop_if_exists
)

PARALOG_ABBREVS = {"DBs": ("DBLMSP", "DBLMSP2")}


class Stats(MetricsRecorder):
    """
    FPMK: Fraction Positionally matching kmers
    FPPK: Fraction Positionally private kmers
    """

    _headers = [
        "pos",
        "FPMK",
        "FPPK",
        "kmer_size",
        "model",
        "seq_region",
        "geo_region",
    ]
    _required_headers = _headers
    TABLE_HEADERS = "(" + ",".join(_headers) + ")"

    def get_elems(self):
        return [self[key] for key in self._headers]


@click.command()
@click.argument("msa_fname")
@click.option("-s", "--sqlite_db_fpath", required=True)
@click.option(
    "--kmer_sizes",
    type=str,
    help="Comma-delimited set of kmer sizes",
    default="10",
    show_default=True,
)
@click.option("--seqtype", type=click.Choice(VALID_SEQTYPES), required=True)
@click.option(
    "--seq_region",
    type=str,
    help=REGION_CLICK_HELP + "By default, processes whole alignment.",
)
@click.option(
    "--geo_region",
    type=str,
    help="Regexp for geographical site or country to restrict sample analysis to",
)
@click.option(
    "--gene_name",
    type=str,
    help=f"The gene set to analyse. Supported are: {PARALOG_ABBREVS}",
    default=list(PARALOG_ABBREVS.keys())[0],
)
def main(
    msa_fname,
    sqlite_db_fpath,
    kmer_sizes,
    seqtype,
    seq_region,
    geo_region,
    gene_name,
):
    matching_ids = get_sample_ids_matching_geo_region(sqlite_db_fpath, geo_region)
    kmer_matcher = KmerMatcher(msa_fname, gene_name, matching_ids)
    result = list()
    reg_name = "full" if seq_region is None else seq_region
    geo_name = (
        "all"
        if geo_region is None
        else geo_region.replace("|", "_").replace(" ", ID_DELIM)
    )
    all_positional_matches = list()
    all_gene_fraction_matches = list()
    for kmer_size in map(int, kmer_sizes.split(",")):
        input_table_name = get_sqlite_table_name(
            msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[1]
        )
        kmer_dict = get_kmer_dict(sqlite_db_fpath, input_table_name, kmer_size)
        kmer_matcher.kmer_mapping = kmer_sharing_mapping(kmer_dict)
        for model_name, pairing_function in zip(
            ["empirical_linkage", "random_linkage"],
            [kmer_matcher.linked_pairing, kmer_matcher.random_pairing],
        ):
            positional_matches = kmer_matcher.compute_positional_gene_stats(
                kmer_size,
                seq_region,
                pairing_function,
            )
            for pos, match_fraction, private_fraction in positional_matches:
                all_positional_matches.append(
                    Stats(
                        pos=pos,
                        FPMK=match_fraction,
                        FPPK=private_fraction,
                        kmer_size=kmer_size,
                        model=model_name,
                        seq_region=reg_name,
                        geo_region=geo_name,
                    )
                )

    con = sqlite3.connect(sqlite_db_fpath)
    db_cursor = con.cursor()
    output_table_name = get_sqlite_table_name(
        msa_fname, gene_name, seqtype, metric_info=VALID_SQLITE_DATA[4]
    )
    drop_if_exists(sqlite_db_fpath, db_cursor, output_table_name)
    db_cursor.execute(
        f"create table {output_table_name} {Stats.TABLE_HEADERS}"
    )
    for stats in all_positional_matches:
        db_cursor.execute(
            f"insert into {output_table_name} values (?,?,?,?,?,?,?)",
            stats.get_elems()
            )
    con.commit()
    con.close()


class KmerMatcher:
    def __init__(self, msa_fname: str, gene_name: str, matching_ids: Set[str]):
        self.gene_names = PARALOG_ABBREVS[gene_name]
        for gene_name in self.gene_names:
            setattr(self, gene_name, dict())
        self._init_sample_dict(msa_fname, matching_ids)

    @property
    def first_gene(self):
        return self.gene_names[0]

    @property
    def first_sample_seqs(self):
        return self.get_sample_seqs(self.first_gene)

    @property
    def first_sequence(self):
        return next(iter(self.first_sample_seqs.values()))

    def get_sample_seqs(self, gene_name: str):
        return getattr(self, gene_name)

    def _init_sample_dict(self, msa_fname: str, matching_ids: Set[str]):
        """
        Produces a mapping {sample_name: sequence} such that the same sample
        names appear in each gene_name entry, and sample_name filtered to `matching_ids`
        """
        sample_dict = defaultdict(list)
        # TODO: log skipping non-matching ids (due to filtering by geography)
        for record in AlignIO.read(open(msa_fname), "fasta"):
            sample_id = get_sample_id(record)
            if len(matching_ids) > 0 and sample_id not in matching_ids:
                continue
            gene_name = get_gene_name(record)
            if gene_name not in self.gene_names:
                continue
            sample_dict[sample_id].append(record)

        filtered_singletons = 0
        for sample_id, records in sample_dict.items():
            if len(records) == 1:
                filtered_singletons += 1
            elif len(records) == 2:
                for record in records:
                    gene_name = get_gene_name(record)
                    self.get_sample_seqs(gene_name)[sample_id] = str(record.seq)
            else:
                raise ValueError(f"More than two paralogs for {sample_id}")
        all_lens = list(
            map(lambda elem: len(self.get_sample_seqs(elem)), self.gene_names)
        )
        assert (
            len(set(all_lens)) == 1
        ), f"Error: not all sample sequences have each of the paralog represented"
        logger.info(
            f"Filtered out {filtered_singletons} singletons (gene sequence with no paralog in sample genome)"
        )
        logger.info(f"Have {all_lens[0]} linked genes")

    def random_pairing(self, gene_name, sid) -> Tuple[str, str]:
        sample_seqs = self.get_sample_seqs(gene_name)
        key_pool = list(sample_seqs.keys())
        selected_matched_sid = rand_choice(key_pool)
        return (selected_matched_sid, sample_seqs[selected_matched_sid])

    def linked_pairing(self, gene_name, sid) -> Tuple[str, str]:
        return (sid, self.get_sample_seqs(gene_name)[sid])

    def compute_positional_gene_stats(
        self, kmer_size, seq_region, pairing_function
    ) -> Iterable[Tuple[int, float, float]]:
        """
        Computes match statistics per aligned paralog position across all samples
        """
        pseudo_alignment = [self.first_sequence]
        all_positions = list(get_pos_tuple(kmer_size, pseudo_alignment, seq_region))
        match_array = np_zeros(len(all_positions))
        no_sharing_array = np_zeros(len(all_positions))
        num_sids = 0
        for sid, seq in self.first_sample_seqs.items():
            num_sids += 1
            linked_genes = [seq]
            for gene_name in self.gene_names[1:]:
                _, selected_seq = pairing_function(gene_name, sid)
                linked_genes.append(selected_seq)
            match_array += compute_matching_array(all_positions, linked_genes)
            no_sharing_array += self.compute_no_sharing_array(
                all_positions, linked_genes
            )
        match_array /= num_sids
        no_sharing_array /= num_sids
        assert len(all_positions) == len(match_array) == len(no_sharing_array)
        return zip(
            map(lambda elem: elem[1], all_positions), match_array, no_sharing_array
        )

    def compute_no_sharing_array(self, pos_tuple, seqlist: List[str]):
        pre_result = list()
        for start, center_pos, stop in pos_tuple:
            kmers = starmap(extract_kmer, zip(seqlist, repeat(start), repeat(stop)))
            keys = [f"{center_pos}{ID_DELIM}{kmer}" for kmer in kmers]
            shared_kmers = [self.kmer_mapping[key] for key in keys]
            if sum(shared_kmers) == 0:
                pre_result.append(1)
            else:
                pre_result.append(0)
        return np_array(pre_result)


def compute_matching_array(pos_tuple, seqlist: List[str]):
    """
    Computes whether all kmers in `seqlist` are identical at positions given by
    `pos_tuple`
    """
    pre_result = list()
    for start, _, stop in pos_tuple:
        kmers = starmap(extract_kmer, zip(seqlist, repeat(start), repeat(stop)))
        if len(set(kmers)) == 1:
            pre_result.append(1)
        else:
            pre_result.append(0)
    return np_array(pre_result)


if __name__ == "__main__":
    main()
