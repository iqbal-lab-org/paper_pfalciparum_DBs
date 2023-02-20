from unittest.mock import patch
from itertools import product

from pytest import fixture
from numpy import ones as np_ones, array as np_array

from test_utils import make_alignment
from common_utils.sequences import get_pos_tuple
from seq_stats.paralog_positional_kmer_sharing import (
        KmerMatcher, PARALOG_ABBREVS,
        compute_matching_array)

PARALOG_ABBREV = "DBs"
genes = PARALOG_ABBREVS["DBs"]


@patch("Bio.AlignIO.read")
@patch("seq_stats.paralog_positional_kmer_sharing.open")
class TestGeneSampleDict:
    @classmethod
    def setup_class(cls):
        cls.samples = ["sample1", "sample2"]
        cls.all_ids = list(map(lambda elems: "_".join(elems), product(genes, cls.samples)))
        cls.all_seqs = ["AA"] * 2 + ["TT"] * 2
        cls.alignment = make_alignment(cls.all_seqs, ids=cls.all_ids)

    def test_GivenLinkedGenesAllPresent_AllLinkedGenesInOutput(
        self, patched_open, patched_AlignIO
    ):
        patched_AlignIO.return_value = self.alignment
        kmer_matcher = KmerMatcher("", PARALOG_ABBREV, [])
        result = {gene_name: kmer_matcher.get_sample_seqs(gene_name)
                for gene_name in kmer_matcher.gene_names}
        expected = {
            genes[0]: {self.samples[0]: self.all_seqs[0], self.samples[1]: self.all_seqs[1]},
            genes[1]: {self.samples[0]: self.all_seqs[2], self.samples[1]: self.all_seqs[3]},
        }
        assert result == expected

    def test_GivenNotAllLinkedGenesPresent_OnlyLinkedGenesInOutput(
        self, patched_open, patched_AlignIO
    ):
        # So now we lose 'gene2:sample2', so we expect only sample1 genes in output.
        patched_AlignIO.return_value = self.alignment[:-1]
        kmer_matcher = KmerMatcher("", PARALOG_ABBREV, [])
        result = {gene_name: kmer_matcher.get_sample_seqs(gene_name)
                for gene_name in kmer_matcher.gene_names}
        expected = {
            genes[0]: {self.samples[0]: self.all_seqs[0]},
            genes[1]: {self.samples[0]: self.all_seqs[2]},
        }
        assert result == expected

class TestComputeMatchingArray:
    def test_all_matches_returns_all_ones(self):
        pos_tuple = get_pos_tuple(kmer_size=5,alignment=["A"*10],region=None)
        num_positions = len(list(get_pos_tuple(kmer_size=5,alignment=["A"*10],region=None)))
        sequences = ["A"*10 for _ in range(5)]
        result = compute_matching_array(pos_tuple, sequences)
        expected = np_ones(num_positions)
        assert all(result == expected)

    def test_matches_and_mismatches_returns_correct_ones_and_zeros(self):
        pos_tuple = get_pos_tuple(kmer_size=2,alignment=["ATCA"],region=None)
        sequences = ["ATCA","AGCA"]
        result = compute_matching_array(pos_tuple, sequences)
        expected = np_array([0,0,1]) # The last kmer, "CA", is identical across sequences
        assert all(result == expected)
