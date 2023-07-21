from unittest.mock import patch
from collections import Counter

from test_utils import make_alignment
from plasmo_paralogs.seq_stats.compute_homozygosity import compute_column_counts, compute_prob_of_same
from plasmo_paralogs.common_utils.msa import ID_DELIM

@patch("Bio.SeqIO.parse")
class TestColumnCounts:
    def test_all_seqs_identical(self, patched_SeqIO):
        self.alignment = make_alignment(
                ["AA","AA","AA"],
                [f"gene1{ID_DELIM}{i}" for i in range(1,4)]
                )
        patched_SeqIO.return_value = self.alignment
        result = compute_column_counts("")
        assert len(result) == 1
        counts = result["gene1"]
        assert counts == [Counter(A=1.0),Counter(A=1.0)]

    def test_seqs_differ(self, patched_SeqIO):
        self.alignment = make_alignment(
                ["AG","AA","CA","CA"],
                [f"gene1{ID_DELIM}{i}" for i in range(1,5)]
                )
        patched_SeqIO.return_value = self.alignment
        result = compute_column_counts("")
        assert len(result) == 1
        counts = result["gene1"]
        assert counts == [Counter(A=0.5,C=0.5),Counter(A=0.75,G=0.25)]

    def test_seqs_differ_multiple_genes(self, patched_SeqIO):
        self.alignment = make_alignment(
                ["AG","AA","CA","CA"],
                [f"gene1{ID_DELIM}{i}" for i in range(1,3)] + 
                [f"gene2{ID_DELIM}{i}" for i in range(3,5)]
                )
        patched_SeqIO.return_value = self.alignment
        result = compute_column_counts("")
        assert len(result) == 2
        assert result["gene1"] == [Counter(A=1),Counter(A=0.5,G=0.5)]
        assert result["gene2"] == [Counter(C=1),Counter(A=1)]

class TestProbOfSame:
    @classmethod
    def setup_class(cls):
        cls.counters = [Counter(A=0.4,G=0.6),Counter(A=1),Counter()]
    def test_empty_collection_return_None(self):
        assert (
                compute_prob_of_same(self.counters[0],self.counters[-1]) is
                compute_prob_of_same(self.counters[-1],self.counters[0]) is
                None)
    def test_size_one_identical_collections_return_one(self):
        result = compute_prob_of_same(self.counters[1], self.counters[1])
        assert result == 1
    def test_size_two_identical_collections_returns_correct_homozygosity(self):
        result = compute_prob_of_same(self.counters[0], self.counters[0])
        assert result == (0.4**2 + 0.6**2)
    def test_different_collections_returns_correct_homozygosity(self):
        result = compute_prob_of_same(self.counters[0], self.counters[1])
        assert result == 0.4



