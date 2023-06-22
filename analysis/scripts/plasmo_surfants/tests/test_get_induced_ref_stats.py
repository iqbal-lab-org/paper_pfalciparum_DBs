from plasmo_surfants.eval_varcalls.get_induced_ref_stats import (
    remove_indels,
    remove_read_ends,
    majority_pileup_is_non_ref,
    get_fraction_leq_limit,
)


class TestFractionBelowDepth:
    def test_all_zerodepth(self):
        depths = [0, 0]
        assert get_fraction_leq_limit(depths, 0) == 1
        assert get_fraction_leq_limit(depths, 1) == 1

    def test_mix_of_depths(self):
        depths = [0, 1, 2, 3]
        assert get_fraction_leq_limit(depths, 0) == 0.25
        assert get_fraction_leq_limit(depths, 2) == 0.75


class TestProcessPileup:
    def test_pileup_with_indels(self):
        pileup_string = ",.-2Aa...G"
        result = remove_indels(pileup_string)
        assert result == ",....G"

    def test_pileup_with_indels_next_to_mismatches(self):
        pileup_string = ",.-2AAAAAA+1Tt."
        result = remove_indels(pileup_string)
        assert result == ",.AAAAt."

    def test_pileup_with_indels_and_read_starts_and_ends(self):
        pileup_string = ",.-2AA.G^a$.^]AT"
        result = remove_read_ends(remove_indels(pileup_string))
        assert result == ",..G.AT"


class TestMajorityPileup:
    def test_equal_match_mismatch(self):
        pileup_string = ",.AAa,"
        assert not majority_pileup_is_non_ref(pileup_string)

    def test_majority_match(self):
        pileup_string = ",....A,"
        assert not majority_pileup_is_non_ref(pileup_string)

    def test_majority_mismatch(self):
        pileup_string = ".AAatcAA,"
        assert majority_pileup_is_non_ref(pileup_string)
