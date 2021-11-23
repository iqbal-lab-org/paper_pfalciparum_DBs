from quality_control.compute_read_stats import mean_and_std
import numpy as np

class TestMeanStd:
    def test_empty_list(self):
        depths = []
        mean, std = mean_and_std(depths)
        assert mean == 0
        assert std == 0

    def test_non_empty_list(self):
        depths = [1,2,3,4]
        mean, std = mean_and_std(depths)
        assert mean == np.mean(depths)
        assert std == np.std(depths)
