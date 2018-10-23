#!/usr/bin/env python
"""This script provides methods to test clean_outliers methods."""

import unittest
import numpy as np
from hrvanalysis.preprocessing import (remove_outlier, interpolate_nan_values,
                                       remove_ectopic_beats, get_nn_intervals)


class CleanOutliersTestCase(unittest.TestCase):
    """Class for UniTests of different methods in clean_outliers module"""

    def test_high_low_outlier(self):
        rri_list = [700, 600, 2300, 200, 1000, 230, 1200]
        self.assertAlmostEqual(remove_outlier(rr_intervals=rri_list),
                               [700, 600, np.nan, np.nan, 1000, np.nan, 1200])

    def test_interpolate_cleaned_outlier(self):
        rri_list = [10, 11, np.nan, 13, 15, np.nan, 16, 17, np.nan, np.nan, 20, np.nan]
        interpolated_list = list(interpolate_nan_values(rri_list))
        expected_list = [10, 11, 12, 13, 15, 15.5, 16, 17, 18, 19, 20, 20]
        self.assertAlmostEqual(interpolated_list, expected_list)

    def test_1_successive_outlier_malik(self):
        rri_list = [100, 110, 100, 130, 100, 100, 70, 100, 120, 100]
        self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="malik"),
                         [100, 110, 100, np.nan, 100, 100, np.nan, 100, 120, 100])

    def test_1_successive_outlier_kamath(self):
        rri_list = [101, 110, 100, 140, 100, 100, 70, 100, 130, 115, 100, 78]
        self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="kamath"),
                         [101, 110, 100, np.nan, 100, 100, np.nan, 100, 130, 115, 100, 78])

    def test_1_successive_outlier_karlsson(self):
        rri_list = [110, 100, 125, 100, 100, 70, 100, 130, 105, 100, 78, 100]
        self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="karlsson"),
                         [110, 100, np.nan, 100, 100, np.nan, 100, np.nan, 105, 100, np.nan, 100])

    def test_1_successive_outlier_Acer(self):
        rri_list = [100, 100, 100, 100, 100, 100, 100, 100, 110, 930, 110, 100, 10]
        self.assertEqual((remove_ectopic_beats(rr_intervals=rri_list, method="acar")),
                         [100, 100, 100, 100, 100, 100, 100, 100, 110, np.nan, 110, 100, np.nan])

    # def test_2_succesive_outliers_malik(self):
    #     rri_list = [103, 110, 100, 150, 50, 100, 100, 130, 115, 100, 78, 70, 100, 100]
    #     self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="Malik"),
    #                      [103, 110, 100, np.nan, np.nan, 100, 100,
    #                       np.nan, 115, 100, np.nan, np.nan, 100, 110])
    #
    # def test_2_succesive_outliers_kamath(self):
    #     rri_list = [104, 110, 100, 130, 115, 100, 100, 75, 140, 100, 78, 70, 92, 100, 140, 140,
    #                 100, 140, 130, 100]
    #     self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="Kamath"),
    #                      [104, 110, 100, 130, 115, 100, 100, np.nan, np.nan, 100, 78, 70, 92,
    #                       100, np.nan, np.nan, 100, np.nan, 130, 110])
    #
    # def test_2_succesive_outliers_karlsson(self):
    #     rri_list = [105, 110, 100, 130, 115, 100, 100, 75, 140, 100, 78, 70, 92, 100, 140, 140,
    #                 100, 140, 130, 100]
    #     self.assertEqual(remove_ectopic_beats(rr_intervals=rri_list, method="Karlsson"),
    #                      [105, 110, 100, np.nan, 115, 100, 100, np.nan, np.nan, 100, np.nan,
    #                       np.nan, 92, 100, np.nan, np.nan, 100, np.nan, np.nan, 110])

    def test_if_get_nn_intervals_creates_right_nn_intervals(self):
        rri_list = [700, 600, 2300, 1000, 1000, 230, 1200]
        expected_rri_list = [700, 600, 800, 1000, 1000, 1100, 1200]
        self.assertEqual(get_nn_intervals(rri_list), expected_rri_list)


if __name__ == '__main__':
    unittest.main()
