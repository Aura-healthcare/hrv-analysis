import unittest
from hrvanalysis import *
import numpy as np


class MyTestCase(unittest.TestCase):

    def test_high_low_outlier(self):
        rri_list = [700, 600, 2300, 200, 1000, 230, 1200]
        self.assertEqual(clean_outlier(rr_intervals=rri_list), 
                         [700, 600, np.nan, np.nan, 1000, np.nan, 1200])

    def test_interpolate_cleaned_outlier(self):
        rri_list = [10, 11, np.nan, 13, 15, np.nan, np.nan, 17, 17]
        self.assertAlmostEqual(interpolate_cleaned_outlier(rri_list), (10, 11, 12, 13, 15, 15.67, 16.33, 17, 17))

    def test_1_successive_outlier_Malik(self):
        rri_list = [100, 110, 100, 130, 100, 100, 70, 100, 120, 100]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Malik"),
                         [100, 110, 100, np.nan, 100, 100, np.nan, 100, 120, 100])

    def test_1_successive_outlier_Kamath(self):
        rri_list = [101, 110, 100, 140, 100, 100, 70, 100, 130, 115, 100, 78]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Kamath"),
                         [101, 110, 100, np.nan, 100, 100, np.nan, 100, 130, 115, 100, 78])

    def test_1_successive_outlier_Karlsson(self):
        rri_list = [102, 110, 100, 125, 100, 100, 70, 100, 130, 105, 100, 78, 100]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Karlsson"),
                         [102, 110, 100, np.nan, 100, 100, np.nan, 100, np.nan, 105, 100, np.nan, 100])

    def test_1_successive_outlier_mean_last_9_elt(self):
        rri_list = [100, 100, 100, 100, 100, 100, 100, 100, 110, 930, 110, 100, 10, 100, 100]
        self.assertEqual((clean_ectopic_beats(rr_intervals=rri_list, method="mean_last9")),
                         [100, 100, 100, 100, 100, 100, 100, 100, 110, np.nan, 110, 100, np.nan, 100, 100])

    def test_2_succesive_outliers_Malik(self):
        rri_list = [103, 110, 100, 150, 50, 100, 100, 130, 115, 100, 78, 70, 100, 100]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Malik"),
                         [103, 110, 100, np.nan, np.nan, 100, 100, np.nan, 115, 100, np.nan, np.nan, 100, 110])

    def test_2_succesive_outliers_Kamath(self):
        rri_list = [104, 110, 100, 130, 115, 100, 100, 75, 140, 100, 78, 70, 92, 100, 140, 140, 100, 140, 130, 100]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Kamath"),
                         [104, 110, 100, 130, 115, 100, 100, np.nan, np.nan, 100, 78, 70, 92, 100, np.nan, np.nan,
                          100, np.nan, 130, 110])

    def test_2_succesive_outliers_Karlsson(self):
        rri_list = [105, 110, 100, 130, 115, 100, 100, 75, 140, 100, 78, 70, 92, 100, 140, 140, 100, 140, 130, 100]
        self.assertEqual(clean_ectopic_beats(rr_intervals=rri_list, method="Karlsson"),
                         [105, 110, 100, np.nan, 115, 100, 100, np.nan, np.nan, 100, np.nan, np.nan, 92, 100, np.nan,
                          np.nan, 100, np.nan, np.nan, 110])


if __name__ == '__main__':
    unittest.main()
