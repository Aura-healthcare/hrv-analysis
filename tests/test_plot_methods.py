#!/usr/bin/env python
"""This script provides methods to test extract_features methods."""

import os
import unittest
from hrvanalysis.plot import (plot_timeseries, plot_distrib, plot_poincare, plot_psd)


TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), 'test_nn_intervals.txt')


def load_test_data(path):
    # Load test rr_intervals data
    with open(path, "r") as text_file:
        lines = text_file.readlines()
    nn_intervals = list(map(lambda x: int(x.strip()), lines))
    return nn_intervals


class ExtractFeaturesTestCase(unittest.TestCase):
    """Class for UniTests of different methods in extract_features module"""

    # def test_if_plot_timseries_does_not_send_an_error(self):
    #     nn_intervals = get_rr_interval_list_from_txt_file(TEST_DATA_FILENAME)
    #     plot_timeseries(nn_intervals, normalize=True)
    #
    # def test_if_plot_distrib_does_not_send_an_error(self):
    #     nn_intervals = get_rr_interval_list_from_txt_file(TEST_DATA_FILENAME)
    #     plot_distrib(nn_intervals)
    #
    # def test_if_plot_poincare_does_not_send_an_error(self):
    #     nn_intervals = get_rr_interval_list_from_txt_file(TEST_DATA_FILENAME)
    #     plot_poincare(nn_intervals)
    #
    # def test_if_plot_psd_does_not_send_an_error(self):
    #     nn_intervals = get_rr_interval_list_from_txt_file(TEST_DATA_FILENAME)
    #     plot_psd(nn_intervals, method="lomb")


if __name__ == '__main__':
    unittest.main()
