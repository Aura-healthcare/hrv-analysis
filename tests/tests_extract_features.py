#!/usr/bin/env python
"""This script provides methods to test extract_features methods."""

import unittest
from hrvanalysis.extract_features import (get_time_domain_features, get_geometrical_features,
                                          create_time_info, create_interpolation_time,
                                          get_csi_cvi_features, get_poincare_plot_features,
                                          get_sampen)


class ExtractFeaturesTestCase(unittest.TestCase):
    """Class for UniTests of different methods in extract_features module"""

    def test_if_time_domain_features_are_correct(self):
        # TO DO
        nn_intervals = []
        function_time_domain_features = get_time_domain_features(nn_intervals)
        real_function_time_domain_features = []
        self.assertAlmostEqual(function_time_domain_features, real_function_time_domain_features)

    def test_if_geometrical_domain_features_are_correct(self):
        # TO DO
        nn_intervals = []
        function_geometrical_domain_features = get_geometrical_features(nn_intervals)
        real_function_geometrical_domain_features = []
        self.assertAlmostEqual(function_geometrical_domain_features,
                               real_function_geometrical_domain_features)

    def test_if_time_info_created_is_correct(self):
        # TO DO
        nn_intervals = [1000, 900, 1200, 950, 850, 1100]
        is_equal = create_time_info(nn_intervals) == [0, 0.9, 2.1, 3.05, 3.9, 5]
        self.assertTrue(is_equal)

    def test_if_interpolated_time_created_is_correct(self):
        # TO DO
        nn_intervals = []
        nni_interpolation_tmstp = create_interpolation_time(nn_intervals)
        real_interpolation_tmstp = []
        self.assertEqual(nni_interpolation_tmstp, real_interpolation_tmstp)

    def test_if_csi_cvi_features_are_correct(self):
        # TO DO
        nn_intervals = []
        function_csi_cvi_features = get_csi_cvi_features(nn_intervals)
        real_csi_cvi_features = []
        self.assertAlmostEqual(function_csi_cvi_features, real_csi_cvi_features)

    def test_if_pointcare_plot_features_features_are_correct(self):
        # TO DO
        nn_intervals = []
        function_pointcare_plot_features = get_poincare_plot_features(nn_intervals)
        real_pointcare_plot_features = []
        self.assertAlmostEqual(function_pointcare_plot_features, real_pointcare_plot_features)

    def test_if_sampen_feature_is_correct(self):
        # TO DO
        nn_intervals = []
        function_sampen_features = get_sampen(nn_intervals)
        sampen_plot_features = []
        self.assertAlmostEqual(function_sampen_features, sampen_plot_features)
