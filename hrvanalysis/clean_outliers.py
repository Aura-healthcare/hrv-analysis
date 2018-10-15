#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides several methods to clean abnormal and ectopic Rr Intervals."""

import pandas as pd
import numpy as np

# Static name for methods params
MALIK_RULE = "Malik"
KARLSSON_RULE = "Karlsson"
KAMATH_RULE = "Kamath"
MEAN_LAST_9 = "mean_last9"
CUSTOM_RULE = "custom"

# ----------------- ClEAN OUTlIER / ECTOPIC BEATS ----------------- #
# TODO / ONGOING


def remove_outlier(rr_intervals, print_outliers=True, low_rri=300, high_rri=2000):
    """
    Function that replace RR Interval outlier by nan

    Parameters
    ---------
    rr_intervals : list
        raw signal extracted.
    low_rri : int
        lowest RrInterval to be considered plausible.
    high_rri : int
        highest RrInterval to be considered plausible.
    print_outliers : bool
        Print information about deleted outliers.

    Returns
    ---------
    rr_intervals_cleaned : list
        list of RR Intervals without outliers

    """

    # Conversion RrInterval to Heart rate ==> rri (ms) =  1000 / (bpm / 60)
    # rri 2000 => bpm 30 / rri 300 => bpm 200
    rr_intervals_cleaned = [rri if high_rri >= rri >= low_rri else np.nan for rri in rr_intervals]

    if print_outliers:
        outliers_list = []
        for rri in rr_intervals:
            if high_rri >= rri >= low_rri:
                pass
            else:
                outliers_list.append(rri)

        nan_count = sum(np.isnan(rr_intervals_cleaned))
        print("{} outlier(s) have been deleted. The outliers values are : {}".format(nan_count,
                                                                                     outliers_list))

    return rr_intervals_cleaned


def remove_ectopic_beats(rr_intervals, method="Malik", custom_rule=None):
    """
    RR intervals differing by more than the removing_rule from the one proceeding it are removed.

    Parameters
    ---------
    rr_intervals : list
        list of Rr Intervals
    method : str
        method to use to clean outlier. Malik, Kamath, Karlsson, mean_last9 or Custom.
    custom_rule : int
        percentage criteria of difference with previous Rr Interval at which we consider
        that it is abnormal.

    Returns
    ---------
    nn_intervals : list
        list of NN Interval
    """

    # set first element in list
    nn_intervals = [rr_intervals[0]]
    outlier_count = 0
    previous_outlier = False

    if method == "Karlsson":
        for i in range(len(rr_intervals)):
            # Condition to go out of loop at limits of list
            if i == len(rr_intervals)-2:
                nn_intervals.append(rr_intervals[i + 1])
                break

            mean_prev_next_rri = (rr_intervals[i] + rr_intervals[i+2]) / 2
            if abs(mean_prev_next_rri - rr_intervals[i+1]) < 0.2 * mean_prev_next_rri:
                nn_intervals.append(rr_intervals[i+1])
            else:
                nn_intervals.append(np.nan)
                outlier_count += 1
        return nn_intervals

    if method == MEAN_LAST_9:
        nn_intervals = []
        for i, rr_interval in enumerate(rr_intervals):

            if i < 9:
                nn_intervals.append(rr_interval)
                continue

            mean_last_9_elt = np.nanmean(nn_intervals[-9:])
            if abs(mean_last_9_elt - rr_interval) < 0.3 * mean_last_9_elt:
                nn_intervals.append(rr_interval)
            else:
                nn_intervals.append(np.nan)
                outlier_count += 1

    else:
        for i, rr_interval in enumerate(rr_intervals[:-1]):

            if previous_outlier:
                nn_intervals.append(rr_intervals[i + 1])
                previous_outlier = False
                continue

            if is_outlier(rr_interval, rr_intervals[i + 1], method=method, custom_rule=custom_rule):
                nn_intervals.append(rr_intervals[i + 1])
            else:
                nn_intervals.append(np.nan)
                outlier_count += 1

    print("{} ectopic beat(s) have been deleted with {} rule.".format(outlier_count, method))

    return nn_intervals


def interpolate_nan_values(rr_intervals, interpolation_method="linear"):
    """
    Function that interpolate Nan values with linear interpolation

    Parameters
    ---------
    rr_intervals : list
        RrIntervals list.
    interpolation_method : str
        Method used to interpolate Nan values of series.

    Returns
    ---------
    interpolated_rr_intervals : list
        new list with outliers replaced by interpolated values.
    """
    series_rr_intervals_cleaned = pd.Series(rr_intervals)
    # Interpolate nan values and convert pandas object to list of values
    interpolated_rr_intervals = series_rr_intervals_cleaned.interpolate(method=interpolation_method)
    return interpolated_rr_intervals.values.tolist()


def get_nn_intervals(rr_intervals, low_rri=300, high_rri=2000, interpolation_method="linear",
                     ectopic_beats_removal_method=KAMATH_RULE, print_outliers=True):
    """
    Function that computes NN Intervals from RR Intervals.

    Parameters
    ---------
    rr_intervals : list
        RrIntervals list.
    interpolation_method : str
        Method used to interpolate Nan values of series.
    ectopic_beats_removal_method : str
        method to use to clean outlier. Malik, Kamath, Karlsson, mean_last9 or Custom.
    low_rri : int
        lowest RrInterval to be considered plausible.
    high_rri : int
        highest RrInterval to be considered plausible.
    print_outliers : bool
        Print information about deleted outliers.

    Returns
    ---------
    interpolated_nn_intervals : list
        list of NN Interval interpolated
    """
    rr_intervals_cleaned = remove_outlier(rr_intervals, low_rri=low_rri, high_rri=high_rri,
                                          print_outliers=print_outliers)
    interpolated_rr_intervals = interpolate_nan_values(rr_intervals_cleaned, interpolation_method)
    nn_intervals = remove_ectopic_beats(interpolated_rr_intervals,
                                        method=ectopic_beats_removal_method)
    interpolated_nn_intervals = interpolate_nan_values(nn_intervals, interpolation_method)
    return interpolated_nn_intervals


def is_outlier(rr_interval, next_rr_interval, method="Malik", custom_rule=None):
    """
    Test if the rr_interval is an outlier

    Parameters
    ----------
    rr_interval : int
        RrInterval
    next_rr_interval : int
        consecutive RrInterval
    method : str
        method to use to clean outlier. Malik, Kamath, Karlsson, mean_last9 or Custom
    custom_rule : int
        percentage criteria of difference with previous Rr Interval at which we consider
        that it is abnormal

    Returns
    ----------
    bool
        True if RrInterval is valid, False if not
    """
    if method not in [MALIK_RULE, KAMATH_RULE, KARLSSON_RULE, MEAN_LAST_9, CUSTOM_RULE]:
        raise ValueError("Not a valid method. Please choose Malik, Kamath or custom.")

    if method == MALIK_RULE:
        return abs(rr_interval - next_rr_interval) <= 0.2 * rr_interval
    elif method == KAMATH_RULE:
        return 0 <= (next_rr_interval - rr_interval) <= 0.325 * rr_interval or 0 <= \
               (rr_interval - next_rr_interval) <= 0.245 * rr_interval
    else:
        return abs(rr_interval - next_rr_interval) <= custom_rule * rr_interval


def is_valid_sample(nn_intervals, outlier_count, removing_rule=0.04):
    """
    Test if the sample meet the condition to be used for analysis

    Parameters
    ----------
    nn_intervals : list
        list of Normal to Normal Interval
    outlier_count : int
        count of outliers or ectopic beats removed from the interval
    removing_rule : str
        rule to follow to determine whether the sample is valid or not

    Returns
    ----------
    bool
        True if sample is valid, False if not
    """
    result = True
    if outlier_count / len(nn_intervals) > removing_rule:
        print("Too much outlier for analyses ! You should descard the sample")
        result = False
    if len(nn_intervals) < 240:
        print("Not enough Heart beat for Nyquist criteria ! ")
        result = False
    return result
