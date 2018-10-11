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


def clean_outlier(rr_intervals, print_outliers=True, low_rri=300, high_rri=2000):
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


def interpolate_nan_values(rr_intervals, method="linear"):
    """
    Function that interpolate Nan values with linear interpolation

    Parameters
    ---------
    rr_intervals : list
        RrIntervals list.
    method : str
        Method used to interpolate Nan values of series.

    Returns
    ---------
    interpolated_rr_intervals : list
        new list with outliers replaced by interpolated values.
    """
    series_rr_intervals_cleaned = pd.Series(rr_intervals)
    # Interpolate nan values and convert pandas object to list of values
    interpolated_rr_intervals = series_rr_intervals_cleaned.interpolate(method=method)
    return interpolated_rr_intervals.values.tolist()


def clean_ectopic_beats(rr_intervals, method="Malik", custom_rule=None):
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

    # if method == "Karlsson":
    #     if i == len(rr_intervals)-2:
    #         break
    #     mean_prev_next_rri = (rr_interval + rr_intervals[i + 2]) / 2
    #     if abs(mean_prev_next_rri - rr_intervals[i+1]) < 0.2 * mean_prev_next_rri:
    #         nn_intervals.append(rr_intervals[i+1])
    #     else:
    #         nn_intervals.append(np.nan)
    #         outlier_count += 1
    #         previous_outlier = True
    #
    # return nn_intervals

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
                # previous_outlier = True
    else:
        for i, rr_interval in enumerate(rr_intervals[:-1]):

            if previous_outlier:
                nn_intervals.append(rr_intervals[i + 1])
                previous_outlier = False
                continue

            # TO DO pour v2 ... Check si plusieurs outliers consécutifs. Quelle règle appliquer ?
            # while previous_outlier:
            #   j += 1
            #  if is_outlier(rr_interval, rr_intervals[i+1+j]):
            #     nn_intervals.append(np.nan)
            #    continue
            # else:
            #    previous_outlier = False

            if is_outlier(rr_interval, rr_intervals[i + 1], method=method, custom_rule=custom_rule):
                nn_intervals.append(rr_intervals[i + 1])
            else:
                # A débattre, Comment remplacer les outliers ?
                nn_intervals.append(np.nan)
                outlier_count += 1
                previous_outlier = True

    print("{} ectopic beat(s) have been deleted with {} rule.".format(outlier_count, method))

    return nn_intervals


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
    if method == MALIK_RULE:
        return abs(rr_interval - next_rr_interval) <= 0.2 * rr_interval
    elif method == KAMATH_RULE:
        return 0 <= (next_rr_interval - rr_interval) <= 0.325 * rr_interval or 0 <= \
               (rr_interval - next_rr_interval) <= 0.245 * rr_interval
    elif method == CUSTOM_RULE:
        return abs(rr_interval - next_rr_interval) <= custom_rule * rr_interval
    else:
        raise ValueError("Not a valid method. Please choose Malik, Kamath or custom.")


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
