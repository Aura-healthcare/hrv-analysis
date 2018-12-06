#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides several methods to clean abnormal and ectopic RR-intervals."""

from typing import Tuple
from typing import List
import pandas as pd
import numpy as np

# Static name for methods params
MALIK_RULE = "malik"
KARLSSON_RULE = "karlsson"
KAMATH_RULE = "kamath"
ACAR_RULE = "acar"
CUSTOM_RULE = "custom"

__all__ = ["remove_outliers", "remove_ectopic_beats", "interpolate_nan_values", "get_nn_intervals"]

# ----------------- ClEAN OUTlIER / ECTOPIC BEATS ----------------- #


def remove_outliers(rr_intervals: List[float], verbose: bool = True, low_rri: int = 300,
                    high_rri: int = 2000) -> List[float]:
    """
    Function that replace RR-interval outlier by nan.

    Parameters
    ---------
    rr_intervals : list
        raw signal extracted.
    low_rri : int
        lowest RrInterval to be considered plausible.
    high_rri : int
        highest RrInterval to be considered plausible.
    verbose : bool
        Print information about deleted outliers.

    Returns
    ---------
    rr_intervals_cleaned : list
        list of RR-intervals without outliers

    References
    ----------
    .. [1] O. Inbar, A. Oten, M. Scheinowitz, A. Rotstein, R. Dlin, R.Casaburi. Normal \
    cardiopulmonary responses during incremental exercise in 20-70-yr-old men.

    .. [2] W. C. Miller, J. P. Wallace, K. E. Eggert. Predicting max HR and the HR-VO2 relationship\
    for exercise prescription in obesity.

    .. [3] H. Tanaka, K. D. Monahan, D. R. Seals. Age-predictedmaximal heart rate revisited.

    .. [4] M. Gulati, L. J. Shaw, R. A. Thisted, H. R. Black, C. N. B.Merz, M. F. Arnsdorf. Heart \
    rate response to exercise stress testing in asymptomatic women.
    """

    # Conversion RrInterval to Heart rate ==> rri (ms) =  1000 / (bpm / 60)
    # rri 2000 => bpm 30 / rri 300 => bpm 200
    rr_intervals_cleaned = [rri if high_rri >= rri >= low_rri else np.nan for rri in rr_intervals]

    if verbose:
        outliers_list = []
        for rri in rr_intervals:
            if high_rri >= rri >= low_rri:
                pass
            else:
                outliers_list.append(rri)

        nan_count = sum(np.isnan(rr_intervals_cleaned))
        if nan_count == 0:
            print("{} outlier(s) have been deleted.".format(nan_count))
        else:
            print("{} outlier(s) have been deleted.".format(nan_count))
            print("The outlier(s) value(s) are : {}".format(outliers_list))
    return rr_intervals_cleaned


def remove_ectopic_beats(rr_intervals: List[float], method: str = "malik",
                         custom_removing_rule: float = 0.2) -> List[float]:
    """
    RR-intervals differing by more than the removing_rule from the one proceeding it are removed.

    Parameters
    ---------
    rr_intervals : list
        list of RR-intervals
    method : str
        method to use to clean outlier. malik, kamath, karlsson, acar or custom.
    custom_removing_rule : int
        Percentage criteria of difference with previous RR-interval at which we consider
        that it is abnormal. If method is set to Karlsson, it is the percentage of difference
        between the absolute mean of previous and next RR-interval at which  to consider the beat
        as abnormal.

    Returns
    ---------
    nn_intervals : list
        list of NN Interval
    outlier_count : int
        Count of outlier detected in RR-interval list

    References
    ----------
    .. [5] Kamath M.V., Fallen E.L.: Correction of the Heart Rate Variability Signal for Ectopics \
    and Miss- ing Beats, In: Malik M., Camm A.J.

    .. [6] Geometric Methods for Heart Rate Variability Assessment - Malik M et al
    """
    if method not in [MALIK_RULE, KAMATH_RULE, KARLSSON_RULE, ACAR_RULE, CUSTOM_RULE]:
        raise ValueError("Not a valid method. Please choose between malik, kamath, karlsson, acar.\
         You can also choose your own removing critera with custom_rule parameter.")

    if method == KARLSSON_RULE:
        nn_intervals, outlier_count = _remove_outlier_karlsson(rr_intervals=rr_intervals,
                                                               removing_rule=custom_removing_rule)

    elif method == ACAR_RULE:
        nn_intervals, outlier_count = _remove_outlier_acar(rr_intervals=rr_intervals)

    else:
        # set first element in list
        outlier_count = 0
        previous_outlier = False
        nn_intervals = [rr_intervals[0]]
        for i, rr_interval in enumerate(rr_intervals[:-1]):

            if previous_outlier:
                nn_intervals.append(rr_intervals[i + 1])
                previous_outlier = False
                continue

            if is_outlier(rr_interval, rr_intervals[i + 1], method=method,
                          custom_rule=custom_removing_rule):
                nn_intervals.append(rr_intervals[i + 1])
            else:
                nn_intervals.append(np.nan)
                outlier_count += 1
                previous_outlier = True

    print("{} ectopic beat(s) have been deleted with {} rule.".format(outlier_count, method))

    return nn_intervals


def is_outlier(rr_interval: int, next_rr_interval: float, method: str = "malik",
               custom_rule: float = 0.2) -> bool:
    """
    Test if the rr_interval is an outlier

    Parameters
    ----------
    rr_interval : int
        RrInterval
    next_rr_interval : int
        consecutive RrInterval
    method : str
        method to use to clean outlier. malik, kamath, karlsson, acar or custom
    custom_rule : int
        percentage criteria of difference with previous RR-interval at which we consider
        that it is abnormal

    Returns
    ----------
    outlier : bool
        True if RrInterval is valid, False if not
    """
    if method == MALIK_RULE:
        outlier = abs(rr_interval - next_rr_interval) <= 0.2 * rr_interval
    elif method == KAMATH_RULE:
        outlier = 0 <= (next_rr_interval - rr_interval) <= 0.325 * rr_interval or 0 <= \
                  (rr_interval - next_rr_interval) <= 0.245 * rr_interval
    else:
        outlier = abs(rr_interval - next_rr_interval) <= custom_rule * rr_interval
    return outlier


def _remove_outlier_karlsson(rr_intervals: List[float], removing_rule: float = 0.2) -> Tuple[List[float], int]:
    """
    RR-intervals differing by more than the 20 % of the mean of previous and next RR-interval
    are removed.

    Parameters
    ---------
    rr_intervals : list
        list of RR-intervals
    removing_rule : float
        Percentage of difference between the absolute mean of previous and next RR-interval at which \
    to consider the beat as abnormal.

    Returns
    ---------
    nn_intervals : list
        list of NN Interval

    References
    ----------
    .. [7]  Automatic filtering of outliers in RR-intervals before analysis of heart rate \
    variability in Holter recordings: a comparison with carefully edited data - Marcus Karlsson, \
    Rolf Hörnsten, Annika Rydberg and Urban Wiklund
    """
    # set first element in list
    nn_intervals = [rr_intervals[0]]
    outlier_count = 0

    for i in range(len(rr_intervals)):
        # Condition to go out of loop at limits of list
        if i == len(rr_intervals)-2:
            nn_intervals.append(rr_intervals[i + 1])
            break
        mean_prev_next_rri = (rr_intervals[i] + rr_intervals[i+2]) / 2
        if abs(mean_prev_next_rri - rr_intervals[i+1]) < removing_rule * mean_prev_next_rri:
            nn_intervals.append(rr_intervals[i+1])
        else:
            nn_intervals.append(np.nan)
            outlier_count += 1
    return nn_intervals, outlier_count


def _remove_outlier_acar(rr_intervals: List[float], custom_rule=0.2) -> Tuple[List[float], int]:
    """
    RR-intervals differing by more than the 20 % of the mean of last 9 RrIntervals
    are removed.

    Parameters
    ---------
    rr_intervals : list
        list of RR-intervals
    custom_rule : int
        percentage criteria of difference with mean of  9 previous RR-intervals at
        which we consider that RR-interval is abnormal. By default, set to 20 %

    Returns
    ---------
    nn_intervals : list
        list of NN Interval

    References
    ----------
    .. [8] Automatic ectopic beat elimination in short-term heart rate variability measurements \
    Acar B., Irina S., Hemingway H., Malik M.:
    """
    nn_intervals = []
    outlier_count = 0
    for i, rr_interval in enumerate(rr_intervals):
        if i < 9:
            nn_intervals.append(rr_interval)
            continue
        acar_rule_elt = np.nanmean(nn_intervals[-9:])
        if abs(acar_rule_elt - rr_interval) < custom_rule * acar_rule_elt:
            nn_intervals.append(rr_interval)
        else:
            nn_intervals.append(np.nan)
            outlier_count += 1
    return nn_intervals, outlier_count


def interpolate_nan_values(rr_intervals: List[float], interpolation_method: str = "linear", limit=1):
    """
    Function that interpolate Nan values with linear interpolation

    Parameters
    ---------
    rr_intervals : list
        RrIntervals list.
    interpolation_method : str
        Method used to interpolate Nan values of series.
    limit: int
        TODO
    Returns
    ---------
    interpolated_rr_intervals : list
        new list with outliers replaced by interpolated values.
    """
    series_rr_intervals_cleaned = pd.Series(rr_intervals)
    # Interpolate nan values and convert pandas object to list of values
    interpolated_rr_intervals = series_rr_intervals_cleaned.interpolate(method=interpolation_method,
                                                                        limit=limit,
                                                                        limit_area="inside")
    return interpolated_rr_intervals.values.tolist()


def get_nn_intervals(rr_intervals: List[float], low_rri: int = 300, high_rri: int = 2000,
                     interpolation_method: str = "linear", ectopic_beats_removal_method: str = KAMATH_RULE,
                     verbose: bool = True) -> List[float]:
    """
    Function that computes NN Intervals from RR-intervals.

    Parameters
    ---------
    rr_intervals : list
        RrIntervals list.
    interpolation_method : str
        Method used to interpolate Nan values of series.
    ectopic_beats_removal_method : str
        method to use to clean outlier. malik, kamath, karlsson, acar or custom.
    low_rri : int
        lowest RrInterval to be considered plausible.
    high_rri : int
        highest RrInterval to be considered plausible.
    verbose : bool
        Print information about deleted outliers.

    Returns
    ---------
    interpolated_nn_intervals : list
        list of NN Interval interpolated
    """
    rr_intervals_cleaned = remove_outliers(rr_intervals, low_rri=low_rri, high_rri=high_rri,
                                           verbose=verbose)
    interpolated_rr_intervals = interpolate_nan_values(rr_intervals_cleaned, interpolation_method)
    nn_intervals = remove_ectopic_beats(interpolated_rr_intervals,
                                        method=ectopic_beats_removal_method)
    interpolated_nn_intervals = interpolate_nan_values(nn_intervals, interpolation_method)
    return interpolated_nn_intervals


def is_valid_sample(nn_intervals: List[float], outlier_count: int, removing_rule: float = 0.04) -> bool:
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
        print("Too much outlier for analyses ! You should descard the sample.")
        result = False
    if len(nn_intervals) < 240:
        print("Not enough Heart beat for Nyquist criteria ! ")
        result = False
    return result
