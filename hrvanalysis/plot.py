#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides several methods to plot RR / NN Intervals."""

from typing import List
import matplotlib.pyplot as plt
from matplotlib import style
from hrvanalysis.extract_features import _get_freq_psd_from_nn_intervals
from collections import namedtuple
import numpy as np

# Named Tuple for different frequency bands
VlfBand = namedtuple("Vlf_band", ["low", "high"])
LfBand = namedtuple("Lf_band", ["low", "high"])
HfBand = namedtuple("Hf_band", ["low", "high"])


def plot_timeseries(nn_intervals: List[int], normalize: bool = True):
    """
    Function plotting the NN-intervals timeseries.

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    normalize : bool
        Set to True to plot X axis as a cumulative sum of Time.
        Set to False to plot X axis using x as index array 0, 1, ..., N-1.
    """

    style.use("seaborn-darkgrid")
    plt.figure(figsize=(12, 8))
    plt.title("Rr Interval time series")
    plt.ylabel("Rr Interval", fontsize=15)

    if normalize:
        plt.xlabel("Time (s)", fontsize=15)
        plt.plot(np.cumsum(nn_intervals) / 1000, nn_intervals)
    else:
        plt.xlabel("RR-interval index", fontsize=15)
        plt.plot(nn_intervals)
    plt.show()


def plot_distrib(nn_intervals: List[int], bin_length: int = 8):
    """
    Function plotting histogram distribution of the NN Intervals. Useful for geometrical features.

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    bin_length : int
        size of the bin for histogram in ms, by default = 8.
    """

    max_nn_i = max(nn_intervals)
    min_nn_i = min(nn_intervals)

    style.use("seaborn-darkgrid")
    plt.figure(figsize=(12, 8))
    plt.title("Distribution of Rr Intervals", fontsize=20)
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Number of Rr Interval per bin", fontsize=15)
    plt.hist(nn_intervals, bins=range(min_nn_i - 10, max_nn_i + 10, bin_length))
    plt.show()


def plot_psd(nn_intervals: List[int], method: str = "welch", sampling_frequency: int = 7,
             interpolation_method: str = "linear", vlf_band: namedtuple = VlfBand(0.0033, 0.04),
             lf_band: namedtuple = LfBand(0.04, 0.15), hf_band: namedtuple = HfBand(0.15, 0.40)):
    """
    Function plotting the power spectral density of the NN Intervals.

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    method : str
        Method used to calculate the psd. Choice are Welch's FFT (welch) or Lomb method (lomb).
    sampling_frequency : int
        frequence at which the signal is sampled. Common value range from 1 Hz to 10 Hz, by default
        set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        kind of interpolation as a string, by default "linear". No need to specify if lomb method is
        used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.
    """

    freq, psd = _get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method=method,
                                                sampling_frequency=sampling_frequency,
                                                interpolation_method=interpolation_method)

    # Calcul of indices between desired frequency bands
    vlf_indexes = np.logical_and(freq >= vlf_band[0], freq < vlf_band[1])
    lf_indexes = np.logical_and(freq >= lf_band[0], freq < lf_band[1])
    hf_indexes = np.logical_and(freq >= hf_band[0], freq < hf_band[1])

    frequency_band_index = [vlf_indexes, lf_indexes, hf_indexes]
    label_list = ["VLF component", "LF component", "HF component"]

    # Plot parameters
    style.use("seaborn-darkgrid")
    plt.figure(figsize=(12, 8))
    plt.xlabel("Frequency (Hz)", fontsize=15)
    plt.ylabel("PSD (s2/ Hz)", fontsize=15)

    if method == "lomb":
        plt.title("Lomb's periodogram", fontsize=20)
        for band_index, label in zip(frequency_band_index, label_list):
            plt.fill_between(freq[band_index], 0, psd[band_index] / (1000 * len(psd[band_index])), label=label)
        plt.legend(prop={"size": 15}, loc="best")

    elif method == "welch":
        plt.title("FFT Spectrum : Welch's periodogram", fontsize=20)
        for band_index, label in zip(frequency_band_index, label_list):
            plt.fill_between(freq[band_index], 0, psd[band_index] / (1000 * len(psd[band_index])), label=label)
        plt.legend(prop={"size": 15}, loc="best")
        plt.xlim(0, 0.5)
    else:
        raise ValueError("Not a valid method. Choose between 'lomb' and 'welch'")

    plt.show()


def plot_poincare(nn_intervals: List[int]):
    """
    Poincaré / Lorentz Plot of the NN Intervals.

    Arguments
    ---------
    nn_intervals : list
        list of NN intervals.

    Notes
    ---------
    The transverse axis (T) reflects beat-to-beat variation.
    the longitudinal axis (L) reflects the overall fluctuation.
    """

    ax1 = nn_intervals[:-1]
    ax2 = nn_intervals[1:]

    style.use("seaborn-darkgrid")
    plt.figure(figsize=(12, 8))
    plt.title("Poincaré / Lorentz Plot", fontsize=20)
    plt.xlabel('NN_n (s)', fontsize=15)
    plt.ylabel('NN_n+1 (s)', fontsize=15)
    plt.scatter(ax1, ax2, c='r', s=12)
    plt.show()
