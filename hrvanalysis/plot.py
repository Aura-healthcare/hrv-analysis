#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides several methods to plot RR / NN Intervals."""

from typing import List
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib.patches import Ellipse
from hrvanalysis.extract_features import _get_freq_psd_from_nn_intervals
from hrvanalysis.extract_features import get_poincare_plot_features
from collections import namedtuple
import numpy as np

# Named Tuple for different frequency bands
VlfBand = namedtuple("Vlf_band", ["low", "high"])
LfBand = namedtuple("Lf_band", ["low", "high"])
HfBand = namedtuple("Hf_band", ["low", "high"])

# TODO : ajouter argument pour passer paramètre xmin / max & ymin / ymax plot OU Auto scale


def plot_timeseries(nn_intervals: List[float], normalize: bool = True,
                    autoscale: bool = True, y_min: float = None, y_max: float = None):
    """
    Function plotting the NN-intervals timeseries.
    TODO

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    normalize : bool
        Set to True to plot X axis as a cumulative sum of Time.
        Set to False to plot X axis using x as index array 0, 1, ..., N-1.
    autoscale : bool
    y_min : float
    y_max : float
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

    if not autoscale:
        plt.ylim(y_min, y_max)
    plt.show()


def plot_distrib(nn_intervals: List[float], bin_length: int = 8):
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


def plot_psd(nn_intervals: List[float], method: str = "welch", sampling_frequency: int = 7,
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
        plt.xlim(0, hf_band[1])
    else:
        raise ValueError("Not a valid method. Choose between 'lomb' and 'welch'")

    plt.show()


def plot_poincare(nn_intervals: List[float], plot_sd_features=True):
    """
    Pointcare / Lorentz Plot of the NN Intervals

    Arguments
    ---------
    nn_intervals : list
        list of NN intervals
    plot_sd_features : bool

    Notes
    ---------
    The transverse axis (T) reflects beat-to-beat variation
    the longitudinal axis (L) reflects the overall fluctuation
    """
    # For Lorentz / poincaré Plot
    ax1 = nn_intervals[:-1]
    ax2 = nn_intervals[1:]

    # compute features for ellipse's height, width and center
    dict_sd1_sd2 = get_poincare_plot_features(nn_intervals)
    sd1 = dict_sd1_sd2["sd1"]
    sd2 = dict_sd1_sd2["sd2"]
    mean_nni = np.mean(nn_intervals)

    # Plot options and settings
    style.use("seaborn-darkgrid")
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    plt.title("Poincaré / Lorentz Plot", fontsize=20)
    plt.xlabel('NN_n (s)', fontsize=15)
    plt.ylabel('NN_n+1 (s)', fontsize=15)
    plt.xlim(min(nn_intervals) - 10, max(nn_intervals) + 10)
    plt.ylim(min(nn_intervals) - 10, max(nn_intervals) + 10)

    # Poincaré Plot
    ax.scatter(ax1, ax2, c='b', s=2)

    if plot_sd_features:
        # Ellipse plot settings
        ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2 + 1,
                       height=2 * sd1 + 1, angle=45, linewidth=2,
                       fill=False)
        ax.add_patch(ells)

        ells = Ellipse(xy=(mean_nni, mean_nni), width=2 * sd2,
                       height=2 * sd1, angle=45)
        ells.set_alpha(0.05)
        ells.set_facecolor("blue")
        ax.add_patch(ells)

        # Arrow plot settings
        sd1_arrow = ax.arrow(mean_nni, mean_nni, -sd1 * np.sqrt(2) / 2, sd1 * np.sqrt(2) / 2,
                             linewidth=3, ec='r', fc="r", label="SD1")
        sd2_arrow = ax.arrow(mean_nni, mean_nni, sd2 * np.sqrt(2) / 2, sd2 * np.sqrt(2) / 2,
                             linewidth=3, ec='g', fc="g", label="SD2")

        plt.legend(handles=[sd1_arrow, sd2_arrow], fontsize=12, loc="best")
    plt.show()
