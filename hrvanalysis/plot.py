#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script provides several methods to plot RR / NN Intervals."""

import matplotlib.pyplot as plt
from matplotlib import style
from hrvanalysis.extract_features import get_freq_psd_from_nn_intervals


def plot_timeseries(nn_intervals):
    """
    Function plotting the NN Intervals timeseries

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    """

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Rr Interval time series")
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Rr Interval", fontsize=15)
    plt.plot(nn_intervals)
    plt.show()


def plot_distrib(nn_intervals, bin_length=8):
    """
    Function plotting histogram distribution of the NN Intervals. Useful for geometrical features

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    bin_length : int
        size of the bin for histogram in ms, by default = 8
    """

    max_nn_i = max(nn_intervals)
    min_nn_i = min(nn_intervals)

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Distribution of Rr Intervals", fontsize=20)
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Number of Rr Interval per bin", fontsize=15)
    plt.hist(nn_intervals, bins=range(min_nn_i - 10, max_nn_i + 10, bin_length))
    plt.show()


def plot_psd(nn_intervals, method="Welch", sampling_frequency=7, interpolation_method="linear"):
    """
    Function plotting the power spectral density of the NN Intervals

    Arguments
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.
    method : str
        Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency : int
        frequence at which the signal is sampled. Common value range from 1 Hz to 10 Hz, by default
        set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        kind of interpolation as a string, by default "linear". No need to specify if lomb method is
        used.
    """

    freq, psd = get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method=method,
                                               sampling_frequency=sampling_frequency,
                                               interpolation_method=interpolation_method)
    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.xlabel("Frequency (Hz)", fontsize=15)
    plt.ylabel("PSD (s2/ Hz)", fontsize=15)

    if method == "Lomb":
        plt.title("Lomb's periodogram", fontsize=20)
        plt.plot(freq, psd / (1000 * len(psd)))
    elif method == "Welch":
        plt.title("FFT Spectrum : Welch's periodogram", fontsize=20)
        plt.plot(freq, psd / (1000 * len(psd)))
        plt.xlim(0, 0.5)
    else:
        raise ValueError("Not a valid method. Choose between 'Lomb' and 'Welch'")

    plt.show()


def plot_poincare(nn_intervals):
    """
    Pointcare / Lorentz Plot of the NN Intervals

    Arguments
    ---------
    nn_intervals : list
        list of NN intervals

    Notes
    ---------
    The transverse axis (T) reflects beat-to-beat variation
    the longitudinal axis (L) reflects the overall fluctuation
    """

    ax1 = nn_intervals[:-1]
    ax2 = nn_intervals[1:]

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Pointcare / Lorentz Plot", fontsize=20)
    plt.xlabel('NN_n (s)', fontsize=15)
    plt.ylabel('NN_n+1 (s)', fontsize=15)
    plt.scatter(ax1, ax2, c='r', s=12)
    plt.show()
