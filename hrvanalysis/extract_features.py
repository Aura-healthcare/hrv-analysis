#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This script provides several methods to extract features from Normal to Normal Intervals
 for heart rate variability analysis."""

from typing import List, Tuple
from collections import namedtuple
import numpy as np
import nolds
from scipy import interpolate
from scipy import signal
from astropy.stats import LombScargle

# limit functions that user might import using "from hrv-analysis import *"
__all__ = ['get_time_domain_features', 'get_frequency_domain_features',
           'get_geometrical_features', 'get_poincare_plot_features',
           "get_csi_cvi_features", "get_sampen"]

# Frequency Methods name
WELCH_METHOD = "welch"
LOMB_METHOD = "lomb"

# Named Tuple for different frequency bands
VlfBand = namedtuple("Vlf_band", ["low", "high"])
LfBand = namedtuple("Lf_band", ["low", "high"])
HfBand = namedtuple("Hf_band", ["low", "high"])

# ----------------- TIME DOMAIN FEATURES ----------------- #


def get_time_domain_features(nn_intervals: List[float]) -> dict:
    """
    Returns a dictionary containing time domain features for HRV analysis.
    Mostly used on long term recordings (24h) but some studies use some of those features on
    short term recordings, from 2 to 5 minutes window.

    Parameters
    ----------
    nn_intervals : list
        list of Normal to Normal Interval

    Returns
    -------
    time_domain_features : dict
        dictionary containing time domain features for HRV analyses. There are details
        about each features below.

    Notes
    -----
    Here are some details about feature engineering...

    - **mean_nni**: The mean of RR-intervals.

    - **sdnn** : The standard deviation of the time interval between successive normal heart beats \
    (i.e. the RR-intervals).

    - **sdsd**: The standard deviation of differences between adjacent RR-intervals

    - **rmssd**: The square root of the mean of the sum of the squares of differences between \
    adjacent NN-intervals. Reflects high frequency (fast or parasympathetic) influences on hrV \
    (*i.e.*, those influencing larger changes from one beat to the next).

    - **median_nni**: Median Absolute values of the successive differences between the RR-intervals.

    - **nni_50**: Number of interval differences of successive RR-intervals greater than 50 ms.

    - **pnni_50**: The proportion derived by dividing nni_50 (The number of interval differences \
    of successive RR-intervals greater than 50 ms) by the total number of RR-intervals.

    - **nni_20**: Number of interval differences of successive RR-intervals greater than 20 ms.

    - **pnni_20**: The proportion derived by dividing nni_20 (The number of interval differences \
    of successive RR-intervals greater than 20 ms) by the total number of RR-intervals.

    - **range_nni**: difference between the maximum and minimum nn_interval.

    - **cvsd**: Coefficient of variation of successive differences equal to the rmssd divided by \
    mean_nni.

    - **cvnni**: Coefficient of variation equal to the ratio of sdnn divided by mean_nni.

    - **mean_hr**: The mean Heart Rate.

    - **max_hr**: Max heart rate.

    - **min_hr**: Min heart rate.

    - **std_hr**: Standard deviation of heart rate.

    References
    ----------
    .. [1] Heart rate variability - Standards of measurement, physiological interpretation, and \
    clinical use, Task Force of The European Society of Cardiology and The North American Society \
    of Pacing and Electrophysiology, 1996
    """

    diff_nni = np.diff(nn_intervals)
    length_int = len(nn_intervals)

    # Basic statistics
    mean_nni = np.mean(nn_intervals)
    median_nni = np.median(nn_intervals)
    range_nni = max(nn_intervals) - min(nn_intervals)

    sdsd = np.std(diff_nni)
    rmssd = np.sqrt(np.mean(diff_nni ** 2))

    nni_50 = sum(np.abs(diff_nni) > 50)
    pnni_50 = 100 * nni_50 / length_int
    nni_20 = sum(np.abs(diff_nni) > 20)
    pnni_20 = 100 * nni_20 / length_int

    # Feature found on github and not in documentation
    cvsd = rmssd / mean_nni

    # Features only for long term recordings
    sdnn = np.std(nn_intervals, ddof=1)  # ddof = 1 : unbiased estimator => divide std by n-1
    cvnni = sdnn / mean_nni

    # Heart Rate equivalent features
    heart_rate_list = np.divide(60000, nn_intervals)
    mean_hr = np.mean(heart_rate_list)
    min_hr = min(heart_rate_list)
    max_hr = max(heart_rate_list)
    std_hr = np.std(heart_rate_list)

    time_domain_features = {
        'mean_nni': mean_nni,
        'sdnn': sdnn,
        'sdsd': sdsd,
        'nni_50': nni_50,
        'pnni_50': pnni_50,
        'nni_20': nni_20,
        'pnni_20': pnni_20,
        'rmssd': rmssd,
        'median_nni': median_nni,
        'range_nni': range_nni,
        'cvsd': cvsd,
        'cvnni': cvnni,
        'mean_hr': mean_hr,
        "max_hr": max_hr,
        "min_hr": min_hr,
        "std_hr": std_hr,
    }

    return time_domain_features


def get_geometrical_features(nn_intervals: List[float]) -> dict:
    """
    Returns a dictionary containing geometrical time domain features for HRV analyses.
    Must use this function on recordings from 20 minutes to 24 Hours window.

    Parameters
    ---------
    nn_intervals : list
        list of Normal to Normal Interval.

    Returns
    ---------
    geometrical_features : dict
        Dictionary containing geometrical time domain features for HRV analyses.
        There are details about each features below.

    Notes
    ----------
    Details about feature engineering...

    - **triangular_index**: The HRV triangular index measurement is the integral of the density \
    distribution (= the number of all NN-intervals) divided by the maximum of the density \
    distribution.

    - **tinn**: The triangular interpolation of NN-interval histogram (TINN) is the baseline width \
     of the distribution measured as a base of a triangle, approximating the NN-interval \
     distribution

    References
    ----------
    .. [2] Heart rate variability - Standards of measurement, physiological interpretation, and \
    clinical use, Task Force of The European Society of Cardiology and The North American Society \
    of Pacing and Electrophysiology, 1996

    """

    triang_idx = len(nn_intervals) / max(np.histogram(nn_intervals, bins=range(300, 2000, 8))[0])
    # TODO
    tinn = None

    geometrical_features = {
        "triangular_index": triang_idx,
        "tinn": tinn
    }

    return geometrical_features


# ----------------- FREQUENCY DOMAIN FEATURES ----------------- #


def get_frequency_domain_features(nn_intervals: List[float], method: str = WELCH_METHOD,
                                  sampling_frequency: int = 7, interpolation_method: str = "linear",
                                  vlf_band: namedtuple = VlfBand(0.0033, 0.04),
                                  lf_band: namedtuple = LfBand(0.04, 0.15),
                                  hf_band: namedtuple = HfBand(0.15, 0.40)) -> dict:
    """
    Returns a dictionary containing frequency domain features for HRV analyses.
    Must use this function on short term recordings, from 2 to 5 minutes window.

    Parameters
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    method : str
        Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency : int
        Frequency at which the signal is sampled. Common value range from 1 Hz to 10 Hz,
        by default set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        kind of interpolation as a string, by default "linear". No need to specify if Lomb
        method is used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.

    Returns
    ---------
    frequency_domain_features : dict
        Dictionary containing frequency domain features for HRV analyses. There are details
        about each features below.

    Notes
    ---------
    Details about feature engineering...

    - **total_power** : Total power density spectral

    - **vlf** : variance ( = power ) in HRV in the Very low Frequency (.003 to .04 Hz by default). \
    Reflect an intrinsic rhythm produced by the heart which is modulated primarily by sympathetic \
    activity.

    - **lf** : variance ( = power ) in HRV in the low Frequency (.04 to .15 Hz). Reflects a \
    mixture of sympathetic and parasympathetic activity, but in long-term recordings, it reflects \
    sympathetic activity and can be reduced by the beta-adrenergic antagonist propanolol.

    - **hf**: variance ( = power ) in HRV in the High Frequency (.15 to .40 Hz by default). \
    Reflects fast changes in beat-to-beat variability due to parasympathetic (vagal) activity. \
    Sometimes called the respiratory band because it corresponds to HRV changes related to the \
    respiratory cycle and can be increased by slow, deep breathing (about 6 or 7 breaths per \
    minute) and decreased by anticholinergic drugs or vagal blockade.

    - **lf_hf_ratio** : lf/hf ratio is sometimes used by some investigators as a quantitative \
    mirror of the sympatho/vagal balance.

    - **lfnu** : normalized lf power

    - **hfnu** : normalized hf power

    References
    ----------
    .. [3] Heart rate variability - Standards of measurement, physiological interpretation, and \
    clinical use, Task Force of The European Society of Cardiology and The North American Society \
    of Pacing and Electrophysiology, 1996

    .. [4] Signal Processing Methods for Heart Rate Variability - Gari D. Clifford, 2002

    """

    # ----------  Compute frequency & Power of signal  ---------- #
    freq, psd = _get_freq_psd_from_nn_intervals(nn_intervals=nn_intervals, method=method,
                                                sampling_frequency=sampling_frequency,
                                                interpolation_method=interpolation_method,
                                                vlf_band=vlf_band, hf_band=hf_band)

    # ---------- Features calculation ---------- #
    freqency_domain_features = _get_features_from_psd(freq=freq, psd=psd,
                                                      vlf_band=vlf_band,
                                                      lf_band=lf_band,
                                                      hf_band=hf_band)

    return freqency_domain_features


def _get_freq_psd_from_nn_intervals(nn_intervals: List[float], method: str = WELCH_METHOD,
                                    sampling_frequency: int = 7,
                                    interpolation_method: str = "linear",
                                    vlf_band: namedtuple = VlfBand(0.0033, 0.04),
                                    hf_band: namedtuple = HfBand(0.15, 0.40)) -> Tuple:
    """
    Returns the frequency and power of the signal.

    Parameters
    ---------
    nn_intervals : list
        list of Normal to Normal Interval
    method : str
        Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency : int
        Frequency at which the signal is sampled. Common value range from 1 Hz to 10 Hz,
        by default set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method : str
        Kind of interpolation as a string, by default "linear". No need to specify if Lomb
        method is used.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.

    Returns
    ---------
    freq : list
        Frequency of the corresponding psd points.
    psd : list
        Power Spectral Density of the signal.

    """

    timestamps = _create_time_info(nn_intervals)

    if method == WELCH_METHOD:
        # ---------- Interpolation of signal ---------- #
        funct = interpolate.interp1d(x=timestamps, y=nn_intervals, kind=interpolation_method)

        timestamps_interpolation = _create_interpolation_time(nn_intervals, sampling_frequency)
        nni_interpolation = funct(timestamps_interpolation)

        # ---------- Remove DC Component ---------- #
        nni_normalized = nni_interpolation - np.mean(nni_interpolation)

        #  ----------  Compute Power Spectral Density  ---------- #
        freq, psd = signal.welch(x=nni_normalized, fs=sampling_frequency, window='hann',
                                 nfft=4096)

    elif method == LOMB_METHOD:
        freq, psd = LombScargle(timestamps, nn_intervals,
                                normalization='psd').autopower(minimum_frequency=vlf_band[0],
                                                               maximum_frequency=hf_band[1])
    else:
        raise ValueError("Not a valid method. Choose between 'lomb' and 'welch'")

    return freq, psd


def _create_time_info(nn_intervals: List[float]) -> List[float]:
    """
    Creates corresponding time interval for all nn_intervals

    Parameters
    ---------
    nn_intervals : list
        List of Normal to Normal Interval.

    Returns
    ---------
    nni_tmstp : list
        list of time intervals between first NN-interval and final NN-interval.
    """
    # Convert in seconds
    nni_tmstp = np.cumsum(nn_intervals) / 1000

    # Force to start at 0
    return nni_tmstp - nni_tmstp[0]


def _create_interpolation_time(nn_intervals: List[float], sampling_frequency: int = 7) -> List[float]:
    """
    Creates the interpolation time used for Fourier transform's method

    Parameters
    ---------
    nn_intervals : list
        List of Normal to Normal Interval.
    sampling_frequency : int
        Frequency at which the signal is sampled.

    Returns
    ---------
    nni_interpolation_tmstp : list
        Timestamp for interpolation.
    """
    time_nni = _create_time_info(nn_intervals)
    # Create timestamp for interpolation
    nni_interpolation_tmstp = np.arange(0, time_nni[-1], 1 / float(sampling_frequency))
    return nni_interpolation_tmstp


def _get_features_from_psd(freq: List[float], psd: List[float], vlf_band: namedtuple = VlfBand(0.0033, 0.04),
                           lf_band: namedtuple = LfBand(0.04, 0.15),
                           hf_band: namedtuple = HfBand(0.15, 0.40)) -> dict:
    """
    Computes frequency domain features from the power spectral decomposition.

    Parameters
    ---------
    freq : array
        Array of sample frequencies.
    psd : list
        Power spectral density or power spectrum.
    vlf_band : tuple
        Very low frequency bands for features extraction from power spectral density.
    lf_band : tuple
        Low frequency bands for features extraction from power spectral density.
    hf_band : tuple
        High frequency bands for features extraction from power spectral density.

    Returns
    ---------
    freqency_domain_features : dict
        Dictionary containing frequency domain features for HRV analyses. There are details
        about each features given below.
    """

    # Calcul of indices between desired frequency bands
    vlf_indexes = np.logical_and(freq >= vlf_band[0], freq < vlf_band[1])
    lf_indexes = np.logical_and(freq >= lf_band[0], freq < lf_band[1])
    hf_indexes = np.logical_and(freq >= hf_band[0], freq < hf_band[1])

    # STANDARDS

    # Integrate using the composite trapezoidal rule
    lf = np.trapz(y=psd[lf_indexes], x=freq[lf_indexes])
    hf = np.trapz(y=psd[hf_indexes], x=freq[hf_indexes])

    # total power & vlf : Feature often used for  "long term recordings" analysis
    vlf = np.trapz(y=psd[vlf_indexes], x=freq[vlf_indexes])
    total_power = vlf + lf + hf

    lf_hf_ratio = lf / hf
    lfnu = (lf / (lf + hf)) * 100
    hfnu = (hf / (lf + hf)) * 100

    freqency_domain_features = {
        'lf': lf,
        'hf': hf,
        'lf_hf_ratio': lf_hf_ratio,
        'lfnu': lfnu,
        'hfnu': hfnu,
        'total_power': total_power,
        'vlf': vlf
    }

    return freqency_domain_features


# ----------------- NON lINEAR DOMAIN FEATURES ----------------- #


def get_csi_cvi_features(nn_intervals: List[float]) -> dict:
    """
    Returns a dictionary containing 3 features from non linear domain for hrV analyses.
    Must use this function on short term recordings, for 30 , 50, 100 RR-intervals (or
    seconds) window.

    Parameters
    ---------
    nn_intervals : list
        Normal to Normal Intervals.

    Returns
    ---------
    csi_cvi_features : dict
        Dictionary containing non linear domain features for hrV analyses. There are  details about
        each features are given below.

    Notes
    ---------
    - **csi** : Cardiac Sympathetic Index.

    - **cvi** : Cadiac Vagal Index.

    - **Modified_csi** : Modified CSI is an alternative measure in research of seizure detection.

    References
    ----------
    .. [5] Using Lorenz plot and Cardiac Sympathetic Index of heart rate variability for detecting \
    seizures for patients with epilepsy, Jesper Jeppesen et al, 2014

    """

    # Measures the width and length of poincare cloud
    poincare_plot_features = get_poincare_plot_features(nn_intervals)
    T = 4 * poincare_plot_features['sd1']
    L = 4 * poincare_plot_features['sd2']

    csi = L / T
    cvi = np.log10(L * T)
    modified_csi = L ** 2 / T

    csi_cvi_features = {
        'csi': csi,
        'cvi': cvi,
        'Modified_csi': modified_csi
    }

    return csi_cvi_features


def get_poincare_plot_features(nn_intervals: List[float]) -> dict:
    """
    Function returning a dictionary containing 3 features from non linear domain
    for hrV analyses.
    Must use this function on short term recordings, from 5 minutes window.

    Parameters
    ---------
    nn_intervals : list
        Normal to Normal Interval

    Returns
    ---------
    poincare_plot_features : dict
        dictionary containing non linear domain features for hrV analyses. There
        are details about each features are given below.

    Notes
    ---------
    - **sd1** : The standard deviation of projection of the Poincaré plot on the line \
    perpendicular to the line of identity.

    - **sd2** : SD2 is defined as the standard deviation of the projection of the Poincaré \
    plot on the line of identity (y=x).

    - **ratio_sd2_sd1** : Ratio between SD2 and SD1.

    References
    ----------
    .. [6] Pre-ictal heart rate variability assessment of epileptic seizures by means of linear \
    and non- linear analyses, Soroor Behbahani, Nader Jafarnia Dabanloo et al - 2013

    """
    diff_nn_intervals = np.diff(nn_intervals)
    # measures the width of poincare cloud
    sd1 = np.sqrt(np.std(diff_nn_intervals, ddof=1) ** 2 * 0.5)
    # measures the length of the poincare cloud
    sd2 = np.sqrt(2 * np.std(nn_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_nn_intervals, ddof=1) ** 2)
    ratio_sd2_sd1 = sd2 / sd1

    poincare_plot_features = {
        'sd1': sd1,
        'sd2': sd2,
        'ratio_sd2_sd1': ratio_sd2_sd1
    }

    return poincare_plot_features


def get_sampen(nn_intervals: List[float]) -> dict:
    """
    Function computing the sample entropy of the given data.
    Must use this function on short term recordings, from 1 minute window.

    Parameters
    ---------
    nn_intervals : list
        Normal to Normal Interval

    Returns
    ---------
    sampen : float
        The sample entropy of the data

    References
    ----------
    .. [7] Physiological time-series analysis using approximate entropy and sample entropy, \
    JOSHUA S. RICHMAN1, J. RANDALL MOORMAN - 2000

    """

    sampen = nolds.sampen(nn_intervals, emb_dim=2)
    return {'sampen': sampen}
