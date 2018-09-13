import numpy as np
import nolds
from scipy import interpolate
from scipy import signal
from astropy.stats import LombScargle
import pandas as pd

# ----------------- ClEAN OUTlIER / ECTOPIC BEATS ----------------- #

# TO DO / ONGOING ...


def clean_outlier(rr_intervals, low_rri=300, high_rri=2000):
    """
    Function that replace RR Interval outlier by nan

    Arguments
    ---------
    rr_intervals - raw signal extracted
    low_rri - lowest RrInterval to be considered plausible
    high_rri - highest RrInterval to be considered plausible

    Returns
    ---------
    rr_intervals_cleaned - list of RR Intervals without outliers

    """

    # Conversion RrInterval to Heart rate ==> rri (ms) =  1000 / (bpm / 60)
    # rri 2000 => bpm 30 / rri 300 => bpm 200
    rr_intervals_cleaned = [x if high_rri >= x >= low_rri else np.nan for x in rr_intervals]
    nan_count = sum(np.isnan(rr_intervals_cleaned))
    print("{} outlier(s) have been deleted.".format(nan_count))
    return rr_intervals_cleaned


def interpolate_cleaned_outlier(rr_intervals_cleaned):
    """
    Function that interpolate Nan values with linear interpolation

    Arguments
    ---------
    rr_intervals_cleaned - RrIntervals without outliers values

    Returns
    ---------
    rr_intervals_interpolated - new list with outliers replaced by interpolated values
    """
    s = pd.Series(rr_intervals_cleaned)
    rr_intervals_interpolated = s.interpolate(method="linear")
    return rr_intervals_interpolated


def clean_ectopic_beats(rr_intervals, method="Malik", custom_rule=None):
    """
    RR intervals differing by more than the removing_rule from the one proceeding it are removed.

    Arguments
    ---------
    rr_intervals - list of Rr Intervals
    method - method to use to clean outlier. Malik, Kamath, Karlsson, mean_last9 or Custom
    custom_rule - percentage criteria of difference with previous Rr
    Interval at which we consider that it is abnormal

    Returns
    ---------
    nn_intervals - list of NN Interval

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

    if method == "mean_last9":
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
    if method == "Malik":
        return abs(rr_interval - next_rr_interval) <= 0.2 * rr_interval
    elif method == "Kamath":
        return 0 <= (next_rr_interval - rr_interval) <= 0.325 * rr_interval or 0 <= (rr_interval - next_rr_interval) \
               <= 0.245 * rr_interval
    elif method == "custom":
        return abs(rr_interval - next_rr_interval) <= custom_rule * rr_interval
    else:
        raise ValueError("Not a valid method. Please choose Malik or Kamath")


def is_valid_sample(nn_intervals, outlier_count, removing_rule=0.04):
    """
    Test if the sample meet the condition to be used for analysis

    Arguments
    ----------
    nn_intervals - list of Normal to Normal Interval
    outlier_count - count of outliers or ectopic beats removed from the interval
    removing_rule - rule to follow to determine whether the sample is valid or not

    Returns
    ----------
    Boolean - True if sample is valid, False if not
    """
    if outlier_count / len(nn_intervals) > removing_rule:
        print("Too much outlier for analyses ! You should descard the sample")
        return False
    if len(nn_intervals) < 240:
        print("Not enough Heart beat for Nyquist criteria ! ")
        return False
    return True


# ----------------- TIME DOMAIN FEATURES ----------------- #


def get_time_domain_features(nn_intervals):
    """
    Function returning a dictionnary containing time domain features for 
    HRV analysis.

    Mostly used on long term recordings (24h) but some studies use this 
    features on short term recordings, from 2 to 5 minutes window.
    
    Arguments
    ----------
    nn_intervals - list of Normal to Normal Interval
    
    Returns
    ----------
    time_domain_features - dictionnary containing time domain features for 
    HRV analyses. There are details about each features below.
    
    Notes
    ----------
    Details about feature engineering...
    
    - **mean_nni**: The mean of RR intervals
    
    - **sdnn**: The standard deviation of the time interval between successive normal heart 
    beats (i.e. the RR intervals)
    
    - **sdsd**: The standard deviation of differences between adjacent RR intervals
    
    - **rmssd**: The square root of the mean of the sum of the squares of differences between
    adjacent NN intervals. Reflects high frequency (fast or parasympathetic) influences on hrV 
    (*i.e.*, those influencing larger changes from one beat to the next).
    
    - **median_nni**: Median Absolute values of the successive differences between the RR intervals.

    - **nni_50**: Number of interval differences of successive RR intervals greater than 50 ms.

    - **pnni_50**: The proportion derived by dividing nni_50 (The number of interval differences of 
    successive RR intervals greater than 50 ms) by the total number of RR intervals.

    - **nni_20**: Number of interval differences of successive RR intervals greater than 20 ms.

    - **pnni_20**: The proportion derived by dividing nni_20 (The number of interval differences of 
    successive RR intervals greater than 20 ms) by the total number of RR intervals.

    - **cvsd**: The coefficient of variation of successive differences (van Dellen et al., 1985), 
    the rmssd divided by mean_nni.
    
    - **range_nni**: différence between the maximum and minimum nn_interval.
    
    
    References
    ----------
    - Signal Processing Methods for Heart Rate Variability - Gari D. Clifford, 2002
    - Heart rate variability - Standards of measurement, physiological interpretation, and clinical use, Task Force of
    The European Society of Cardiology and The North American Society of Pacing and Electrophysiology, 1996
    
    """
    
    diff_nni = np.diff(nn_intervals)
    length_int = len(nn_intervals)
    
    mean_nni = np.mean(nn_intervals)
    sdsd = np.std(diff_nni)
    rmssd = np.sqrt(np.mean(diff_nni ** 2))
    median_nni = np.median(nn_intervals)
    
    nni_50 = sum(abs(diff_nni) > 50)
    pnni_50 = 100 * nni_50 / length_int
    nni_20 = sum(abs(diff_nni) > 20)
    pnni_20 = 100 * nni_20 / length_int
    
    range_nni = max(nn_intervals) - min(nn_intervals)
    
    # Feature found on github et non in documentation
    cvsd = rmssd / mean_nni
    
    # Features only for long term recordings
    sdnn = np.std(nn_intervals, ddof=1)  # ddof = 1 : unbiased estimator => divide std by n-1
    cvnn = sdnn / mean_nni
    
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
        'cvnn': cvnn
    }
    
    return time_domain_features


# TO DO -> GEOMETRICAl FEATURES

# ----------------- FREQUENCY DOMAIN FEATURES ----------------- #

def get_frequency_domain_features(nn_intervals, method="Welch", sampling_frequency=7, interpolation_method="linear",
                                  vlf_band=(0.0033, 0.04), lf_band=(0.04, 0.15), hf_band=(0.15, 0.40), plot=0):
    
    """
    Function returning a dictionnary containing frequency domain features for HRV analyses.
    Must use this function on short term recordings, from 2 to 5 minutes window.
    
    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    method - Method used to calculate the psd. Choice are Welch's FFT or Lomb method.
    sampling_frequency - frequence at which the signal is sampled. Common value range from 
    1 Hz to 10 Hz, by default set to 7 Hz. No need to specify if Lomb method is used.
    interpolation_method - kind of interpolation as a string, by default "linear". No need to 
    specify if lomb method is used
    vlf_band - Very low frequency band for features extraction from power spectral density
    lf_band - Low frequency band for features extraction from power spectral density
    hf_band - High frequency band for features extraction from power spectral density

    Returns
    ---------
    frequency_domain_features - dictionnary containing frequency domain features for HRV analyses.
    There are details about each features given in "get_features_from_psd" function.

    """
    timestamps = create_time_info(nn_intervals)

    # ---------- Interpolation of signal ---------- #
    if method == "Welch":
        funct = interpolate.interp1d(x=timestamps, y=nn_intervals, kind=interpolation_method)

        timestamps_interpolation = create_interpolation_time(nn_intervals, sampling_frequency)
        nni_interpolation = funct(timestamps_interpolation)
        
        # ---------- Remove DC Component ---------- #
        nni_normalized = nni_interpolation - np.mean(nni_interpolation)
    
    #  ----------  Calcul Power Spectral Density  ---------- #
    # Describes the distribution of power into frequency components composing that signal.
        freq, psd = signal.welch(x=nni_normalized, fs=sampling_frequency, window='hann')
    
    elif method == "Lomb":
        freq, psd = LombScargle(timestamps, nn_intervals, normalization='psd').autopower(minimum_frequency=vlf_band[0],
                                                                                         maximum_frequency=hf_band[1])
    else:
        raise ValueError("Not a valid method. Choose between 'Lomb' and 'Welch'")
        
    # ----------  Calcul features  ---------- #
    freqency_domain_features = get_features_from_psd(freq=freq, psd=psd,
                                                     vlf_band=vlf_band,
                                                     lf_band=lf_band,
                                                     hf_band=hf_band)
    
    # TO DO 
    # Plotting options
    if plot == 1:
        pass
    
    return freqency_domain_features


def create_time_info(nn_intervals):
    """
    Function creating time interval of all nn_intervals
    
    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    
    Returns
    ---------
    nni_tmstp - time interval between first NN Interval and final NN Interval
    """
    # Convert in seconds
    nni_tmstp = np.cumsum(nn_intervals) / 1000.0 
    
    # Force to start at 0
    return nni_tmstp - nni_tmstp[0]


def create_interpolation_time(nn_intervals, sampling_frequency=7):
    """
    Function creating the interpolation time used for Fourier transform's method
    
    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    sampling_frequency - frequence at which the signal is sampled
    
    Returns
    ---------
    nni_interpolation_tmstp - timestamp for interpolation
    """
    time_nni = create_time_info(nn_intervals)
    # Create timestamp for interpolation
    nni_interpolation_tmstp = np.arange(0, time_nni[-1], 1 / float(sampling_frequency))
    return nni_interpolation_tmstp


def get_features_from_psd(freq, psd, vlf_band=(0, 0.04), lf_band=(0.04, 0.15), hf_band=(0.15, 0.40)):
    """
    Function computing frequency domain features from the power spectral decomposition.
    
    Arguments
    ---------
    freq - Array of sample frequencies
    psd - Power spectral density or power spectrum
    
    Returns
    ---------
    freqency_domain_features - dictionnary containing frequency domain features for HRV analyses.
    There are details about each features given below.

    Notes
    ---------
    Details about feature engineering...
    
    - **total_power** : Total power density spectra
    
    - **vlf** : variance ( = power ) in HRV in the Very low Frequency (.003 to .04 Hz by default).
    Reflect an intrinsic rhythm produced by the heart which is modulated primarily by sympathetic activity.
    
    - **lf** : variance ( = power ) in HRV in the low Frequency (.04 to .15 Hz). Reflects a  mixture of sympathetic
    and parasympathetic activity, but in long-term recordings like ours, it reflects sympathetic activity and can be
    reduced by the beta-adrenergic antagonist propanolol
    
    - **hf**: variance ( = power ) in HRV in the High Frequency (.15 to .40 Hz by default).
    Reflects fast changes in beat-to-beat variability due to parasympathetic (vagal) activity. 
    Sometimes called the respiratory band because it corresponds to HRV changes related to the
    respiratory cycle and can be increased by slow, deep breathing (about 6 or 7 breaths per 
    minute) and decreased by anticholinergic drugs or vagal blockade.
    
    - **lf_hf_ratio** : lf/hf ratio is sometimes used by some investigators as a quantitative mirror 
    of the sympatho/vagal balance.
    
    - **lfnu** : normalized lf power
    
    - **hfnu** : normalized hf power
    
    References
    ----------
    - Signal Processing Methods for Heart Rate Variability - Gari D. Clifford, 2002
    - Heart rate variability - Standards of measurement, physiological interpretation, and clinical use, Task Force of
    The European Society of Cardiology and The North American Society of Pacing and Electrophysiology, 1996

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

    lf_hr_ratio = lf / hf
    lfnu = (lf / (lf + hf)) * 100
    hfnu = (hf / (lf + hf)) * 100

    freqency_domain_features = {
        'lf': lf,
        'hf': hf,
        'lf_hr_ratio': lf_hr_ratio,
        'lfnu': lfnu,
        'hfnu': hfnu,
        'total_power': total_power,
        'vlf': vlf
    }
    
    return freqency_domain_features

# ----------------- NON lINEAR DOMAIN FEATURES ----------------- #


def get_csi_cvi_features(nn_intervals):
    """
    Function returning a dictionnary containing 3 features from non linear domain  
    for hrV analyses.
    Must use this function on short term recordings, for 30 , 50, 100 Rr Intervals 
    or seconds window.

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Intervals

    Returns
    ---------
    csi_cvi_features - dictionnary containing non linear domain features 
    for hrV analyses. Thera are details about each features are given below. 
    
    Notes
    ---------
    - **csi** : 
    - **cvi** : 
    - **Modified_csi** : 
    
    References
    ----------
    - Using Lorenz plot and Cardiac Sympathetic Index of heart rate variability for detecting seizures for patients
    with epilepsy, Jesper Jeppesen et al, 2014
    """

    # Measures the width and length of poincare cloud
    poincare_plot_features = get_poincare_plot_features(nn_intervals)
    T = 4 * poincare_plot_features['sd1']
    L = 4 * poincare_plot_features['sd2']
    
    csi = L / T
    cvi = np.log10(L * T)
    modified_csi = L**2 / T
    
    csi_cvi_features = {
        'csi': csi,
        'cvi': cvi,
        'Modified_csi': modified_csi
    }
    
    return csi_cvi_features


def get_poincare_plot_features(nn_intervals):
    """
    Function returning a dictionnary containing 3 features from non linear domain
    for hrV analyses.
    Must use this function on short term recordings, from 5 minutes window.

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval

    Returns
    ---------
    poincare_plot_features - dictionnary containing non linear domain features 
    for hrV analyses. Thera are details about each features are given below. 
    
    Notes
    ---------
    - **sd1** : 
    - **sd2** : 
    - **ratio_sd1_sd2** : 
    
    References
    ----------
    - Pre-ictal heart rate variability assessment of epileptic seizures by means of linear and non- linear analyses,
    Soroor Behbahani, Nader Jafarnia Dabanloo et al - 2013

    """
    diff_nn_intervals = np.diff(nn_intervals)
    # measures the width of poincare cloud
    sd1 = np.sqrt(np.std(diff_nn_intervals, ddof=1) ** 2 * 0.5)
    # measures the length of the poincare cloud
    sd2 = np.sqrt(2 * np.std(nn_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_nn_intervals, ddof=1) ** 2)
    ratio_sd1_sd2 = sd1 / sd2
    
    poincare_plot_features = {
        'sd1': sd1,
        'sd2': sd2,
        'ratio_sd1_sd2': ratio_sd1_sd2
    }

    return poincare_plot_features


def get_sampen(nn_intervals):
    """
    Function computing the sample entropy of the given data.
    Must use this function on short term recordings, from 1 minute window.

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval

    Returns
    ---------
    sampen - the sample entropy of the data

    References
    ----------
    - Physiological time-series analysis using approximate entropy and sample entropy,
    JOSHUA S. RICHMAN1, J. RANDALL MOORMAN - 2000
    """

    sampen = nolds.sampen(nn_intervals, emb_dim=2)
    return {'sampen': sampen}
