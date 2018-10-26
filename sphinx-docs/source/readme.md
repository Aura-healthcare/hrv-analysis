# Getting started

## hrv-analysis

[![PyPI version](https://badge.fury.io/py/hrv-analysis.svg)](https://badge.fury.io/py/hrv-analysis)
[![Build Status](https://travis-ci.com/robinchampseix/hrvanalysis.svg?branch=master)](https://travis-ci.com/robinchampseix/hrvanalysis)
[![codecov](https://codecov.io/gh/robinchampseix/hrvanalysis/branch/master/graphs/badge.svg)](https://codecov.io/gh/robinchampseix/hrvanalysis)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**hrvanalysis** is a Python module for Heart Rate Variability analysis of RR-intervals built on top of SciPy, AstroPy, Nolds and NumPy and distributed under the GPLv3 license.

The development of this library started in July 2018 as part of Aura Healthcare project and is currently maintained by Robin Champseix.

![Plot image](timeserie_distrib_plot.png)

## Installation / Prerequisites

#### Dependencies

hrvanalysis requires the following:
- Python (>= 3.6)
- astropy >= 3.0.4
- future >= 0.16.0
- nolds >= 0.4.1
- numpy >= 1.15.1
- scipy >= 1.1.0


#### User installation

The easiest way to install hrv-analysis is using ``pip`` :

    $ pip install hrv-analysis

you can also clone the repository:

    $ git clone https://github.com/robinchampseix/hrvanalysis.git
    $ python setup.py install


## Quick examples

### Outliers and ectopic beats filtering methods

This package provides methods to remove outliers and ectopic beats from signal for further analysis. Those methods are useful to get Normal to Normal Interval (NN interval) from Rr Interval.
Please use this methods carefully as they might have a huge impact on features calculation.

```python
from hrvanalysis.clean_outliers import clean_outlier, clean_ectopic_beats

# rr_intervals is a list containing integers value of Rr Intervals
cleaned_rr_intervals = clean_outlier(rr_intervals=rr_intervals,  low_rri=300, high_rri=2000) # This remove outliers from signal

nn_interval = clean_ectopic_beats(rr_intervals=cleaned_rr_intervals, method="Malik") # This remove ectopic beats from signal
```

You can find how to use the following methods, references and more details in the Tutorial part (*https://robinchampseix.github.io/hrvanalysis/tutorial.html*):
- clean_outlier
- clean_ectopic_beats

### Features calculation 

There are 4 types of features you can get from NN Intervals: 

> Time domain features : **Mean_NNI, SDNN, SDSD, NN50, pNN50, NN20, pNN20, RMSSD, Median_NN, Range_NN, CVSD, CV_NNI, Mean_HR, Max_HR, Min_HR, STD_HR**

> Geometrical domain features : **Triangular_index, TINN**

> Frequency domain features : **LF, HF, VLF, LH/HF ratio, LFnu, HFnu, Total_Power**

> Non Linear domain features : **CSI, CVI, Modified_CSI, SD1, SD2, SD1/SD2 ratio, SampEn**

As an exemple, what you can compute to get Time domain analysis is :

```python
from hrvanalysis.extract_features import get_time_domain_features
 
# nn_intervals is a list containing integers value of NN Intervals
time_domain_features = get_time_domain_features(nn_intervals)

>>> time_domain_features
{'mean_nni': 718.248,
 'sdnn': 43.113,
 'sdsd': 19.519,
 'nni_50': 24,
 'pnni_50': 2.4,
 'nni_20': 225,
 'pnni_20': 22.5,
 'rmssd': 19.519,
 'median_nni': 722.5,
 'range_nni': 249,
 'cvsd': 0.0272,
 'cvnni': 0.060,
 'mean_hr': 83.847,
 'max_hr': 101.694,
 'min_hr': 71.513,
 'std_hr': 5.196}
```

You can find how to use the following methods, references and more details in the Tutorial part (*https://robinchampseix.github.io/hrvanalysis/tutorial.html*):
- get_time_domain_features
- get_geometrical_features
- get_frequency_domain_features
- get_csi_cvi_features
- get_poincare_plot_features
- get_sampen

### Plot functions

There are several plot functions that allow you to see, for example, the Power spectral density for frequency domain features :

```python
from hrvanalysis.plot import plot_psd, plot_distrib

# nn_intervals is a list containing integers value of NN Intervals
plot_psd(nn_intervals, method="Welch")
plot_distrib(nn_intervals)
```

![alt text](https://github.com/robinchampseix/hrvanalysis/blob/master/figures/lomb_density_plot.png)

You can find how to use the following methods, references and more details in the Tutorial part (*https://robinchampseix.github.io/hrvanalysis/tutorial.html*):
- plot_distrib
- plot_timeseries
- plot_psd
- plot_poincare
