# Heart Rate Variability analysis

**hrvanalysis** is a Python module for Heart Rate Variability analysis of RrIntervals built on top of SciPy, AstroPy, Nolds and NumPy and distributed under the GPLv3 license.

The development of this library started in July 2018 as part of Aura Healthcare project and is maintained by Robin Champseix.

Website : https://www.aura.healthcare 

Github : https://github.com/Aura-healthcare

version : 0.0.1


### Installion / Prerequisites

#### Dependencies

hrvanalysis requires the following:
- Python (>= 3.4)
- astropy = 3.0.4
- future = 0.16.0
- nolds = 0.4.1
- numpy = 1.15.1
- scipy = 1.1.0

#### User installation

The easiest way to install hrvanalysis is using ``pip`` :

    $ pip install -U hrvanalysis

you can also clone the repository:

    $ git clone https://github.com/robinchampseix/hrvanalysis.git
    $ python setup.py install

### Getting started 

There are 3 types of features you can get from NN Intervals: 

> Time domain features : **Mean_NN, SDNN, SDSD, NN50, pNN50, NN20, pNN20, RMSSD, Median_NN, Range_NN**

> Frequency domain features : **LF, HF, VLF, LH/HF ratio, LFnu, HFnu, Total_Power**

> Non Linear domain features : **CSI, CVI, Modified_CSI, SD1, SD2, SD1/SD2 ratio, SampEn**

As an exemple, what you can compute to get Time domain analysis is :

```python
from hrvanalysis import get_time_domain_features
nn_intervals = 
time_domain_features = get_time_domain_features(nn_intervals)
>>> time_domain_features
{'mean_nn': 718.248,
 'sdnn': 43.113074968427306,
 'sdsd': 19.519367520775713,
 'nn50': 24,
 'pnn50': 2.4,
 'nn20': 225,
 'pnn20': 22.5,
 'rmssd': 19.519400785039664,
 'Median_nn': 722.5,
 'Range_nn': 249}
```

You can find how to use methods, references and details about each feature in the documentation of each function:
- get_time_domain_features
- get_frequency_domain_features
- get_csi_cvi_features
- get_poincare_plot_features
- get_sampen

## References

Here are the main references used to compute the set of features and for signal processing methods:

> Heart rate variability - Standards of measurement, physiological interpretation, and clinical use, Task Force of The European Society of Cardiology and The North American Society of Pacing and Electrophysiology, 1996
    
> Signal Processing Methods for Heart Rate Variability - Gari D. Clifford, 2002

> Physiological time-series analysis using approximate entropy and sample entropy, Joshua S. Richman, J. Randall Moorman - 2000
    
> Using Lorenz plot and Cardiac Sympathetic Index of heart rate variability for detecting seizures for patients with epilepsy, Jesper Jeppesen et al, 2014

## Authors

**Robin Champseix** - *Initial work* - (https://github.com/robinchampseix)

## License

This project is licensed under the *GNU GENERAL PUBLIC License* - see the [LICENSE.md](https://github.com/robinchampseix/hrv_library/LICENSE) file for details

## Acknowledgments

I hereby thank Laurent Ribière and Clément Le Couedic, my coworkers who gave me time to Open Source this library.
I also thank Fabien Arcellier for his advices on to how build a library in PyPi.