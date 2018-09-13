import matplotlib.pyplot as plt
from matplotlib import style


def plot_timeseries(nn_intervals):
    """
    Function plotting the NN Intervals timeseries

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    """

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Rr Interval through time")
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Rr Interval", fontsize=15)
    plt.plot(nn_intervals)
    plt.show()


def plot_distrib(nn_intervals, bin_length=8):
    """
    Function plotting histogram distribution of the NN Intervals. Useful for geometrical features

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    bin_length - size of the bin for histogram in ms, by default = 8
    """

    max_nn_i = max(nn_intervals)
    min_nn_i = min(nn_intervals)

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Density of Rr Intervals")
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Number of Rr Interval per bin", fontsize=15)
    plt.hist(nn_intervals, bins=range(min_nn_i - 10, max_nn_i + 10, bin_length))
    plt.show()


def plot_psd(nn_intervals):
    """
    Function plotting the power spectral density of the NN Intervals

    Arguments
    ---------
    nn_intervals - list of Normal to Normal Interval
    """

    # TO DO
    pass
    return None


def plot_poincare(nn_intervals):
    """
    Pointcare / Lorentz Plot of the NN Intervals

    Arguments
    ---------
    nn_intervals - list of NN intervals

    Notes
    ---------
    The transverse axis (T) reflects beat-to-beat variation
    the longitudinal axis (L) reflects the overall fluctuation
    """

    ax1 = nn_intervals[:-1]
    ax2 = nn_intervals[1:]

    style.use('ggplot')
    plt.figure(figsize=(12, 8))
    plt.title("Pointcare / Lorentz Plot", )
    plt.xlabel('NN_n (s)', fontsize=15)
    plt.ylabel('NN_n+1 (s)', fontsize=15)
    plt.scatter(ax1, ax2, c='r', s=12)
    plt.show()
