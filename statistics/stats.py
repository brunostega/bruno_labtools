import numpy as np 
import scipy as sp

def autocorrelation_function(data):
    """
    Calculate the autocorrelation of a given dataset.

    Parameters:
    data (array-like): A one-dimensional array or list of numerical data.

    Returns:
    numpy.ndarray: An array of autocorrelation values for different time lags,
                   normalized by the variance of the input data.

    The function calculates the autocorrelation for lags ranging from 0 to
    half the length of the input data. Autocorrelation measures how the data
    points in the series are related to each other at different time lags.
    """
    Ac = []
    mean = np.mean(data)

    for i in range(int(len(data)/2)):
        Ac_appo = np.mean((data-mean)*(np.roll(data,i)-mean))
        Ac.append(Ac_appo)
    return Ac/np.var(data)


def block_analysis(data, block_size_step=2):
    """
    Calculate the standard deviation of block means for different block sizes.

    Parameters:
    data (array-like): A one-dimensional array or list of numerical data.
    block_size_step (int): Step size to increase the block size.

    Returns:
    tuple: A tuple containing two numpy arrays:
           - sigma_blocks: Standard deviations of the block means for different block sizes.
           - block_sizes: Corresponding block sizes.

    The function computes the standard deviation of block means for block sizes
    ranging from 2 to one-fifth of the length of the input data, increasing by
    the specified block size step. The block size is the number of data points
    in each block.
    """

    bs_sizes = np.arange(2,len(data)//5, block_size_step)
    blocks = np.array(len(data)/ bs_sizes, dtype = int)

    data_size = int(len(data))
    sigma_blocks = []

    for b in blocks:
        sigma_blocks.append( np.std(np.mean( np.array(np.split( data[:int(data_size/b)*b], b) ) , axis=1 ))/(np.sqrt(b)))

    sigma_blocks = np.array(sigma_blocks)
    block_sizes = np.array(len(data)/blocks, dtype = int)

    return blocks,block_sizes, np.mean(data), sigma_blocks

def cdf_poisson(t,tau):
    """
    Cumulative distribution function (CDF) for the Poisson distribution.

    Parameters:
    t (array-like): A one-dimensional array of time values.
    tau (float): The mean time between events (tau) for the Poisson distribution.

    Returns:
    array: The Poisson cumulative distribution evaluated at the given time values.
    """
    return 1 - np.exp(-t/tau)

def poisson_fitting(times):
    """
    Fit a Poisson distribution to the given times and perform a Kolmogorov-Smirnov test.

    Parameters:
    times (array-like): A one-dimensional array or list of times (event occurrences).

    Returns:
    tuple: A tuple containing:
           - popt (array): Optimal parameter for the Poisson distribution fit.
           - ks_test (tuple): Results of the Kolmogorov-Smirnov test comparing
                              the empirical cumulative distribution function (CDF)
                              with the fitted Poisson CDF.

    The function performs the following steps:
    1. Calculates the empirical cumulative distribution of the given times.
    2. Fits the empirical data to a Poisson cumulative distribution function (CDF).
    3. Conducts a Kolmogorov-Smirnov test to compare the empirical CDF with the fitted Poisson CDF.
    """

    #calculate cumulative 
    ordered_times = np.sort(times)
    cum = np.array([ i / (len(ordered_times) - 1) for i in range(len(ordered_times)) ])

    #fit to poisson distribution
    popt, pcov = sp.optimize.curve_fit(cdf_poisson, ordered_times, cum, p0 = np.mean(ordered_times))
    fit_data = cdf_poisson(ordered_times, popt[0])

    #kolmogorov-smirnov test
    ks_test = sp.stats.kstest(cum, fit_data)

    return popt[0], ks_test


