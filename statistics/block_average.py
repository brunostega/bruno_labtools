import numpy as np 

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