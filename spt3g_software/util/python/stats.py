"""
Contains useful statistical functions, like gaussian moments, median absolute
deviation, outlier filtering, etc.
"""
import numpy as np
from spt3g import core

@core.usefulfunc
def gaussian2d_moments(data):
    '''Returns (height, amp, x_0, y_0, sigma_x, sigma_y, theta), the gaussian
       parameters of a 2D distribution by calculating its moments.
    '''
    total = np.abs(data).sum()
    Y,X = np.indices(data.shape)
    y_0 = np.argmax((X*np.abs(data)).sum(axis=1)/total)
    x_0 = np.argmax((Y*np.abs(data)).sum(axis=0)/total)

    col = data[int(y_0),:]
    sigma_x = np.sqrt(np.abs((np.arange(col.size)-y_0)*col).sum() /
                      np.abs(col).sum())

    row = data[:, int(x_0)]
    sigma_y = np.sqrt(np.abs((np.arange(row.size)-x_0)*row).sum() /
                      np.abs(row).sum())

    height = np.median(data.ravel())
    if np.abs(data.max()) > np.abs(data.min()):
        amp = data.max() - height
    else:
        amp = data.min() - height
       
    if ( np.isnan(sigma_x) or np.isnan(sigma_y) or
        np.isnan(height) or np.isnan(amp) ):
        raise ValueError('Something is nan...')

    return [amp,x_0,y_0,sigma_x,sigma_y,0.,height]

@core.usefulfunc
def where_not_outlier(data, nsigma=2., frac_cut=0.05, max_iter=15,
                    target_stddev=0.):
    """
    Returns an array containing indices of elements of the input data
    array which are not outliers.. We use an iterative method to
    determine outliers; we recalculate the standard deviation after each
    iteration and check if any new points are outliers.
        
    Intended to replace IDL procedure where_not_outlier.pro.
    
    Parameters:
    -----------
    nsigma [2.0]: (float) A value must be within this many standard deviations
        of the data's median in order to not be an outlier.

    frac_cut [0.05]: (float) End the iteration if the standard deviation
        decreases by this fraction or less of the previous standard deviation.

    max_iter [15]: (int) Halt iteration after this many iterations.

    target_stddev [0.]: (float) Stop iteration if the standard deviation is 
        ever less than this target number.
          
    Returns:
    --------
    An array containing indices of non-outlier data.
    """
    # Start by defining non-outliers as all data which is finite.
    good = np.where(np.isfinite(data))[0]
    
    # Initialize variables for the loop
    n_iterations=0
    good_data = data[good]
    this_stddev, new_stddev = np.std(good_data), 0.
    while n_iterations < max_iter and this_stddev > target_stddev:
        if len(good)<2:
            raise ArithmeticError("Fewer than 2 good points remaining! "
                                  "These data are pathological.")
        
        # Remove indices of data points which are too far from the median.
        good_data -= np.median(good_data, False)
        good = good[(abs(good_data) <= nsigma*this_stddev)]

        # Recompute the standard deviation and check if the fractional
        # change was small enough to let us leave the loop.
        good_data = data[good]
        new_stddev = np.std(good_data)
        
        if (this_stddev - new_stddev)/this_stddev < frac_cut:
            break
        
        # Prepare for the next iteration.
        this_stddev = new_stddev
        n_iterations += 1

    return good

@core.usefulfunc
def float_mode(x, bin_size=None, n_bins=None):
    """
    Finds the mode of an array of values. Intended to be used with floats.

    Parameters:
    -----------
    x : The NumPy array of values.
        It is treated as a 1D array, regardless of dimension.

    bin_size [None]: If supplied, the size of a bin within which numbers are
        considered to be the same value.

    n_bins [None]: If supplied, the number of bins into which the input
        numbers will be sorted. Overrides any bin_size input. If neither
        bin_size nor n_bins are supplied, n_bins will be taken to be 100.

    Returns:
    --------
    The most common number in the input array, as a float.
    This will be the center of the bin with the most entries.
    """

    # First, make sure that we have a 1-D array.
    flat_arr = x.flatten()

    # Next, find the number of bins to use. Default to 100 if we
    # weren't given a bin size or number of bins.
    if n_bins is not None:
        pass
    elif bin_size is not None:
        n_bins = int( float(flat_arr.max() - flat_arr.min()) / bin_size ) + 1
    else:
        n_bins = 100

    # Histogram the input values.
    hist, bins = np.histogram(flat_arr, bins=n_bins)

    # Return the center of the bin with the most entries.
    mode = (bins[np.argmax(hist)] + bins[np.argmax(hist)+1])/2.

    return mode

@core.usefulfunc
def weighted_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    
    AUTHOR
        Stephen Hoover, 12 November 2013 (but stolen from the web)
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights) 
    return np.sqrt(variance)

@core.usefulfunc
def robust_sigma(array, nsthresh=3.):

    """
    Returns standard deviation of an array, iteratively throwing out > N*sigma
    outliers. Default is N = 3.
    """

    ftol = 1e-4
    npts = np.size(np.asarray(array))
    fnpts = np.float(npts)
    med0 = np.median(np.asarray(array))
    sd0 = np.sqrt(np.sum((np.asarray(array)-med0)**2)/(fnpts-1.))
    sd = sd0
    md = med0
    array_cut = np.asarray(array).copy()
    tol = 1e6
    while tol > ftol:
        array_cut = array_cut[np.abs(array_cut - md) <= nsthresh * sd]
        if len(array_cut) == 0:
            break
        newmd = np.mean(array_cut)
        newsd = np.std(array_cut)
        tol = np.abs(newsd - sd)/sd
        sd = newsd
        md = newmd

    return sd

