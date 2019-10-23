"""
Contains various useful math functions which cannot be found
in other packages. These are generic tools which do not use any 3G specific
data structures.
"""
import numpy as np
import numexpr as ne
import scipy
from scipy.optimize import curve_fit, fmin_slsqp, leastsq
from spt3g import core
from spt3g.util.genericutils import shift

@core.usefulfunc
def gaussian2d(amp, y_0, x_0, sigma_y, sigma_x, theta, height=0.):
    '''Returns a 2D gaussian with the given parameters.'''
    sigma_x = np.float(sigma_x)
    sigma_y = np.float(sigma_y)

    a = np.cos(theta)**2./(2.*sigma_x**2.) + np.sin(theta)**2./(2.*sigma_y**2.)
    b = np.sin(2.*theta)/(4.*sigma_x**2.) - np.sin(2.*theta)/(4.*sigma_y**2.)
    c = np.sin(theta)**2./(2.*sigma_x**2.) + np.cos(theta)**2./(2.*sigma_y**2.)
    return lambda x,y: height + amp*np.exp(-(a*(x-x_0)**2. +
                                             2.*b*(x-x_0)*(y-y_0) +
                                             c*(y-y_0)**2.))

def gauss2d_deriv(x, constant=1., mean_x=0., sigma_x=1., mean_y=0., sigma_y=1.,
                  offset=0., rotation=0.):
    """
    Partial derivative of the 2D Gaussian.
    """
    try:
        if x.shape[1]==2: y, x = x[...,0], x[...,1]
        else: y, x = x[0], x[1]
    except AttributeError:        
        y, x = x[0], x[1]
        
    xp = (x-mean_x)*np.cos(rotation) - (y-mean_y)*np.sin(rotation)
    yp = (x-mean_x)*np.sin(rotation) + (y-mean_y)*np.cos(rotation)
    exp_term = (np.exp( -( (xp/sigma_x)**2 + (yp/sigma_y)**2) / 2 ) /
                (2*np.pi*sigma_x*sigma_y))
    doffset = [1.0]*len(x)
    dconstant = exp_term
    
    exp_term *= constant
    dmean_x = exp_term*(xp*np.cos(rotation)/sigma_y + 
                        yp*np.sin(rotation)/sigma_x)
    dmean_y = exp_term*(-xp*np.sin(rotation)/sigma_y + 
                        yp*np.cos(rotation)/sigma_x)
    dsigma_x = exp_term * yp**2 / sigma_x
    dsigma_y = exp_term * xp**2 / sigma_y
    drotation = -exp_term * xp*yp * (sigma_y/sigma_x - sigma_x/sigma_y)

    return dmean_x, dsigma_x, dmean_y, dsigma_y, dconstant, doffset, drotation

def find_clustering(_input, fof_radius):
    """
    Takes an input list of points and uses a simple friends-of-friends
    algorithm to sort them into individual clusters of points.
    
    Not the most efficient implementation. Could be speeded up if necessary.
    
    Parameters:
    -----------
    _input : (Iterable of 2-element iterables)
        A list or array of coordinates of points in the plane.

    fof_radius : (float) If two points are closer than this (in whatever
        coordinate system the points are in), then combine them into the
        same "cluster" of points.
            
    Returns:
    --------
    A three-tuple:  cluster_center (a list of points which are the mean of 
                                    all points in the cluster)
                    clusters (a list of 2-D arrays which contain all the
                              points in a cluster)
                    n_cluster_points (a list of integer cluster sizes)
    """
    # Filter out the NaNs from the input array. This also copies the input
    # so that it won't be changed.
    unassigned_points = np.ascontiguousarray(
        [ pt for pt in _input if np.isfinite(pt).all()])

    # Set up the output lists.
    clusters = []
    cluster_center = []
    n_cluster_points = []
    
    # Go through the points and assign them to clusters. When a cluster is 
    # complete, the next remaining point nucleates a new cluster. Once we've    
    # started a cluster, go through each point in the input list one at a time.
    # Compute its distance to each cluster member. If it's close enough to any
    # of the cluster members, add it to the cluster. If we've added any
    # point(s) to the cluster, go through remaining points once more, to see
    # if the new point(s) are close enough to them.
    while len(unassigned_points)>0:
        seed = unassigned_points[0]
        unassigned_points = np.delete(unassigned_points, 0, axis=0)
        
        clusters.append(np.ascontiguousarray( [ seed ]))

        cluster_has_grown = True
        while cluster_has_grown:
            i_pts_assigned = []
            for i_pt, pt in enumerate(unassigned_points):
                dist = np.sqrt((clusters[-1][:,0]-pt[0])**2 +
                               (clusters[-1][:,1]-pt[1])**2)
                if np.min(dist) < fof_radius:
                    clusters[-1] = np.append(
                        clusters[-1], np.ascontiguousarray([pt]), axis=0)
                    i_pts_assigned.append(i_pt)
            if len(i_pts_assigned)!=0:
                unassigned_points = np.delete(
                    unassigned_points, np.ascontiguousarray(i_pts_assigned),
                    axis=0)
            cluster_has_grown = (len(i_pts_assigned)!=0)
            
        cluster_center.append( np.mean(clusters[-1], axis=0))
        n_cluster_points.append(len(clusters[-1]))

    return cluster_center, clusters, n_cluster_points

@core.usefulfunc
def distance_array(nside, radius=1.):
    """
    Returns an nside[0] x nside[1] array in which each element is the distance 
    from the center of the array, normalized such that the length of each
    side of the array is radius[0], radius[1] units.
    
    Parameters:
    -----------
    nside: (int or 2-tuple of ints) 
        The number of pixels on each side of the array.
        If the input is a scalar, the array will be square.

    radius [1.]: (float or 2-tuple of floats) 
        The half-length of each side of the array.
        If the input is a scalar, each side will have the same length.
        
    Returns:
    --------
    An np.ndarray with shape [nside[0],nside[1]] where each element is
        the distance from the center.
        
    AUTHOR
        Stephen Hoover, 17 May 2013
    """
    try:
        nside_y, nside_x = nside
    except TypeError:
        nside_y = nside_x = nside
        
    try:
        radius_y, radius_x = radius
    except TypeError:
        radius_y = radius_x = radius
    
    ys, xs = np.meshgrid( np.linspace(-radius_x, radius_x, nside_x),
                         np.linspace(-radius_y, radius_y, nside_y) )
    center_dist = np.sqrt(xs**2 + ys**2)

    return center_dist

@core.usefulfunc
def sum_circle2d(map, center=None, r_cutoff=None):
    """
    Sums all values in a 2D array up to r_cutoff indices away from
    a specified center.

    Parameters:
    -----------
    map: 2D array to perform circular sum over.

    center [None]: [y,x] position of center of circle.  Defaults to
                center of map.

    r_cutoff [None]: Radius of circle.  Defaults to just returning the
                   value at center.

    Returns:
    --------
    sum of all map pixels within r_cutoff
    """
    if r_cutoff is None:
        #if not specified, r_cutoff is less than one index
        #so that sumCircle2D returns only value at center
        r_cutoff = 0.9

    if center is None:
        #if not specified, default to center of array
        center = np.ascontiguousarray([np.floor(map.shape[0]/2),
                                       np.floor(map.shape[1]/2)])
        
    xs, ys = np.meshgrid(
        np.ascontiguousarray(np.arange(-center[0], map.shape[0]-center[0])), 
        np.ascontiguousarray(np.arange(-center[1], map.shape[1]-center[1])) )
    dr = np.sqrt(xs**2 + ys**2)
    xxx = np.sum(map[dr <= r_cutoff])
    
    return xxx

@core.usefulfunc
def recenter_map(map0, xcenter, ycenter, widthx=None, widthy=None):
    """
    Move the center of map0 to (xcenter, ycenter) by cutting pixels 
    at the edges or, if widthx or widthy is not None and is greater
    than the width resulting from the cut, zero pad the edge(s). 
    If widthx/y is None, cut the minimum possible number of pixels.
    
    Note that (xcenter, ycenter) are in numpy coordinates-- if you 
    plot map0 with imshow, the x-axis is vertical and y-axis is 
    horizontal, so are the opposite of numpy array conventions.

    Required Options:
    ----------------
    map0 : 2-d array
        Array of data to recenter
    xcenter : int
        x location at which to center array
    ycenter : int
        y location at which to center array

    Optional Options:
    ----------------
    widthx : int
        Desired width in x of output array. If longer than there
        are pixels in the recentered map, zero-pad. If None, cut
        as few pixels as possible.
    widthy : int
        Desired width in y of output array. If longer than there
        are pixels in the recentered map, zero-pad. If None, cut
        as few pixels as possible.
    
    Returns
    -------
    map_out : 2-d array
        New array centered at (x, y)
    """
    xcenter = int(xcenter)
    ycenter = int(ycenter)

    map_dims = map0.shape
    if widthx is None:
        widthx = np.min([2 * xcenter, 2 * np.abs(map_dims[0] - xcenter)])
    if widthy is None:
        widthy = np.min([2 * ycenter, 2 * np.abs(map_dims[1] - ycenter)])

    # zero pad to make sure map0 is big enough to slice given desired widths
    map_out = np.pad(np.asarray(map0), ((widthx // 2, widthx // 2), 
                                        (widthy // 2, widthy // 2)),
                     mode='constant',
                     constant_values=((0, 0), (0, 0)))
    map_out = map_out[xcenter:xcenter + widthx, ycenter:ycenter + widthy]
    return map_out

# =============================================================================
# Generic function forms for fits
# -----------------------------------------------------------------------------
def poly(x, *coeffs):
    """
    A polynomial function.
    
    The order of the polynomial is determined by the number
    of coefficients passed to this function as arguments.
    Repeated tries are unable to do better than the speed
    of np.polyval, at least on large inputs.
    """
    return np.polyval(coeffs[::-1], x)

def poly_deriv(x,*coeffs):
    """
    The derivative of polynomial function, used for fitting functions.
    
    The order of the polynomial is determined by the number
    of coefficients passed to this function as arguments.
    There must be at least two coefficients.
    """
    if len(coeffs)==2:
        y = np.empty(len(x))
        y.fill(coeffs[1])
    else:
        y = (2. * x * coeffs[2]) + coeffs[1]
        if (len(coeffs) > 3):
            t = 2.*x.copy()
            for order, coeff in enumerate(coeffs[3:]):
                t *= (order+3)*x
                y += coeff * t
    return y


def elliptic_parabaloid(a, b, x0, y0, offset, rotation, x, y):
    """
    Parameters:
    -----------
    a: Semi-major axis
    b: Semi-minor axis
    x0: x-coordinate center
    y0: y-coordinate center
    offset : Constant offset
    rotation : In degrees, a rotation of the parabaloid around its center.

    x: Array of x-coordinates
    y: Array of y-coordinates
    """
    rotation *= np.pi/180.
    x_rot = np.cos(rotation)*(x-x0) - np.sin(rotation)*(y-y0)
    y_rot = np.sin(rotation)*(x-x0) + np.cos(rotation)*(y-y0)
    return (x_rot/a)**2 + (y_rot/b)**2 + offset


def square_wave(x, frequency, phase=0., amplitude=1.): 
    """
    Returns square wave array between 0 and amplitude
    """
    return amplitude*np.ascontiguousarray(
        np.sin(2*np.pi*frequency * x + phase*np.pi/180.) > 0., dtype=np.float32)
    

# =============================================================================
# Derivatives
# -----------------------------------------------------------------------------

def deriv(x, y=None):
    """
    Take the derivative of y. If only one input is given, then assume that the
        points are equally spaced.

    Converted from the IDL deriv.pro.
    http://idlastro.gsfc.nasa.gov/ftp/ittvislib/deriv.pro
    
    df/dx = y0*(2x-x1-x2)/(x01*x02) + y1*(2x-x0-x2)/(x10*x12) + 
            y2*(2x-x0-x1)/(x20*x21)
    Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
    """
    if y is not None:
        if len(x)!=len(y):
            raise ValueError(
                "The x array must be the same length as the y array!")
        
        x12 = x - np.roll(x,-1)
        x01 = np.roll(x,1) - x
        x02 = np.roll(x,1) - np.roll(x,-1)
        
        # Middle points
        d = (np.roll(y,1) * (x12 / (x01*x02)) + y * (1./x12 - 1./x01)
             - np.roll(y,-1) * (x01 / (x02 * x12)))

        # Formulae for the first and last points:
        d[0] = (y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) -
                y[1] * x02[1]/(x01[1]*x12[1]) +
                y[2] * x01[1]/(x02[1]*x12[1]))
        d[-1] = (-y[-3] * x12[-2]/(x01[-2]*x02[-2]) +
                 y[-2] * x02[-2]/(x01[-2]*x12[-2]) -
                 y[-1] * (x02[-2]+x12[-2]) / (x02[-2]*x12[-2]))
    else: # Assume points are equally spaced.
        d = (np.roll(x,-1) - np.roll(x,1))/2.
        d[0] = (-3.0*x[0] + 4.0*x[1] - x[2])/2.
        d[-1] = (3.*x[-1] - 4.*x[-2] + x[-3])/2.

    return d


def partial_deriv(grid, resolution, axis, range=1):
    """
    Takes the partial derivative of an input array.
    The grid entries are assumed to be equally spaced.
    Uses the central finite difference method.
    
    Parameters:
    -----------
    grid: (ndarray) The array of which we want to take the derivative. Can be
      of any dimensionality.

    resolution: (float) The separation of grid elements in the derivative's
        direction.

    axis: (int or str) Either the number of the axis along which to take the
        derivative (this corresponds to the position of the index), or the
        name of the axis. For a 1D, 2D, or 3D array, the 'x' axis is the last
        index, the 'y' axis (in a 2D or 3D array) is the second-to-last index,
        and the 'z' axis (in a 3D array only) is the first (0th) index.
        Note that positive index numbers (0, 1, 2, ...) count from the leftmost
        index, and negative index numbers (-1, -2, -3, ...) count from the
        rightmost index.

    range [1]: (int) Use grid elements up to this many steps away to form 
        the finite difference.

    Returns:
    --------
    An array equal to d / d(axis) of the input array.
    """

    # Check for string-ness and, if stringy, assign the proper
    # axis number (counting from the back).
    this_axis = axis
    try:
        if axis.lower()=='x': this_axis=-1
        elif axis.lower()=='y': this_axis=-2
        elif axis.lower()=='z': this_axis=-3
        else:
            raise ValueError(
                "I don't recognize "+str(axis)+" as the name of an axis.")
    except AttributeError:
        pass

    if range==1:
        deriv = (np.roll(grid,-1,axis=this_axis) -
                 np.roll(grid,1,axis=this_axis))/(2.*resolution)
    elif range==2:
        deriv = (-np.roll(grid,-2,axis=this_axis) +
                 8*np.roll(grid,-1,axis=this_axis) -
                 8*np.roll(grid,1,axis=this_axis) +
                 np.roll(grid,2,axis=this_axis))/(12.*resolution)
    else:
        raise NotImplementedError
    
    return deriv


def ddx(grid, resolution, range=1):
    """
    Alias for the first order partial derivative along the 'x' axis.
    """
    return partial_deriv(grid, resolution, 'x', range=range)


def ddy(grid, resolution, range=1):
    """
    Alias for the first order partial derivative along the 'y' axis.
    """
    return partial_deriv(grid, resolution, 'y', range=range)


def partial_second_deriv(grid, resolution, axis, range=1):
    """
    Takes the second order partial derivative along a single axis of an input
    array. The grid entries are assumed to be equally spaced.
    Uses the central finite difference method.
    
    Parameters:
    -----------
    grid: (ndarray) The array for which we want to take the derivative. Can be
        of any dimensionality.

    resolution: (float) The separation of grid elements in the derivative's
        direction.

    axis: (int or str) Either the number of the axis along which to take the
        derivative (this corresponds to the position of the index), or the name
        of the axis. For a 1D, 2D, or 3D array, the 'x' axis is the last index,
        the 'y' axis (in a 2D or 3D array) is the second-to-last index,
        and the 'z' axis (in a 3D array only) is the first (0th) index.
        Note that positive index numbers (0, 1, 2, ...) count from the leftmost
        index, and negative index numbers (-1, -2, -3, ...) count from the
        rightmost index.

    range [1]: (int) Use grid elements up to this many steps away to form 
        the finite difference.

    Returns:
    --------
    An array (of the same dimension as "grid") equal to d^2 / d(axis)^2 of
    the input array.
    """
    # Check for string-ness and, if stringy, assign the proper
    # axis number (counting from the back).
    this_axis = axis
    try:
        if axis.lower()=='x': this_axis=-1
        elif axis.lower()=='y': this_axis=-2
        elif axis.lower()=='z': this_axis=-3
        else:
            raise ValueError(
                "I don't recognize "+str(axis)+" as the name of an axis.")
    except AttributeError:
        pass

    if range==1:
        deriv = (np.roll(grid,-1,axis=this_axis) -
                 2*grid + np.roll(grid,1,axis=this_axis))/resolution**2
    elif range==2:
        deriv = (-np.roll(grid,2,axis=this_axis) +
                 16*np.roll(grid,1,axis=this_axis) -
                 30*grid + 16*np.roll(grid,-1,axis=this_axis) -
                 np.roll(grid,-2,axis=this_axis) ) / (12.*resolution**2)
    else:
        raise NotImplementedError

    return deriv


def d2dx2(grid, resolution, range=1):
    """
    Alias for the second order partial derivative along the 'x' axis.
    """
    return partial_second_deriv(grid, resolution, 'x', range=range)


def d2dy2(grid, resolution, range=1):
    """
    Alias for the second order partial derivative along the 'y' axis.
    """
    return partial_second_deriv(grid, resolution, 'y', range=range)


def d2dxdy(array, resolution, range=1):
    """
    Take the mixed second order partial derivative of an input 2D array.
    Uses the central finite difference method to find derivatives.

    Parameters:
    -----------
    array: (ndarray) The array for which we want to take the derivative.
      Must be 2D.

    resolution: (float or 2-element iterable of floats) The resolution of each 
        grid dimension. If a scalar, the same resolution is used for each axis.
        If an iterable, then element 0 is used for the resolution of axis 0,
        and element 1 for axis 1.

    range [1]: (int) Use grid elements up to this many steps away to form 
        the finite difference.

    Returns:
    --------
    An array with the same dimensions as the input "array" which contains the
    partial derivative of the input.
    """

    if np.isscalar(resolution): res = 2*[resolution]
    else: res = resolution

    if range==1:
        return (shift(array,[1,1]) - shift(array,[1,-1]) -
                shift(array,[-1,1]) + shift(array,[-1,-1])) / (4*res[0]*res[1])
    else:
        # For pixel range > 1 do this by taking ddy then ddx. Could replace
        # with the explicit coefficients if / when I bother to work it out.
        return ddx(ddy(array, resolution, range), resolution, range)


# =============================================================================
# Autocorrelation and autocovariance
# -----------------------------------------------------------------------------  
def fft_cyclic_autocovariance(signal):
    """
    Given a n-dimensional signal, return an estimation of its
    autocovariance function.

    The estimation is made by considering that the input signal
    actually desribes a full period of a wider, cyclic signal. The
    estimation is then the autocovariance of this wider signal.

    Uses the Fast Fourier Transform internally.
    
    Code from http://www.tibonihoo.net/literate_musing/autocorrelations.html
    """
    centered_signal = signal - np.mean(signal)
    ft_signal = np.fft.fftn(centered_signal)
    power_spectral_density = np.abs(ft_signal)**2
    autocovariance = np.fft.ifftn(power_spectral_density)/len(centered_signal)
    return np.real(autocovariance)


def fft_cyclic_autocorrelation(signal):
    """
    Given a n-dimensional signal, return an estimation of its
    autocorrelation function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fft_cyclic_autocovariance.
    """
    autocovariance = fft_cyclic_autocovariance(signal)
    variance = autocovariance.flat[0]
    if variance==0.:
        return np.ascontiguousarray(np.zeros(autocovariance.shape))
    else:
        return (autocovariance / variance)


def fft_autocovariance(signal):
    """
    Compute the autocovariance of the input n-dimensional signal.

    Consider the input signal to be a representative sample of a wider
    signal that has no other pattern that those present on the sample
    (this is what "representative" stands for) and especially no
    pattern whose scale is higer or equal to the input signal's size
    on each of its dimensions (this is for the difference with
    fft_cyclic_autocovariance).

    The autocovariance is computed by a FFT and with a zero padding
    made in such a way that the padded signal is `2**n` bigger than
    the input one (where n is the dimension). However the returned
    function is of the same size as the signal on every dimension.
    
    Code from http://www.tibonihoo.net/literate_musing/autocorrelations.html
    """
    centered_signal = signal - np.mean(signal)
    padded_shape = [2*s+1 for s in centered_signal.shape]
    ft_signal = np.fft.fftn(centered_signal, padded_shape)
    
    pseudo_powerSpectralDensity = np.abs(ft_signal)**2
    pseudo_autocovariance = np.fft.ifftn(pseudo_powerSpectralDensity)
    
    input_domain  = np.ones_like(centered_signal)
    ft_mask = np.fft.fftn(input_domain, padded_shape)
    mask_correction_factors = np.fft.ifftn(np.abs(ft_mask)**2)
    autocovariance = pseudo_autocovariance / mask_correction_factors
    
    crop_slices = [slice(i) for i in signal.shape]
    return np.real(autocovariance[crop_slices])


def fft_autocorrelation(signal):
    """
    Given a n-dimensional signal, return an estimation of its
    autocorrelation function.

    The autocorrelation is obtained by normalizing the autocovariance
    function computed by fft_autocovariance.
    
    Code from http://www.tibonihoo.net/literate_musing/autocorrelations.html
    """
    autocovariance = fft_autocovariance(signal)
    variance = autocovariance.flat[0]
    if variance==0.:
        return np.ascontiguousarray(np.zeros(autocovariance.shape))
    else:
        return (autocovariance / variance)


# =============================================================================
# Smoothing functions
# -----------------------------------------------------------------------------

def gauss_kernel(size, sizey=None):
    """ 
    Returns a normalized 2D gauss kernel array for convolutions 
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()


def blur_image(im, n, ny=None) :
    """ 
    blurs the image by convolving with a gaussian kernel of typical
    size n. The optional keyword argument ny allows for a different
    size in the y direction.
    """
    g = gauss_kernel(n, sizey=ny)
    improc = scipy.signal.convolve(im,g, mode='valid')
    return(improc)


def cos_filter(y,Fs,hp_st=None,hp_wid=None,lp_st=None,lp_wid=None):
    """
    Does a half-cosine filter of the input data (typically an FFT).  Can 
    do high pass, low pass, or both

    Parameters:
    -----------
    y:  timestream to be filtered

    Fs:  Sample rate of the data to be filtered

    hp_str [None]:  The edge (in the scale of 0 to Fs/2) for the high-pass
        portion of the filter.  This is the lowest x value at which the filter
        is 'off' (i.e. 100% transmissive).

    hp_wid [None]:  The width of the high-pass filter turn on (in units of the
        0-Fs/2 scale).  A wider turn-in leads to less ringing.  The filter goes
        from 0 at hp_str-hp_wid to 1 at hp_str.


    lp_str [None]:  The edge (in the scale of 0 to Fs/2) for the low-pass
        portion of the filter.  This is the highest x value at which the filter
        is 'off' (i.e. 100% transmissive).

    lp_wid [None]:  The width of the low-pass filter turn off (in units of the
        Fs/2 scale).  A wider turn-in leads to less ringing.  The filter goes
        from 1 at lp_str to to 0 at lp_str+lp_wid.

    Returns:
    --------
    y_filt:     filtered timestream
    """

    N = len(y) - len(y)%2

    y = y[:N]

    L, Fny = N/2, Fs/2.
    pad_to = 2**(np.ceil(np.log2(N)))
    npad = (pad_to - N)/2
    pad_y = np.ascontiguousarray(np.zeros(pad_to))
    pad_y[npad:(npad+N)] = y
    pad_y[:(npad+1)]=np.mean(y[:10])
    pad_y[(-npad-1):] = np.mean(y[-10:])
    hann_filt = np.hanning(pad_to)
    pad_y *= hann_filt

    yfft = np.fft.fft(pad_y)
    freq = np.linspace(0,1,pad_to/2+1)*Fny
    filt = np.ones(len(freq))
    if (hp_st is not None and hp_wid is not None):
        x1 = hp_st-hp_wid
        if x1 <0:
            print("You chose a high-pass start and width "
                  "too close to zero frequency.")
            x1=0

        x2 = hp_st
        xinds = np.where((freq>=x1)&(freq<=x2))[0]
        filt[freq<x1]=0
        filt[xinds] = np.flipud(0.5*cos((xinds-xinds[0])*(pi) /
                                        (xinds[-1]-xinds[0])) + 0.5)

    if (lp_st is not None and lp_wid is not None):
        x1 = lp_st
        x2 = lp_st+lp_wid
        if x2 >= Fny:
            print("You chose a low-pass start and width "
                  "too close to nyquist frequency.")
            x2 = Fny
        xinds = np.where((freq>=x1)&(freq<=x2))[0]
        filt[freq>x2]=0
        filt[xinds] = 0.5*cos((xinds-xinds[0])*(pi)/(xinds[-1]-xinds[0]))+0.5
        
    filt = np.append(filt,filt[-1:1:-1])
    y_filt_fft = filt * yfft
    y_fit = np.real(np.fft.ifft(y_filt_fft))
    y_out = y_fit[npad:(npad+N)]/hann_filt[npad:(npad+N)]

    return y_out


def smooth(x, window_len=10, window='hanning'):
    """
    Smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    Parameters:
    -----------
    x: The input signal .
    window_len: The dimension of the smoothing window.
    window: The type of window. Must be one of
        'flat', 'hanning', 'hamming', 'bartlett', and 'blackman'.
        The 'flat' window will produce a moving average smoothing.

    Returns:
    --------
    The smoothed signal
        
    Example:

    >>> import numpy as np    
    >>> t = np.linspace(-2,2,0.1)
    >>> x = np.sin(t)+np.random.randn(len(t))*0.1
    >>> y = smooth(x)
    
    See also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead
    of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'uniform', 'hanning', 'hamming',
                      'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    # Create the window function.
    if window == 'flat' or window=='uniform': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)

    # Apply the window function to the data.
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


def hanning( dim1, dim2=1, alpha=0.5 ):
    """
    Returns a numpy array object which contains a Hanning function. 1D or 2D.

    Parameters:
    -----------
    dim1: Size of the first dimension (ie, length of the array for 1D windows).
         Can also be a subscriptable object with the length of both the
         first and second dimensions.

    dim2 [1]: Size of the second dimension. Ignored if dim1 is a subscriptable
        object.

    alpha [0.5]: Width parameter of generalized Hamming window.
                Alpha=0.5 (default) is the "Hanning" window;
                Alpha=0.54 is the "Hamming" window.

    Returns:
    --------
    An array, of the specified shape, with a Hanning window.
    """

    if not 0.5 <= alpha <= 1.0:
        raise ValueError("Alpha must be in [0.5, 1.0]!")

    # If the first argument is actually a list-like object,
    # then unpack it into dim1 and dim2.
    try:
        dim2=dim1[1]
        dim1=dim1[0]
    except TypeError: pass

    window = (alpha-1.) * np.cos(np.arange(dim1)*2*np.pi/dim1) + alpha
    if dim2 > 1:
        row = np.asmatrix((alpha-1.) * np.cos(np.arange(dim2)*2*np.pi/dim2) +
                          alpha)
        window = np.ascontiguousarray(np.asmatrix(window).T * row)

    return window
