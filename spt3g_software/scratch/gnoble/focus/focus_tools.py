import os
import numpy as np

from spt3g.core import G3File
from spt3g.pointing import focus

PI = np.pi
TRTLT = 2.0*np.sqrt(2.0*np.log(2.0)) # ~ 2.355

def rad2arcmin(rad):
    return rad*60.*180./PI

def gaussian_1d(x, x0=0., var=1., norm=None):
    """Evaluate a one-dimensional Gaussian at position `x`."""
    if norm is None:
        norm = 1./np.sqrt(2.*PI*var)
    return norm*np.exp(-((x - x0)**2.)/(2.*var))

def gaussian_2d(y, x, y0, x0, var_y, var_x, theta, norm=1.):
    """Evaluate a two-dimensional Gaussian at position (y, x)."""    
    a = (np.cos(theta)**2.)/(2.*var_x) + (np.sin(theta)**2.)/(2.*var_y)
    b = (-np.sin(2.*theta))/(4.*var_x) + (np.sin(2.*theta))/(4.*var_y)
    c = (np.sin(theta)**2.)/(2.*var_x) + (np.cos(theta)**2.)/(2.*var_y)    
    return norm*np.exp(-(a*((x-x0)**2.) + 2.*b*(x-x0)*(y-y0) + c*((y-y0)**2.)))

def var2fwhm(var_y, var_x=None):
    """Compute the one or two-dimensional Gaussian FWHM."""
    if var_x is None: # 1-dimensional case
        return TRTLT*np.sqrt(var_y)
    return TRTLT*np.sqrt((var_y + var_x)/2.)

def var2ellip(var_y, var_x):
    """Compute the two-dimensional Gaussian ellipticity."""
    return (var_y - var_x)/(var_y + var_x)

def fwhm2var(fwhm, ellip=None):
    """Compute the Gaussian variances(s).
    
    The larger (semi-major axis) variance will always be returned first.
    """
    if ellip is None: # 1-dimensional case
        return (fwhm/TRTLT)**2.

    var_y = (ellip + 1.)*(fwhm/TRTLT)**2.
    var_x = var_y*(1. - ellip)/(1. + ellip)

    return (var_y, var_x)

def annular_indices(shape, inner, outer, i0=None, j0=None):
    """Retrieve array indices that lie in an annulus, centered on (i0, j0).
    
    Annular indices are strictly greater than the inner radius, and not greater
    than to the outer radius, in this definition.
    
    Args:
        shape: array shape to consider
        inner: inner radius of the annulus
        outer: outer radius of the annulus
        i0, j0: index of the center of the annulus, the array center by default
    Returns:
        ndarray of the annular values
    """
    if i0 is None:
        i0 = (shape[0] - 1)/2.
    if j0 is None:
        j0 = (shape[1] - 1)/2.

    # compute Euclidean distances from (i0, j0)
    i, j = np.indices(shape)
    r = np.sqrt((i - i0)**2. + (j - j0)**2.)
    
    return np.all([inner < r, r <= outer], axis=0).nonzero()

def crop(arr, shape, i0=None, j0=None):
    """Crop a 2D array, centered on a specified point.
    
    Note that, to avoid interpolation, when `arr.shape` and `shape` differ by
    an odd number, the cropped region will be offset from the requested center
    by one index in both axes, towards the origin.
    
    Args:
        arr: 2D array to be cropped
        shape: crop shape
        i0, j0: index of the center of the windowed region (assumed to be the
            center by default)
    Returns:
        slice of `arr` with shape `shape`
    """
    if i0 is None:
        i0 = arr.shape[0]/2.0
    if j0 is None:
        j0 = arr.shape[1]/2.0
    
    i_min, j_min = int(round(i0 - shape[0]/2.0)), int(round(j0 - shape[1]/2.0))
    i_max, j_max = i_min + int(shape[0]), int(j_min + shape[1]) 

    if (i_min < 0) or (j_min < 0) or (i_max >= arr.shape[0]) or (j_max >= arr.shape[1]):
        raise ValueError("Chosen center ({}, {}) is too close to the edge for"
                         "the given shape: {} and array size:"
                         "{}".format(i0, j0, shape, arr.shape))
    
    return arr[i_min:i_max, j_min:j_max]

def get_optics_bench_position(obsid, bolodata_dir, optical_coords=False):
    """Retrieve an observation's commanded optics bench position.

    Only the `0000.g3` is examined, making the implicit assumption that bench
    position won't change.
    
    Args:
        obsid: observation id, a str
        bolodata_dir: directory that bolodata lives in, a str or Path
        optical_coords: return the bench coordinates relative to the optical
            axis, if True.
    Returns:
        ndarray of offsets relative to the optical axis; the default system is
        the typical XYZ coordinates, the `optical_coords` are, as viewed from
        the primary: [left/right (+/-), up/down (+/-), and away/towards (+/-)]
    """
    for frame in G3File(os.path.join(bolodata_dir, str(obsid), '0000.g3')):
        if 'BenchCommandedPosition' in frame.keys():
            position = focus.g3bench2xyz(fr['BenchCommandedPosition']) 
            if optical_coords:
                position = np.array(focus.bench2optical(*bench_pos)).flatten()
            return position
    raise ValueError('no bench position found for {} in {}'.format(obsid, bolodata_dir))
