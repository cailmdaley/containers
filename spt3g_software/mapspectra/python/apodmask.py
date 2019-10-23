from spt3g import core, coordinateutils, mapmaker, std_processing, sources
from spt3g.mapspectra import erode, erode_cardinal
from spt3g.util import files
from copy import copy

import numpy as np
import scipy.signal
import scipy.ndimage

def _get_center_of_mass(m):
    '''
    returns the center of mass in the form [y, x]
    '''
    msum = np.sum(m)
    row_inds, col_inds = np.indices(np.shape(m))
    return ( (1.0/msum * np.sum(row_inds * m)), 
             (1.0/msum * np.sum(col_inds * m)) )
      
def _get_square_bounds(m, top_ind, bottom_ind, horiz_center):
    l_cut = np.asarray(m)[bottom_ind:top_ind, :horiz_center]
    l_bound = np.max(np.where( l_cut == 0)[1])
    r_cut = np.asarray(m)[bottom_ind:top_ind, horiz_center:]
    r_bound = horiz_center + np.min(np.where( r_cut == 0)[1])
    return l_bound, r_bound

def _put_a_square_into_something( m, get_bounds = False ):
    y_cm, x_cm = _get_center_of_mass(m)
    y_cm = int(y_cm)
    x_cm = int(x_cm)

    top_bound = y_cm + np.min(np.where(np.asarray(m)[y_cm:, x_cm] == 0)[0])
    bottom_bound = np.max(np.where(np.asarray(m)[:y_cm, x_cm] == 0)[0])

    best_area = 0
    best_top = None
    best_bottom = None

    for tb in range(top_bound + 1, y_cm, -1):
        local_area = 0
        for bb in range(bottom_bound, y_cm, 1):
            l_bound, r_bound = _get_square_bounds(m, tb, bb, x_cm)
            area = (tb - bb - 1) * (r_bound - l_bound - 1)
            if area > best_area:
                best_area = area
                best_top = tb
                best_bottom = bb
            if area > local_area:
                local_area = area
        if local_area < best_area:
            break
    l_bound, r_bound = _get_square_bounds(m, best_top, best_bottom, x_cm)
    out_m = copy(m)
    np.asarray(out_m)[:] = 0
    np.asarray(out_m)[best_bottom:best_top, l_bound:r_bound] = 1
    if get_bounds:
        return out_m, best_bottom, best_top, l_bound, r_bound
    else:
        return out_m

def _threshold_map(w_tt, cutoff = 0.1, mode='percentile'):
    w = copy(w_tt)
    if mode.lower().startswith('per'):
        cutoff = np.percentile(np.asarray(w_tt)[np.where(w_tt > 0)],
                               100.0 * cutoff)
    elif mode.lower().startswith( 'med'):
        cutoff = cutoff * np.median(np.asarray(w)[np.where(np.asarray(w)!=0.)])
    elif mode.lower().startswith('max'):
        cutoff = cutoff * np.asarray(w).max()
    else:
        raise ValueError('Mode %s not recognized'%mode)
    np.asarray(w)[ np.where(w_tt > cutoff)] = 1    
    np.asarray(w)[ np.where(w_tt <= cutoff)] = 0
    return w

def _smooth_via_erosion(m, n_iters = 10): # XXX: hard-coded 10?
    return erode(m, n_iters, 6) # XXX: why is this 6?

def _smooth_via_transform_and_square(m):
    w_rebin = copy(m)
    w_rebin.proj = coordinateutils.MapProjection.ProjCAR
    w_rebin.res *= 2
    w_rebin.x_res *= 2
    coordinateutils.reproj_map(m, w_rebin)
    w_out = _put_a_square_into_something(w_rebin)
    coordinateutils.reproj_map(w_out, m)
    np.asarray(m)[ np.where(m > 0.99) ] =1 
    np.asarray(m)[ np.where(m < 1) ] = 0
    return m

def _smooth_via_convolution(m, convolution_size = 15):
    kernel = np.asarray(np.matrix(np.hanning(convolution_size)).transpose() * 
                        np.matrix(np.hanning(convolution_size)))
    m = copy(m)
    smoothed = scipy.signal.convolve2d(np.asarray(m), kernel, mode='same')
    np.asarray(m)[np.where(smoothed > 0.5)] = 1
    np.asarray(m)[np.where(smoothed <= 0.5)] = 0
    return m

def generate_apodization_mask(weight_tt, edge_kernel, 
                              bad_threshold_fraction = 0.1, 
                              smooth_kernel_size = 20,
                              pre_mask_erosion = 0,
                              post_mask_erosion = 0,
                              use_square_smoothing = False,
                              skip_enforcement_of_bad_values = False,
                              point_source_mask = None,
                              return_debug_info = False
                          ):
    '''
    DEPRECATED -- use make_border_apodization instead

    Generates an apodization mask, hopefully, a reasonable one.  
    
    It consists of four steps:
     1) Determine which pixels are good by cutting the lowest weight values,
         set by bad_threshold_fraction
     2) It smooths these values out to generate an outline of
         good pixel values.
     3) It reduces the size of this good pixel mask by the
         size of the edge_kernel
     4) It smooths the good pixel mask with the edge kernel provided to
         generate an apodization mask.
    
    It's important to note that this procedure may require some hand holding
    and tweaking of values.  For a given map you should look at the
    apodization mask and check if anything crazy is going on.

    The term erosion in the context of this function refers to eating away
    edge pixels. Pixels are removed from the mask via erosion if their
    distance from an edge pixel is the amount specified, where distance in
    this case means the number of discrete steps along the cardinal directions
    or diagonals taken to get to another pixel.

    After generating the mask, if you are using it for polarized maps you are
    encouraged to double check that the mask is reasonable by plotting it and
    with the function: validate_apodization_mask

    Returns:
      the apodization mask

    Arguments:
      weight_tt:  The TT portion of the weight map.

      edge_kernel:  The 2d array you are smoothing the edges of the map with

      bad_threshold_fraction:  The fraction of pixels cut out of the map

      smooth_kernel_size:  The default smoothing of the good pixel values is
          done via real space convolution this is the size of that smoothing
          kernel.

      pre_mask_erosion:  Amount of erosion inbetween steps 1 and 2.

      post_mask_erosion:  Amount of erosion inbetween steps 3 and 4. 

      use_square_smoothing:  If true, uses a slow, silly method of generating
          the apodization mask.  This method involves transforming the mask
          into a 2d space where ra and dec align along the axes and then
          putting a square into the map that maximizes area while not including
          any bad pixels.  It then transforms back to the original projection.
          This is slow and probably cuts more area than you want, but
          damn are the apodization masks pretty that come out.

      skip_enforcement_of_bad_values:  After step 2 we double check that the 
          pixels included are not pixels we have determined to be bad. If one
          of the pixels in the center of the map is determined to be bad this
          process will result in that pixel being cut and a big old hole will
          wind up in your apodization mask. This usually happens with single
          observation maps.  You can not do this process by setting this
          argument to true.  Just understand that this can result in wonky
          maps/power spectrum with poorly conditioned pixels in the center.
          If you shoot your foot, don't blame me for supplying the gun.

      return_debug_info: returns some extra debug info from the function
    '''
    core.log_warn(
        'generate_apodization_mask is deprecated. Use make_border_apodization')
    if not point_source_mask is None and use_square_smoothing:
        raise RuntimeError(
            "Square smoothing cannot work with a point source mask")

    # sorry this function is a bit piecemeal to make testing each individual
    # section easy. erode makes a pixel 0 if any neighboring pixels are 0
    
    # assign 1s for good data and 0s for bad data.
    m_thresh = _threshold_map(weight_tt, cutoff = bad_threshold_fraction)

    if pre_mask_erosion > 0:
        m_thresh = erode(m_thresh, pre_mask_erosion)

    #smooth this threshold map to produces a map with smooth edges, 
    #still maintains binary values  there are other smoothing functions
    #you can plug in here if your heart desires that
    if use_square_smoothing:
        m_smooth = _smooth_via_transform_and_square(m_thresh)
    else:
        m_smooth = _smooth_via_convolution(
            _smooth_via_erosion(m_thresh),
            convolution_size = smooth_kernel_size)
    #this forces any pixel we want to be set as bad to bad
    if not skip_enforcement_of_bad_values:
        m_smooth = m_smooth * m_thresh

    kernel_shape = np.shape(edge_kernel)
    assert(len(kernel_shape) == 2 and kernel_shape[0] == kernel_shape[1])

    #point source mask code
    if not point_source_mask is None:
        assert(np.max(point_source_mask) == 1)
        assert(np.min(point_source_mask) == 0)
        m_smooth = m_smooth - point_source_mask

    #erode the map by the edge kernel amount.
    erosion_amount = kernel_shape[0]//2 + kernel_shape[0]%2
    erosion_amount = erosion_amount + post_mask_erosion
    m_erode = m_smooth
    for i in range(erosion_amount):
        if i % 2 == 1:
            m_erode = erode_cardinal(m_erode, 1)
        else:
            m_erode = erode(m_erode, 1)

    #convoluve the map with the smoothing kernel,
    # making certain to use realspace convolutions
    apod_mask = scipy.signal.convolve2d(np.asarray(m_erode),
                                        edge_kernel, mode='same')
    np.asarray(m_thresh)[:,:] = apod_mask

    m_thresh /= np.max(m_thresh)

    if return_debug_info:
        return m_thresh, m_erode, m_smooth
    else:
        return m_thresh


def validate_apodization_mask(apod_mask, w_pol, cond_thresh = 1e-2):
    '''
    Checks that an apodization mask does not include poorly conditioned pixels.
    
    This uses the statistic:
        determinant_map / (w_pol.TT * w_pol.QQ * w_pol.UU)
    to determine whether or not a pixel is poorly conditioned because in the
    limit of a very well determined pixel this goes to 1.0

    Arguments:
        apod_mask: the apodization mask
        w_pol: the polarized sky map weight
        cond_thresh: the condition threshold below which the map is said to
            contain poorly conditioned pixels.
    '''
    dmap = np.asarray(mapmaker.mapmakerutils.make_determinant_map(w_pol))
    asign = np.sign(apod_mask)
    if np.any(asign < 0):
        return False, 'Apod mask has negative values, which is weird'
    #in the limit of great pixel converage this should converge to 1
    #because the off diagonal elements of the weight matrix converge to 0
    det_conditioning = np.nan_to_num(dmap / (w_pol.TT * w_pol.QQ * w_pol.UU))
    
    min_cond = np.min(np.abs(det_conditioning) + (1 - asign) )
    if (min_cond < cond_thresh):
        return (False,"You have some poorly conditioned pixels still included"+
                "%f"%min_cond)
    return (True, "A ok, I think...")


@core.usefulfunc
def make_border_apodization(input,
                            res=None,
                            apod_type='cos',
                            custom_kernel=None,
                            threshold_type='median',
                            weight_threshold=0.3,
                            radius_arcmin=20.,
                            zero_border_arcmin=0.,
                            smooth_weights_arcmin=10.,
                            verbose=False,
                            use_square_smoothing=False,
                            apod_threshold=1e-8):
    """
    Takes an input weight map and makes an apodization mask based off of it.
    Adapted from sptpol_software.

    Parameters
    ----------
    input : object
        Either a frame containing Wpol or Wunpol G3SkyMapWeights,
        the SkyMapWeights themselves, the TT portion of weights,
        or a 2d numpy array.
    res : float
        The resolution of the input weight map in G3 Units.
        Only need to specify if `input` is a numpy array.
    apod_type : str
        Supported options include 'cosine', 'gaussian', 'tophat' and 'custom'.
    custom_kernel : ndarray
        If `apod_type` is set to 'custom', set this to be the square array
        with which you are smoothing the edges of the map.
    threshold_type : str
        Threshold the weight map based on 'percentile', 'max', or 'median'.
    weight_threshold : float
        Set the mask to 0 when the weight is below this percentile or this
        fraction of the median or maximum weight, set by `threshold_type`
    radius_arcmin : float
        In arcminutes. Use this as the width of the function with which
        we're convolving the mask, if any.
    zero_border_arcmin : float
        Pad the border of the mask with zeros in to this radius, in arcmin
    smooth_weights_arcmin : float
        In arcminutes. Smooth the input weights with a Gaussian having this
        sigma before thresholding.
    verbose : bool
        Extra screen output.
    use_square_smoothing : bool
        If True, uses a slow, silly method of generating
        the apodization mask.  This method involves transforming the mask
        into a 2d space where ra and dec align along the axes and then
        putting a square into the map that maximizes area while not
        including any bad pixels.  It then transforms back to the original
        projection. This is slow and probably cuts more area than you want,
        but damn are the apodization masks pretty that come out.
    apod_threshold : float
        Make the apodization mask below this threshold 0.0, and make pixels
        in the mask this close to 1.0 exactly equal to 1.0.
        We don't want small mask values later
        being multiplied by crazy-high unweighted map pixels, which started
        happening with rounding issues when we switched to fftconvolve.
        Experimentally for one particular mask, the difference between fft
        and convolution mask is less than 2e-15 everywhere, but the the
        lowest non-zero number was 1e-6, so 1e-8 is conservative.

    Returns
    -------
    If input is an 2d array, returns apodization mask as a 2d array.
    Otherwise returns apodization mask as a FlatSkyMap.

    Stolen from sptpol_software. DPD
    """
    # Sort out the inputs
    if isinstance(input, str):
        if input.endswith('.pkl'):
            input = files.load_pickle(input)
        elif input.endswith('.g3') or input.endswith('.gz'):
            for frame in core.G3File(input):
                if frame.type == core.G3FrameType.Map:
                    input = frame
                    break
        else:
            raise TypeError("I don't know how to load %s"%input)

    if isinstance(input, core.G3Frame):
        if 'Wunpol' in input:
            apod_mask = input['Wunpol'].TT.Clone(True)
        elif 'Wpol' in input:
            apod_mask = input['Wpol'].TT.Clone(True)
        else:
            raise KeyError("No weights in frame.")
    elif isinstance(input, coordinateutils.G3SkyMapWeights):
        apod_mask = input.TT.Clone(True)
    elif isinstance(input, coordinateutils.FlatSkyMap):
        apod_mask = input.Clone(True)
    elif isinstance(input, np.ndarray):
        if res is None:
            raise ValueError(
                "If input is a numpy array, you must specify the resolution.")
        if use_square_smoothing:
            core.log_warn(
                "Square smoothing requires a FlatSkyMap. Setting to False.")
            use_square_smoothing = False
        apod_mask = input
    else:
        raise TypeError(
            "Unsupported %s passed to make_border_apodization"%type(input))

    if res is not None:
        reso_arcmin = res/core.G3Units.arcmin
    else:
        reso_arcmin = apod_mask.res/core.G3Units.arcmin

    if apod_type.lower() == 'custom' and custom_kernel is None:
        raise ValueError(
            "If using a custom apod mask, you must specify the custom_kernel.")

    mask = np.asarray(apod_mask)

    # Smooth weights, if desired.
    if smooth_weights_arcmin > 0.:
        mask = scipy.ndimage.gaussian_filter(
            mask, smooth_weights_arcmin / reso_arcmin / (2*np.sqrt(2*np.log(2))))
        # Not certain we need the normalization.

    # Create the mask, and start by masking pixels with too low a weight.
    mask = _threshold_map(mask, cutoff=weight_threshold, mode=threshold_type)

    if use_square_smoothing:
        np.asarray(apod_mask)[:] = mask
        mask = _smooth_via_transform_and_square(apod_mask)
        mask = np.asarray(mask)

    # If we're doing a cosine apodization, adjust the border zero padding
    # so that the smearing doesn't extend the mask into the region
    # which is supposed to be masked out.
    if apod_type.lower().startswith('cos') and (radius_arcmin > 0):
        if verbose:
            print("Increasing zero_border_arcmin by %f" % (radius_arcmin/2.) +
                  " so that the cosine smearing doesn't put us outside the"
                  " boundaries.")
        zero_border_arcmin += radius_arcmin/2.

    # Add extra zeros around the border.
    if zero_border_arcmin > 0:
        if verbose:
            print("Padding the border with zeros to a depth of "
                  "%f arcminutes (%d pixels)."
                  %(zero_border_arcmin, int(zero_border_arcmin / reso_arcmin)))
        distance_transform = scipy.ndimage.distance_transform_edt(
            mask, sampling=reso_arcmin)
        mask[distance_transform <= zero_border_arcmin] = 0

    # Smear the mask.
    if apod_type.lower().startswith('gaus'):
        # Smear by a Gaussian.
        if verbose:
            print("Smearing the mask with a Gaussian having sigma = "
                  "%f arcminutes." % radius_arcmin)
        mask = scipy.ndimage.gaussian_filter(
            mask, radius_arcmin / reso_arcmin / (2*np.sqrt(2*np.log(2))))
        # Not certain we need the normalization.

    elif apod_type.lower().startswith('cos') and (radius_arcmin > 0):
        # Use a cosine function
        npix_cos = int(radius_arcmin / reso_arcmin)

        if verbose:
            print("Applying cosine apodization with width "
            "%f arcminutes (%d pixels). (2D convolution version)"
            % (radius_arcmin, npix_cos))
        # First generate a convolution kernel based on the requested
        # radius_arcmin
        # Start by calculating the distance of each point from the center,
        # then apply the cosine function.
        xs, ys = np.meshgrid(
            np.linspace(-1., 1., npix_cos+1), np.linspace(-1., 1., npix_cos+1))
        kern_cos = np.sqrt(xs**2 + ys**2)
        outside_circle = kern_cos > 1 # Set these to zero in a moment.
        kern_cos = 0.5 - 0.5 * np.cos((1. - kern_cos) * np.pi)
        kern_cos[outside_circle] = 0.

        # Smear the mask by the just-constructed cosine function.
        # But see note about apod_threshold above.
        mask = scipy.signal.fftconvolve(mask, kern_cos, 'same')

        # Normalize the mask.
        mask /= mask.max()

    elif apod_type.lower().startswith('fastcos') and (radius_arcmin > 0):
        # Use a cosine function - this way is faster than the convolution,
        # but less smooth.
        if verbose:
            print("Applying cosine apodization with width "
                  "%f arcminutes. (Quick and dirty version.)" % radius_arcmin)
        distance_transform = scipy.ndimage.distance_transform_edt(
            mask, sampling=reso_arcmin)
        wh_apply_mask = (distance_transform <= radius_arcmin)
        mask[wh_apply_mask] = 0.5 - 0.5 * np.cos(
            distance_transform[wh_apply_mask] / radius_arcmin * np.pi)

    elif apod_type.lower().startswith('sin') and (radius_arcmin > 0):
        # 1 / (sin(theta) + epsilon) : See eq 36 of astro-ph/0511629v2
        raise NotImplementedError()

    elif apod_type.lower() == 'custom':
        # Enforce a square kernel
        assert(len(np.shape(custom_kernel)) == 2 and
               kernel_shape[0] == kernel_shape[1])
        mask = scipy.signal.convolve2d(mask, custom_kernel, mode='same')
        mask /= mask.max()

    else:
        if verbose: print("Tophat mask (no smearing).")

    if apod_threshold != 0.0:
        mask[np.abs(mask) < apod_threshold] = 0.0
        mask[np.abs(mask - 1.0) < apod_threshold] = 1.0
        # now the fftconvolve mask is identical to the
        # old scipy.ndimage.convolve mask everywhere except the transition,
        # where it's only off by a part in 10^15

    np.asarray(apod_mask)[:] = mask
    return apod_mask

@core.usefulfunc
def make_apodized_ptsrc_mask(input,
                             point_source_file,
                             apod_type='cosine',
                             radius_arcmin=10.,
                             zero_border_arcmin=0,
                             apod_threshold=1e-8,
                             verbose=False):
    """
    Make a smoothed point source mask.

    Parameters:
    -----------
    input : FlatSkyMap or G3Frame containing a FlatSkyMap
        Map parameters for the point source mask.
    point_source_file : str
        The location of a point source configuration file. All sources in
        the file are used, and all use the mask radius specified in the file.
    apod_type : str
        Supported options include 'cosine', 'gaussian', and 'tophat'.
    radius_arcmin : float
        The width of the convolving function.
    zero_border_arcmin : float:
        Pad point sources with zeros out to this radius.
    apod_threshold : (float)
        As with make_border_apodization, make pixels in the mask below this
        threshold equal to 0.0, and make pixels in the mask within this
        threshold from 1.0 exactly equal to 1.0.

    Returns:
    --------
    Smoothed point source mask as FlatSkyMap.

    Stolen from sptpol_software. DPD
    """
    # Sort out input
    if isinstance(input, core.G3Frame):
        if 'T' in input:
            ptsrc_mask = input['T'].Clone(False)
        elif 'Wunpol' in input:
            ptsrc_mask = input['Wunpol'].TT.Clone(False)
        elif 'Wpol' in input:
            ptsrc_mask = input['Wpol'].TT.Clone(False)
        else:
            raise KeyError("Didn't find a FlatSkyMap. Is input a Map frame?")
    elif isinstance(input, coordinateutils.FlatSkyMap):
        ptsrc_mask = input.Clone(False)
    elif isinstance(input, str) and '.g3' in input:
        for frame in core.G3File(input):
            if frame.type == core.G3FrameType.Map:
                ptsrc_mask = frame['T'].Clone(False)
                break
    else:
        raise TypeError(
            "Unsupported %s passed to make_apodized_ptsrc_mask." % type(input))

    reso_arcmin = ptsrc_mask.res/core.G3Units.arcmin
    sources.source_utils.make_point_source_map(ptsrc_mask, point_source_file)

    # We want a mask, not a map:
    mask = np.abs(np.asarray(ptsrc_mask) - 1)

    # This is the same thing that is used in make_border_apodization.
    # If we're doing a cosine apodization, adjust the border zero padding
    # so that the smearing doesn't extend the mask into the region
    # which is supposed to be masked out.
    if apod_type.lower().startswith('cos') and (radius_arcmin > 0):
        if verbose:
            print("Increasing zero_border_arcmin by %f" % (radius_arcmin/2.) +
                  " so that the cosine smearing doesn't put us outside the"
                  " boundaries.")
        zero_border_arcmin += radius_arcmin/2.

    # Add extra zeros around the border.
    if zero_border_arcmin > 0:
        if verbose:
            print("Padding the border with zeros to a depth of "
                  "%f arcminutes (%d pixels)."
                  %(zero_border_arcmin, int(zero_border_arcmin / reso_arcmin)))
        distance_transform = scipy.ndimage.distance_transform_edt(
            mask, sampling=reso_arcmin)
        mask[distance_transform <= zero_border_arcmin] = 0

    # Smear the mask.
    if apod_type.lower().startswith('gaus'):
        # Smear by a Gaussian.
        if verbose:
            print("Smearing the mask with a Gaussian having sigma = "
                  "%f arcminutes." % radius_arcmin)
        mask = scipy.ndimage.gaussian_filter(
            mask, radius_arcmin / reso_arcmin / (2*np.sqrt(2*np.log(2))))
        # Not certain we need the normalization.

    elif apod_type.lower().startswith('cos') and (radius_arcmin > 0):
        # Use a cosine function
        npix_cos = int(radius_arcmin / reso_arcmin)

        if verbose:
            print("Applying cosine apodization with width "
                  "%f arcminutes (%d pixels). (2D convolution version)"
                  % (radius_arcmin, npix_cos))
        # First generate a convolution kernel based on requested radius_arcmin.
        # Start by calculating the distance of each point from the center, then
        # apply the cosine function.
        xs, ys = np.meshgrid(
            np.linspace(-1., 1., npix_cos+1), np.linspace(-1., 1., npix_cos+1))
        kern_cos = np.sqrt(xs**2 + ys**2)
        outside_circle = kern_cos > 1 # Set these to zero in a moment.
        kern_cos = 0.5 - 0.5 * np.cos((1.0 - kern_cos) * np.pi)
        kern_cos[outside_circle] = 0.0

        # Smear the mask by the just-constructed cosine function.
        # mask = ndimage.convolve(mask, kern_cos)
        mask = scipy.signal.fftconvolve(mask, kern_cos, 'same')

        if apod_threshold != 0.0:
            mask[np.abs(mask) < apod_threshold] = 0.0
            mask[np.abs(mask-1.0) < apod_threshold] = 1.0
        # Normalize the mask.
        mask /= mask.max()

    elif apod_type.lower().startswith('fastcos') and (radius_arcmin > 0):
        # Use a cosine function - this way is faster than the convolution,
        # but less smooth.
        if verbose:
            print("Applying cosine apodization with width "
                  "%f arcminutes. (Quick and dirty version.)" % radius_arcmin)
        distance_transform = scipy.ndimage.distance_transform_edt(
            mask, sampling=reso_arcmin)
        wh_apply_mask = (distance_transform <= radius_arcmin)
        mask[wh_apply_mask] = 0.5 - 0.5 * np.cos(
            distance_transform[wh_apply_mask] / radius_arcmin * np.pi)

    elif apod_type.lower().startswith('sin') and (radius_arcmin > 0):
        # 1 / (sin(theta) + epsilon) : See eq 36 of astro-ph/0511629v2
        raise NotImplementedError()

    else:
        if verbose: print("Tophat mask (no smearing).")

    np.asarray(ptsrc_mask)[:] = mask
    return ptsrc_mask
