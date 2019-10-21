'''
More advanced methods for calculating E and B from Q and U maps, including
    calculateChi
    smith_pure_ft
among others, and various functions that utilize these estimators.
'''
from spt3g import core
from copy import copy
import numpy as np
import scipy.fftpack

from spt3g.util.maths import ddx, d2dx2, ddy, d2dy2, d2dxdy
from .basicmaputils import (
    make_ellgrid, map_to_ft, av_ps, get_fft_scale_fac,
    qu_to_eb_ft, eb_ft_to_qu, check_arrs_size,
    _get_queb_trans_angle
    )
from spt3g.mapspectra import inpaint_map_laplace

fft2_func = scipy.fftpack.fft2
ifft2_func = scipy.fftpack.ifft2

def chi_scaling(ellgrid):
    '''
    The \chi_B maps have a very blue spectrum, see arXiv:astro-ph/0610059v1
    This is the correction factor applied to the estimated A_lms from the
    chib map
    '''
    return ((ellgrid+2)*(ellgrid+1)*(ellgrid)*(ellgrid-1))**-0.5

def calculateChi(q_map, u_map, res, pixel_radius=2, which_chi='B',
                 apodization_window=None):
    """
    Calculates chi_B or chi_E (The B-mode (or E-mode) map, with each
    mode multiplied by sqrt((l-1)(l)(l+1)(l+2))), given an input PolarizedMap
    object with Q and U maps.

    This function should be treating any apodizations of the Q and U maps
    correctly, but the apodization map should be chosen such that both the
    apodization and its derivative are zero at the boundaries.

    Note that any pixel within pixel_radius pixels of the boundaries is
    definitely nonsense, and it might be necessary to remove a few more
    pixels to get rid of artifacts at the edges.

    We assume a flat sky. Under this assumption, we have
    chi_B = (d_x^2 - d_y^2)U - 2(d_xy^2)Q,
    chi_E = (d_x^2 - d_y^2)Q + 2(d_xy^2)U,
    where d denotes the partial derivative. We use the method of finite
    differences to calculate the partial derivatives. The coefficients for each
    term have been chosen so as to minimize E->B leakage.

    The basis of this method is taken from Smith & Zaldarriaga 2006
    (arXiv:astro-ph/0610059v1), equation A4. The coefficients are given in
    Table II from that paper.
    The correction factor which compensates for the apodization is inspired by
    section II of Zhao & Baskaran 2010 (arXiv:1005.1201v3)

    INPUTS
       q_map: (2D ndarray) Q polarization data.

       u_map: (2D ndarray) U polarization data.

       res: (float) Width of a pixel in one of the above maps.

       pixel_radius [2]: (int) How many pixels out should we look when taking
          the partial derivatives? Valid choices are 1, 2, or 3.
          pixel_radius=1 uses only neighboring pixels, pixel_radius=2 uses
          pixels up to 2 pixels away, etc.
          Using more pixels will reduce E->B leakage. A 2-pixel radius should
          be sufficient for any analysis with a realistic level of the B-mode
          gravitational lensing signal.

       which_chi ['B']: (string) Either 'B' to indicate that we should
          calculate chi_B, or 'E', to indicate that we should calculate chi_E.

       apodization_window [None]: (2D ndarray, same size as q_map and u_map)
          This array multiplies the q_map and u_map. Any effects of the
          weighting should be removed by the corrections factors at the end,
          but the weighting function must be chosen such that both the function
          and its derivative are zero at the borders.

    OUTPUT
       A Map with the calculated chi_B or chi_E. The input "map" is also
       modified to add this Map, as the polarization "chiB" (or "chiE").
       Note that the pixel values on the map's border (to a depth of
       pixel_radius pixels) will be nonsense.

    AUTHOR
        Stephen Hoover, March 2013

    CHANGES
        1 Nov 2013: Rename "calculateChi". If input mask is None, then don't do
            any calculations that involve it. SH
        11 Nov 2013: BUGFIX: Switch minus sign on the Q/U swap when calculating
            chiE. SH
    """

    reso_rad = res/core.G3Units.rad


    if pixel_radius==3:
        # Leave this enabled for testing purposes right now.
        core.log_warn("Pixel radius 3 not yet working. The E and B mode maps will be incorrect!", RuntimeWarning)

    if which_chi!='B' and which_chi!='E':
        raise ValueError('You must pick either chi_B or chi_E! '+
                         'Please set which_chi="B" or which_chi="E".')

    # Pixel weights from the Smith & Zaldarriaga (2006) paper (Table 2).
    if pixel_radius==1:
        w = np.array([1., 1./2, 0, 0, 0, 0])
    elif pixel_radius==2:
        w = np.array([4./3, 2./3, -1./12, -1./24, 0, 0])
    elif pixel_radius==3:
        w = np.array([806./2625, 2243./7875, 1501./10500,
                      -239./63000, 53./2625, 907./15750])
    else:
        raise ValueError(
            'The pixel_radius input must be one of the integers 1, 2, and 3!')
    # Divide out the resolution
    w /= reso_rad**2

    if apodization_window is not None:
        window_normalization = np.sqrt( np.sum(apodization_window**2) /
                                       np.product(apodization_window.shape) )
    else:
        window_normalization = None

    # Define a helper function to simplify the notation.
    # roll( arr, [m, n] ) shifts the x-axis (index 1) of the
    # input array by m pixels and shifts the y-axis (index 0)
    # of the input array by n pixels.
    #   Note that I flip the order of the indices! This is
    # so that I can type in the equation with the x-shift first,
    # but still operate on arrays in which the x-axis is last.
    def roll(array, shift):
        out = array
        if shift[0]:
            out = np.roll( out, shift[0], axis=1 )
        if shift[1]:
            out = np.roll( out, shift[1], axis=0 )
        return out

    # Alias the array part of the Q and U Map objects.
    # Note that in the case of chi_E, the "u" and "q"
    # variables are misnamed! Everything's correct for chi_B,
    # but the naming is wrong for chi_E. This was the easiest
    # way to reuse all of the shared code.
    if which_chi=='B':
        q_noweight = np.nan_to_num(q_map)
        u_noweight = np.nan_to_num(u_map)
        if apodization_window is not None:
            q = q_noweight*apodization_window/window_normalization
            u = u_noweight*apodization_window/window_normalization
        else:
            q, u = q_noweight, u_noweight
    elif which_chi=='E':
        u_noweight = np.nan_to_num(q_map)
        q_noweight = -1*np.nan_to_num(u_map)
        if apodization_window is not None:
            u = u_noweight*apodization_window/window_normalization
            q = q_noweight*apodization_window/window_normalization
        else:
            q, u = q_noweight, u_noweight
    else:
        raise ValueError('You must pick either chi_B or chi_E! '+
                         'Please set which_chi="B" or which_chi="E".')

    # Equation A4. We won't bother to add in any terms with a zero weight.
    chi=( w[0]*( roll(u,[-1,0]) - roll(u,[0,-1])
                 + roll(u,[+1,0]) - roll(u,[0,+1]) )  # Alternate +/- signs, changed from paper.
          -w[1]*( roll(q,[+1,+1]) + roll(q,[-1,-1])
                  - roll(q,[-1,+1]) - roll(q,[+1,-1]) ) )
    if w[2]: # Alternate +/- signs, changed from paper.
        chi += w[2]*( roll(u,[-2,0]) - roll(u,[0,-2])
                      + roll(u,[+2,0]) - roll(u,[0,+2]) )
    if w[3]:
        chi += -w[3]*( roll(q,[+2,+2]) + roll(q,[-2,-2])
                      - roll(q,[-2,+2]) - roll(q,[+2,-2]) )
    if w[4]:
        chi += w[4]*( roll(u,[+2,-1]) + roll(u,[+1,-2]) + roll(u,[-1,+2]) + roll(u,[-2,+1])
                      - roll(u,[+2,+1]) - roll(u,[+1,+2]) - roll(u,[-1,-2]) - roll(u,[-2,-1]) )
    if w[5]:
        chi += w[5]*( roll(q,[+2,+1]) + roll(q,[+2,-1]) + roll(q,[-2,+1]) + roll(q,[-2,-1])
                      - roll(q,[+1,+2]) - roll(q,[-1,-2]) - roll(q,[+1,-2]) - roll(q,[-1,+2]) )


    if apodization_window is not None:
        # We used a non-unity weight. Compensate for that by subtracting off
        # the counterterms.
        # Alias for convenience. Divide off the normalization factor:
        win = apodization_window/window_normalization
        res = reso_rad # Alias for convenience
        correction_term = ( u_noweight*(d2dx2(win,res,range=pixel_radius) - d2dy2(win,res,range=pixel_radius))
                            + 2*ddx(win,res,range=pixel_radius)*(ddx(u_noweight,res,range=pixel_radius) - ddy(q_noweight,res,range=pixel_radius))
                            - 2*ddy(win,res,range=pixel_radius)*(ddy(u_noweight,res,range=pixel_radius) + ddx(q_noweight,res,range=pixel_radius))
                            - 2*d2dxdy(win,res,range=pixel_radius)*q_noweight )
        chi -= correction_term
        chi[win != 0] /= win[win != 0]

    ## Finally, zero the pixels on the borders. Those are nonsense.
#     threshold_weight = 0.005
#     thresholded_weight = np.zeros_like(win)
#     thresholded_weight[win>threshold_weight] = 1
#     chi *= thresholded_weight
#     border_window = ndimage.minimum_filter(
#         thresholded_weight, footprint=np.ones(
#             [2*pixel_radius+1,2*pixel_radius+1]))
#     chi *= border_window

    out_chi = copy(q_map)
    np.asarray(out_chi)[:] = chi

    return out_chi

def smith_pure_ft(q_map, u_map, apodization, res = None,
                  do_e=False, pixel_range=2):
    """
    Given an input Q and U map, this function produces the FT of the
    corresponding B-mode (or E-mode) map, modified to remove all ambiguous
    modes. The output FT will be free of leakage from the E modes (or B modes).

    This function is an implementation of Kendrick Smith's pure B-mode
    estimator, first presented in Phys Rev. D 74, 083002 (2006);
    arXiv:astro-ph/0511629v2 . See also http://arxiv.org/abs/astro-ph/0608662.
    The implementation in this function is translated from the IDL function
    spt_analysis/simulations/cl_flatsky_kendrick.pro (also in
    $SPTPOL_SOFTWARE/idl_functions/cl_flatsky_kendrick.pro) by Tom Crawford.

    Note that we define k-space angles differently here than in
    cl_flatsky_kendrick, so the form of the counterterms is slightly changed.
    We here define k-space angles starting from zero in the +ell_y direction
    and increasing toward +ell_x. In cl_flatsky_kendrick, angles are zero in
    the +ell_x direction and increase toward +ell_y.

    INPUTS
        q_map: (2D array or Map) The Stokes Q map.

        u_map: (2D array or Map) The Stokes U map.

        apodization: (2D array or dictionary) The apodization mask. The mask
            must be such that the mask and all derivatives go to zero at its
            boundaries.
            If you input a 2D array, we will take the derivatives numerically
            (using the "pixel_range" input to determine how many surrounding
            pixels we'll use). If you prefer, you may input a dictionary which
            contains the windown and some or all of its pre-computed
            derivatives. We'll look for the following keys:
            'W','dW/dx','dW/dy','d2W/dx2','d2W/dy2','d2W/dx/dy'.
            Any derivatives not present will be computed numerically and stored
            in the dictionary.

        pixel_range [2]: (int) What pixel radius should we use for any numeric
            derivatives? There should be no reason to alter this from the
            defaults.

    OUTPUT
        The k-space pure B modes (or E modes, if do_e is True) corresponding
        to the input Q and U maps.

    EXECPTIONS
        ValueError if we don't get a reso_rad input and neither q_map nor u_map
        are Map objects.
        ValueError if the "apodization" input is a dictionary, but doesn't
        contain a 'W' (for "window") key.

    AUTHOR
        Stephen Hoover, 8 November 2013 (translated from the IDL
        cl_flatsky_kendrick.pro by Tom Crawford)

    CHANGES
        -Ported to spt3g (past Nick)
    """
    # If the "q_map" and "u_map" are Map objects, pull out the map data arrays.
    if res is None:
        res = q_map.res
    reso_rad = res / core.G3Units.rad

    # If we want to calculate pure E modes, swap Q and U.
    # The math is identical, with this substitution.
    if do_e: q_map, u_map = -u_map, q_map

    # Compute the apodization window derivatives if we don't have them already.
    if not isinstance(apodization, dict):
        apodization = {'W':apodization}
    apodization.setdefault('dW/dx',
                           ddx(apodization['W'], reso_rad, pixel_range) )
    apodization.setdefault('dW/dy',
                           ddy(apodization['W'], reso_rad, pixel_range) )
    apodization.setdefault('d2W/dx2',
                           d2dx2(apodization['W'], reso_rad, pixel_range) )
    apodization.setdefault('d2W/dy2',
                           d2dy2(apodization['W'], reso_rad, pixel_range) )
    apodization.setdefault('d2W/dx/dy',
                           d2dxdy(apodization['W'], reso_rad, pixel_range) )

    # Compute k-space angles
    qu_to_eb_angle = _get_queb_trans_angle(q_map, reso_rad)
    #get our ell grid

    ellgrid = make_ellgrid(res, q_map.shape)
    ellgrid[0,0] = 1e6 # Make safe to invert by removing central zero.

    # For convenience and readability, pre-compute the two factors of ell in
    # our final equation. The factor_1 will be attached to the first
    # derivatives of the apodization window, and factor_2 goes will multiply
    # the second derivatives of the apodization window.
    ellgrid_factor_1 = 1. / np.sqrt( (ellgrid-1)*(ellgrid+2) )
    ellgrid_factor_2 = 1. / np.sqrt( (ellgrid-1)*ellgrid*
                                    (ellgrid+1)*(ellgrid+2) )

    # Calcluate various combinations of FTs.
    scale_fac = 1./get_fft_scale_fac(reso_rad, q_map.shape[0],
                                     q_map.shape[1], apodization['W'] )

    map_fft = {}
    map_fft['Q*W'] = fft2_func(q_map*apodization['W']) * scale_fac
    map_fft['U*W'] = fft2_func(u_map*apodization['W']) * scale_fac
    map_fft['Q*dW/dx'] = fft2_func(q_map*apodization['dW/dx']) * scale_fac
    map_fft['Q*dW/dy'] = fft2_func(q_map*apodization['dW/dy']) * scale_fac
    map_fft['U*dW/dx'] = fft2_func(u_map*apodization['dW/dx']) * scale_fac
    map_fft['U*dW/dy'] = fft2_func(u_map*apodization['dW/dy']) * scale_fac
    map_fft['2Q*d2W/dx/dy + U*(d2W/dy2 - d2W/dx2)'] = fft2_func(2.*q_map*apodization['d2W/dx/dy'] +
                                                                  u_map*(apodization['d2W/dy2']
                                                                         -apodization['d2W/dx2'])) * scale_fac
    # Add all of that with the correct prefactors to get pure B estimates.
    # (This will actually be E if we'd set do_e=True up above.)
    map_fft['B'] = ( -np.sin(2*qu_to_eb_angle)*map_fft['Q*W'] + np.cos(2*qu_to_eb_angle)*map_fft['U*W']
                     - 2j*ellgrid_factor_1 * (np.cos(qu_to_eb_angle)*map_fft['Q*dW/dx'] - np.sin(qu_to_eb_angle)*map_fft['Q*dW/dy'] +
                                              np.cos(qu_to_eb_angle)*map_fft['U*dW/dy'] + np.sin(qu_to_eb_angle)*map_fft['U*dW/dx'])
                     - ellgrid_factor_2 * map_fft['2Q*d2W/dx/dy + U*(d2W/dy2 - d2W/dx2)'] 
                 )
    # Finally, normalize :
    #map_fft['B'] /= (np.sqrt(np.sum(apodization['W']**2)/np.product(q_map.shape)))
    #import pdb; pdb.set_trace()
    return map_fft['B']

def get_chi_b_spectra(ell_bins, q, u, apod_mask, inpaint_mask=None,
                      n_iters = 5000):
    ellgrid = make_ellgrid(apod_mask.res, apod_mask.shape)
    chi_b = calculateChi(q, u, q.res, pixel_radius = 2, which_chi = 'B')
    if not inpaint_mask is None:
        inpaint_map_laplace(inpaint_mask, n_iters, chi_b)
    b_ft = map_to_ft(chi_b, apod_mask) * chi_scaling(ellgrid)
    chi_b_cls = av_ps( np.abs(b_ft)**2, apod_mask.res, ell_bins)
    return chi_b_cls

def get_pure_spectra(ell_bins, q, u, apod_mask):
    ellgrid = make_ellgrid(apod_mask.res, apod_mask.shape)
    b_ft = smith_pure_ft(q, u, apod_mask)
    pure_b_cls = av_ps( np.abs(b_ft)**2, apod_mask.res, ell_bins)
    return pure_b_cls

def qu_to_e_only_qu(q,u, apod_mask, kspace_filt = None):
    q = copy(q)
    u = copy(u)
    res = q.res
    check_arrs_size(q,u,apod_mask)
    e_ft, b_ft = qu_to_eb_ft(q,u, apod_mask, kspace_filt = kspace_filt)
    oq, ou = eb_ft_to_qu(e_ft, 0, res, apod_mask)
    np.asarray(q)[:] = oq
    np.asarray(u)[:] = ou
    return q,u

def _apod_inverse(m, apod, cutoff = 0.05):
    apod = copy(apod)
    bad_inds = np.where(apod < cutoff)
    np.asarray(apod)[bad_inds] = 1e6
    inverted = m/apod
    np.asarray(inverted)[bad_inds] = 0
    return inverted

def fancy_qu_to_e_only_qu(q,u, apod_mask, ell_cutoff = 4000, 
                          apod_inv_co=0.05, use_pure = False):
    '''
    The use_pure option uses the pure estimator to get the e mode map and
    it seems to suck, but it's what I would naively expect to work.  It's
    included as a warning.
    That and if you don't read the doc strings you deserve what you get.
    '''
    q = copy(q)
    u = copy(u)
    res = q.res
    check_arrs_size(q,u,apod_mask)

    e_ft, b_ft = qu_to_eb_ft(q,u, apod_mask)

    ell_grid = make_ellgrid(q.res, q.shape)
    
    ell_co = 1.0 - 1.0 / ( 1.0 + (ell_grid / ell_cutoff)**6.0 )

    if use_pure:
        e_pure = smith_pure_ft(q, u, apod_mask, do_e=True)
        oq, ou = eb_ft_to_qu(e_pure, ell_co * b_ft, res, apod_mask)
    else:
        oq, ou = eb_ft_to_qu(e_ft, ell_co * b_ft, res, apod_mask)
        oq = _apod_inverse(oq, apod_mask, apod_inv_co)
        ou = _apod_inverse(ou, apod_mask, apod_inv_co)
    np.asarray(q)[:] = oq
    np.asarray(u)[:] = ou

    return q,u, np.sign(np.abs(oq))
