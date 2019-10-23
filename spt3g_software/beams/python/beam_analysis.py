"""
This module contains functions used in analyzing maps to make beam products.
"""
import numpy as _np
import healpy as _hp
import os as _os
import astropy.io.fits as _apf
import scipy.optimize as _opt
import scipy.stats as _stats
import scipy.special as _special
from spt3g import core
from spt3g.util import fitting
from spt3g.util import maths
from spt3g.mapspectra import map_analysis as ma
from spt3g.mapspectra import basicmaputils as bmu
from spt3g.mapmaker import mapmakerutils as mmu

def Bl(ell, fwhm, A):
    """
    ell-space Gaussian beam
    
    Arguments
    ---------
    ell : [array-like]
        Array of multipoles
    fwhm : [float]
        FHWM of beam in degrees
    A : [float]
        Amplitude of beam at ell=0

    Returns
    -------
    b_l : [array-like]
        Gaussian beam array of len(ell)
    """
    sigma = _np.deg2rad(fwhm) / (2. * _np.sqrt(2 * _np.log(2.)))
    return A * _np.exp(-0.5 * (ell * (ell + 1)) * sigma ** 2.)

@core.usefulfunc
def Bl_planet(ell, size=15):
    """
    Window function for correcting non-point-sourceness of planet
    For mars size = approx 15" from 9/28-10/11 (from JPL Horizons)
    Planck13 VII had the formula for window function of a circular disk in ell

    Arguments
    ---------
    ell : [array-like]
        Array of multipoles
    size : [float]
        Diameter of planet in arcsec

    Returns
    -------
    B_disk : [array-like]
        Ell-space window function of the disk
    """
    planet_d = _np.deg2rad(size/3600.)
    B_disk = 2 * _special.jv(1, ell * planet_d) / (ell * planet_d)
    return B_disk

@core.usefulfunc
def fit_gauss_lspace(ell, bl_data, p0=[1/60., 0.02]):
    """
    Fit an ell-space gaussian to the data

    Arguments
    ---------
    ell : [array-like]
        Array of multipoles
    bl_data : [array-like]
        Data vector to be fit
    p0 : [list]
        Initial guesses at [fwhm, amplitude]

    Returns
    -------
    popt : [1-d array]
        Best fit parameters [fwhm, amplitude]
    pcov : [2-d array]
        Covariance of popt
    """
    popt, pcov = _opt.curve_fit(Bl, ell, bl_data, p0)
    return popt, pcov

def get_planck_beam(freq, subset='full', lmax=None):
    """
    Get the Planck B_l for a given frequency and Planck subset
    
    Arguments
    ---------
    freq : [int]
        Planck frequency. Options: 100, 143, or 217
    subset : [str]
        Planck data split. Options: full, hm1, hm2
    lmax : [int]
        Max ell to return.

    Returns
    -------
    b_l : [array-like]
        Planck beam window function
    """
    beam_dir = _os.path.join('/spt', 'user', 'agambrel', 'planck_beams')
    if subset == 'full':
        bfile1 = 'Bl_T_R3.01_fullsky_{0}x{0}.fits'.format(freq)
    else:
        bfile1 = 'Bl_T_R3.01_fullsky_{0}{1}x{0}{1}.fits'.format(freq, 
                                                                subset.lower())

    blp = _apf.getdata(_os.path.join(beam_dir, bfile1))
    blp = blp['TEMPERATURE']
    if lmax is not None:
        if len(blp) > lmax + 1:
            blp = blp[:lmax+1] 
        else:
            blp = _np.pad(blp, (0, lmax + 1 - len(blp)), 'constant')
    return blp

def recenter_maps(fr, map_dict, weight_dict, map_ind, min_width_x, min_width_y):
    """
    Fit a 2-d gaussian to a map, then recenter it around the fit center by
    cutting pixels at the edges.
    """

    f = freq_field(fr['Id'])
    Tmap = _np.asarray(
        mmu.remove_weight_t(fr['T'], fr['Wunpol']))

    # just fit to a small section around the center of the map
    x0 = int(Tmap.shape[0]/4)

    mfit = maths.recenter_map(Tmap, Tmap.shape[0]/2, Tmap.shape[1]/2, widthx=x0,
                              widthy=x0)
    p = fitting.fit_gaussian2d(mfit)
    xcen, ycen = p[2], p[1]
    # x and y indices switch between fitting and recentering
    (x, y) = (xcen+(Tmap.shape[0]-x0)/2, ycen+(Tmap.shape[1]-x0)/2)

    map_dict[f][map_ind] = maths.recenter_map(Tmap, x, y)
    weight_dict[f][map_ind] = maths.recenter_map(_np.asarray(fr['Wunpol'].TT), 
                                                 x, y)

    if map_dict[f][map_ind].shape[0] < min_width_x:
        min_width_x = map_dict[f][map_ind].shape[0]
    if map_dict[f][map_ind].shape[1] < min_width_y:
        min_width_y = map_dict[f][map_ind].shape[1]

    return map_dict, weight_dict, min_width_x, min_width_y

def freq_field(fr_id):
    """
    Look for a frequency in the frame Id. Return matching frequency
    """
    freqs = [90, 150, 220]    
    for freq in freqs:
        if str(freq) in fr_id:
            return freq
    return None

def get_field_maps(field_file, tag='', return_mask=False, 
                   return_mask_split_dec=False):
    """
    Read field maps into a data structure keyed by frequency. Optionally
    also return a mask based on the average weights across the frequencies. 
    Additionally, can return a mask for the top and bottom half of the field
    for computing two separate calibration factors.
    """
    freqs = [90, 150, 220]
    starting_maps_field = _np.zeros(len(freqs), dtype='object')
    field_maps = dict(zip(freqs, starting_maps_field.copy()))
    field_weights = dict(zip(freqs, starting_maps_field.copy()))
    fl = core.G3File(field_file)
    for fr in fl:
        if 'T' in fr and fr['Id'] not in ['PointSourceMask', 
                                          'PointSourceMap',
                                          'bsmap']:
            f = freq_field(fr['Id'])
            
            if tag in fr['Id']:
                field_maps[f] = mmu.remove_weight_t(fr['T'], 
                                                    fr['Wpol'])
                field_weights[f] = _np.asarray(fr['Wpol'].TT)
                res_field = fr['T'].res
                w = fr['Wpol'].TT
    if return_mask:
        # Use average weights over three freqs to make mask
        av_weights = _np.mean(_np.asarray(
            [wmap for freqs, wmap in field_weights.items()]), 
                                  axis=0)
        apod_mask_field = ma.apodmask.make_border_apodization(av_weights, 
                                                            res=res_field)
        if not return_mask_split_dec:
            return field_maps, apod_mask_field
        else:
            # get high and low dec masks to get high and low dec cal
            dec_arr = [w.pixel_to_angle(x)[1] for x in 
                       range(len(_np.asarray(w).ravel()))]
            dec_arr = _np.reshape(dec_arr, w.shape) / core.G3Units.degrees
            mask_low = (dec_arr <= -56) * av_weights
            mask_high = (dec_arr > -56) * av_weights
            mask_low = ma.apodmask.make_border_apodization(mask_low, 
                                                         res=res_field)
            mask_high = ma.apodmask.make_border_apodization(mask_high, 
                                                          res=res_field)
            return field_maps, apod_mask_field, mask_low, mask_high
            
    return field_maps

def get_field_xspecs(field_maps_l, field_maps_r, hm1_maps_l, hm1_maps_r, 
                     hm2_maps_l, hm2_maps_r, apod_mask_field, lmin=50, 
                     lmax=20000, dl=50, freqs=[90, 150, 220],
                     out_file=None):
    """
    Create a dictionary of all of the cross spectra of SPT field maps and 
    mock observed Planck maps, including corrections for Planck beams and
    Planck input map pixel window function.
    """
    field_xspecs = {'sxhm1':{}, 'sxhm2': {}, 'sxs': {}, 'hm1xhm2': {},
                    'sxp_ell': []}
    pfreqs = {90: 100, 150: 143, 220: 217}
    ells = _np.arange(lmin, lmax+1, dl)
    bins = bmu.get_reg_spaced_ell_bins(ells)
    ell_bins = [x[0] for x in bins]
    ell_bins.append(bins[-1][-1])

    for i, f in enumerate(freqs):
        planck_freq = pfreqs[f]
        pw_planck = _hp.pixwin(2048)[:lmax+1]
        if len(pw_planck) < lmax + 1:
            pw_planck = _np.pad(pw_planck, (0, lmax + 1 - len(pw_planck)), 
                                'constant')
        blpT1 = get_planck_beam(planck_freq, subset='hm1', lmax=lmax)
        blpT1 *= pw_planck
        blpT2 = get_planck_beam(planck_freq, subset='hm2', lmax=lmax)
        blpT2 *= pw_planck

        blp_binned1 = _stats.binned_statistic(_np.arange(len(blpT1)), blpT1, 
                                             bins=ell_bins)[0]
        blp_binned2 = _stats.binned_statistic(_np.arange(len(blpT2)), blpT2, 
                                             bins=ell_bins)[0]

        sxp_ell, cls_sxhm1_0 = bmu.simple_cls(
            field_maps_l[f], hm1_maps_r[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        _, cls_sxhm1_1 = bmu.simple_cls(
            field_maps_r[f], hm1_maps_l[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        field_xspecs['sxhm1'][f] = _np.mean([cls_sxhm1_0, cls_sxhm1_1], 
                                           axis=0) / blp_binned1
        field_xspecs['sxp_ell'] = sxp_ell

        _, cls_sxhm2_0 = bmu.simple_cls(
            field_maps_l[f], hm2_maps_r[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        _, cls_sxhm2_1 = bmu.simple_cls(
            field_maps_r[f], hm2_maps_l[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        field_xspecs['sxhm2'][f] = _np.mean([cls_sxhm2_0, cls_sxhm2_1], 
                                           axis=0) / blp_binned2

        _, cls_hm1xhm2_0 = bmu.simple_cls(
            hm1_maps_l[f], hm2_maps_r[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        _, cls_hm1xhm2_1 = bmu.simple_cls(
            hm1_maps_r[f], hm2_maps_l[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        field_xspecs['hm1xhm2'][f] = _np.mean([cls_hm1xhm2_0, cls_hm1xhm2_1],
                                             axis=0) / (blp_binned1 * 
                                                        blp_binned2)

        _, field_xspecs['sxs'][f] = bmu.simple_cls(
            field_maps_l[f], field_maps_r[f], apod_mask=apod_mask_field, 
            ell_min=lmin, ell_max=lmax, delta_ell=dl)
        if out_file is not None:
            _np.save(out_file, field_xspecs)
    return field_xspecs

