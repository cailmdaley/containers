"""
This module contains functions to do cross-spectra
with Planck maps for calibration purposes.
Most of it is cribbed from Anne's code.
"""

import os
import numpy as np
from glob import glob
from astropy.io import fits
from spt3g import core
from spt3g.simulations.instrument import get_beams
from spt3g.util import files
from spt3g.mapspectra import map_analysis

from y1_ee import y1_ee_cross_spectra as cs


def get_planck_beam(freq, subset='full', lmax=None, teb='T'):
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
    teb : [str]
        Beam to return. 'T', 'E', or 'B'.

    Returns
    -------
    b_l : [array-like]
        Planck beam window function
    """
    assert subset in ['full', 'hm1', 'hm2']
    assert teb in ['T', 'E', 'B']
    
    beam_dir = '/spt/user/agambrel/planck_beams'
    if subset == 'full' or int(freq) != 100:
        bfile1 = 'Bl_TEB_R3.00_fullsky_{0}x{0}.fits'.format(freq)
    else:
        # Half-mission TEB beams only available for 100 GHz
        bfile1 = 'Bl_TEB_R3.00_fullsky_{0}{1}x{0}{1}.fits'.format(
            freq, subset.lower())

    blp = fits.getdata(os.path.join(beam_dir, bfile1))
    blp = blp[teb]
    if lmax is not None:
        if len(blp) > lmax + 1:
            blp = blp[:lmax+1] 
        else:
            blp = np.pad(blp, (0, lmax + 1 - len(blp)), 'constant')
    return blp


if __name__ == "__main__":
#     splits=['high-az', 'high-resp', 'right',
#             'low-az', 'low-resp', 'left',
#             'even', 'even-az', 'odd', 'odd-az',
#             'first-half', 'second-half']

    spec, cal = 'EE', 'Pcal'
#     spec, cal = 'TT', 'Tcal'
    
    splits = ['full']

    masks = glob('/spt/user/ddutcher/masks/3band_res2_*bundle_total_mask.pkl')
    
    spt_beam = get_beams(fwhm_90=1.7, fwhm_150=1.4, fwhm_220=1.2, lmax=10000-1)
#     spt_beam = get_beams(beamfile='/home/ddutcher/code/spt3g_software/beams/products/beam00.txt',
#                          lmax=10000-1)

    for band, pfreq in [('220GHz', 217), ('90GHz', 100), ('150GHz', 143) ]:
        if band == '90GHz':
            lo = 0.9331
            hi = 0.9934
        elif band == '150GHz':
            lo = 0.8882
            hi = 0.9648
        elif band == '220GHz':
            lo = 0.9400
            hi = 0.9551
        else:
            raise ValueError(band)
        print('Doing '+band)
        
        template = cs.make_planck_calibration_template(high_dec=hi, low_dec=lo)
        
        redo_spt = True
        if redo_spt:
            left = list(core.G3File(
                '/spt/user/ddutcher/coadds/20190917_left_{}.g3.gz'.format(band)))[-1]
            cs.apply_planck_calibration(left, template)
            map_analysis.deproject_tp(left, int(band.strip('GHz')))

            right = list(core.G3File(
                '/spt/user/ddutcher/coadds/20190917_right_{}.g3.gz'.format(band)))[-1]
            cs.apply_planck_calibration(right, template)
            map_analysis.deproject_tp(right, int(band.strip('GHz')))
        else: 
            spt_x_spt = files.load_pickle(
                '/spt/user/ddutcher/xspectra/{}/left_x_right.pkl'.format(band))

        planck_map = list(core.G3File(
            '/spt/user/ddutcher/coadds/20190916_planck_full_{}.g3.gz'.format(band)))[-1]
        planck_beam = get_planck_beam(pfreq, lmax=10000-1, teb=spec[0])

        for mask in masks:
            apod = files.load_pickle(mask)
            if redo_spt:
                spt_x_spt =  map_analysis.calculate_powerspectra(
                    left, input2=right, apod_mask=apod, b_mode_method='basic', calculate_dls=True)
            
            mask_label = mask.split('bundle')[0].split('_')[-1]
            if mask_label == '':
                mask_label = 'full'
            print('***Using %s field***' % mask_label)
        
#             results = np.zeros((12, 10000))
            for i, split in enumerate(splits):
                print('Doing %s haf-depth %s map' % (split, band))
                pth = '/spt/user/ddutcher/coadds/20190917_{}_{}.g3.gz'.format(split, band)
                if not os.path.exists(pth):
                    raise FileNotFoundError(pth)
                map1 = list(core.G3File(pth))[-1]
                cs.apply_planck_calibration(map1, template)
                map_analysis.deproject_tp(map1, int(band.strip('GHz')))

                planck_x_spt = map_analysis.calculate_powerspectra(
                    map1, input2=planck_map, apod_mask=apod, b_mode_method='basic', calculate_dls=True)
                ratio = (spt_beam[int(band.strip('GHz'))] * planck_x_spt[spec]
                         / (planck_beam * spt_x_spt[spec]))
#                 results[i,:] = ratio

            files.save_pickle(
                ratio,
                '/spt/user/ddutcher/planck/planck_{}_{}_{}_fulldepth_tp_corr.pkl'.format(
                    cal, band, mask_label)
            )