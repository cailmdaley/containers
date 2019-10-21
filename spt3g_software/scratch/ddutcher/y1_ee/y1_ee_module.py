'''
y1_ee_module.py

Collection of utilities for the Y1 1500d EE/TE analysis
'''
import os
import argparse
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from spt3g import core, calibration, coordinateutils, util, mapspectra
from spt3g.util import files
from spt3g.mapspectra import map_analysis
import y1_ee.average_cls as av
import plotting as pltng
import y1_ee.null_tests as null

spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
map_dir = '/spt/user/ddutcher/'


def make_bundle_masks(bundle_dir='/spt/user/ddutcher/bundles/150GHz/total',
                      output_dir='/spt/user/ddutcher/bundles/apods',
                      verbose=False):
    """
    Read in all the bundles in bundle_dir and make border apodization
    masks for each of them, saving as .pkl files.
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    bundles = sorted(glob(os.path.join(bundle_dir, '*.g3*')))
    if len(bundles)==0:
        raise FileNotFoundError('No g3 files in %s'%bundle_dir)

    for bundle in bundles:
        if verbose:
            print(bundle)
        root = os.path.basename(bundle).partition('.')[0]
        for fr in core.G3File(bundle):
            if fr.type == core.G3FrameType.Map:
                for k in fr.keys():
                    # Allow 'Wpol' and 'Wunpol'
                    if 'W' in k and 'pol' in k:
                        wgt = fr[k]
                        break
            break

        apod = mapspectra.apodmask.make_border_apodization(
            wgt, radius_arcmin=0., zero_border_arcmin=0,
            smooth_weights_arcmin=5)

        files.save_pickle(np.array(apod),
                          os.path.join(output_dir,root+'_apod.pkl'))


def combine_bundle_masks(bundle_dir='/spt/user/ddutcher/bundles/150GHz/total',
                         apod_dir='/spt/user/ddutcher/bundles/apods',
                         output_dir='/spt/user/ddutcher/masks',
                         tag='',
                         make_ptsrc_mask=True,
                         ptsrc_file = os.path.join(
                             spt3g_software,
                             'sources/1500d_ptsrc_150GHz_50mJy.txt'),
                         radius_arcmin=30.,
                         zero_border_arcmin=10.0,
                         weight_threshold=0.3,
                         ptsrc_radius_arcmin=10,
                         save_masks=True,
                         verbose=False,
                         return_data=False):
    """
    This will read in individual apod masks, or make them if they
    don't exist, and find the intersection of them in order to make
    one combined mask.
    """
    # Need to get map shape and resolution. Bundle apods might just be arrays.
    pth = glob(os.path.join(bundle_dir, '*.g3*'))[0]
    bun = core.G3File(pth)
    for map_frame in bun:
        if map_frame.type == core.G3FrameType.Map:
            break

    if not os.path.exists(apod_dir):
        os.mkdir(apod_dir)
    apods = sorted(glob(os.path.join(apod_dir,'*pkl')))
    if len(apods) == 0:
        print('No apod masks in %s'%apod_dir)
        print('Making apods from bundles now')
        make_bundle_masks(bundle_dir, out_dir = apod_dir,
                          verbose=verbose)
    apods = sorted(glob(os.path.join(apod_dir,'*_apod.pkl')))

    for i,pth in enumerate(apods):
        if verbose:
            print('Loading ', pth)
        now_apod = files.load_pickle(pth)
        if i==0:
            total_apod = now_apod
        else:
            total_apod *= now_apod

    total_apod = mapspectra.apodmask.make_border_apodization(
        total_apod, apod_type = 'cos', weight_threshold = weight_threshold,
        radius_arcmin = radius_arcmin, zero_border_arcmin=zero_border_arcmin,
        smooth_weights_arcmin=0, res = map_frame['T'].res)
    if save_masks:
        files.save_pickle(total_apod,
                          os.path.join(output_dir, tag+'bundle_apod_mask.pkl'))
    if make_ptsrc_mask:
        if verbose:
            print('Making ptsrc mask')
        ptsrc_mask = mapspectra.apodmask.make_apodized_ptsrc_mask(
            map_frame, ptsrc_file, apod_type = 'cos',
            radius_arcmin = ptsrc_radius_arcmin, zero_border_arcmin=0)
        total_apod *= ptsrc_mask
        if save_masks:
            files.save_pickle(ptsrc_mask, os.path.join(
                output_dir, tag+'bundle_ptsrc_mask.pkl'))
            files.save_pickle(total_apod, os.path.join(
                output_dir, tag+'bundle_total_mask.pkl'))
    if return_data:
        return total_apod


def analyze_all_nulls(data_dir='/spt/user/ddutcher/null_xspectra/old',
                      bands=['90GHz','150GHz','220GHz'],
                      tests=['1_2','azimuth','left_right','saturation','moon','wafer'],
                      spectra=['EE','TE'],
                      save_plots=True,
                      plot_dir='/home/ddutcher/plots'):
    plt.switch_backend('Agg')
    ptes = []
    for band in bands:
        print(band)
        for test in tests:
            print(test)
            now_dir = os.path.join(data_dir, band, test)
            results = av.average_cls(now_dir, spectra=spectra,rebin=np.arange(275,3025,50))
            for spec in spectra:
                pte = pltng.plot_ps_with_err(
                    results, spectrum=spec, label_pte = True, save_plots=save_plots,
                    title = ' '.join([test, band, spec]), return_pte=True,
                    filename=os.path.join(
                        plot_dir, '_'.join([spec.lower(), band,test, 'null.png'])),
                    linestyle='none', marker='.')
                plt.close('all')
                ptes.append(pte)
    return np.asarray(ptes)


if __name__=="__main__":
#     make_bundle_masks(bundle_dir = '/spt/user/ddutcher/sim_bundles/planck_full/*/high_el',
#                       output_dir = '/spt/user/ddutcher/bundles/tmp_high-el_apods',
#                       verbose = True)

    combine_bundle_masks(bundle_dir = '/spt/user/ddutcher/bundles/90GHz/total',
                         apod_dir = '/spt/user/ddutcher/bundles/tmp_high-el_apods',
                         output_dir = '/spt/user/ddutcher/masks',
                         tag = '3band_res2_high-el',
                         verbose = True)
