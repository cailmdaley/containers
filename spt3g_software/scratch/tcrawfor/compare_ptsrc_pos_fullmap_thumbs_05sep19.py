import numpy as np
from spt3g import core, std_processing, pointing, mapmaker, coordinateutils
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.pointing.astrometry import check_astrometry_at20
from spt3g.pointing import offline_pointing as op
from spt3g.util.fitting import fit_gaussian2d
import argparse as ap
import sys, glob, pickle
from scipy import ndimage
from astropy.time import Time, TimeDelta

thumbdir = '/spt/user/axf295/Condor_Return/Calibration_Sources/'
fulldir = '/spt/data/onlinemaps/'
tfiles = glob.glob(thumbdir + '*/*.g3*')
tfiles.sort()
ra_dict = {}
dec_dict = {}
npix_cutout = 60
npix_cutout_half = np.int(npix_cutout/2)
for tfile in tfiles:
    print(tfile)
    str1 = tfile.split('/')
    obsid = str1[-1].split('_')[0]
    field = str1[-2]
    ffile = fulldir + field + '/' + obsid + '_90GHz_tonly.g3.gz'
    for thisframe in core.G3File(ffile):
        if thisframe.type is core.G3FrameType.Map:
            if 'GHz' in thisframe['Id']:
                fullframe = thisframe
                mm.RemoveWeightModule(fullframe)
    fwhm_pix = 0.9/(fullframe['T'].res/core.G3Units.arcmin)
    fullmap_sm = ndimage.gaussian_filter(np.asarray(fullframe['T']),fwhm_pix*0.42)
    ra_dict[obsid] = {}
    dec_dict[obsid] = {}
    for thisframe in core.G3File(tfile):
        if thisframe.type is core.G3FrameType.Map:
            if 'GHz' in thisframe['Id']:
                try:
                    thumbframe = thisframe
                    mm.RemoveWeightModule(thumbframe)
                    fwhm_pix = 0.9/(thumbframe['T'].res/core.G3Units.arcmin)
                    thumbmap_sm = ndimage.gaussian_filter(np.asarray(thumbframe['T']),fwhm_pix*0.42)
                    ra_dict[obsid][thumbframe['Id']] = {}
                    dec_dict[obsid][thumbframe['Id']] = {}
                    tmapshape = thumbframe['T'].shape
                    #                tycenter, txcenter = np.unravel_index(np.argmax(np.abs(thumbmap_sm)),[tmapshape[0],tmapshape[1]])
                    #                tradec_t = thumbframe['T'].pixel_to_angle(np.int(txcenter),np.int(tycenter))
                    params_t = fit_gaussian2d(thumbmap_sm,fit_offset=True)
                    tradec_t = thumbframe['T'].xy_to_angle(params_t[1],params_t[2])
                    ra_dict[obsid][thumbframe['Id']]['thumb'] = tradec_t[0]
                    dec_dict[obsid][thumbframe['Id']]['thumb'] = tradec_t[1]
                    fmapshape = fullframe['T'].shape
                    #                tempcenter = np.unravel_index(fullframe['T'].angle_to_pixel(tradec_t[0],tradec_t[1]),[fmapshape[0],fmapshape[1]])
                    tempcenter = np.unravel_index(fullframe['T'].angle_to_pixel(thumbframe['T'].alpha_center,thumbframe['T'].delta_center),[fmapshape[0],fmapshape[1]])
                    fullmap_sm_cut = fullmap_sm[tempcenter[0]-npix_cutout_half:tempcenter[0]+npix_cutout_half,tempcenter[1]-npix_cutout_half:tempcenter[1]+npix_cutout_half]
                    #                fycenter, fxcenter = np.unravel_index(np.argmax(np.abs(fullmap_sm_cut)),[npix_cutout,npix_cutout])
                    #                tradec_f = fullframe['T'].pixel_to_angle(np.int(tempcenter[1]-npix_cutout_half+fxcenter),np.int(tempcenter[0]-npix_cutout_half+fycenter))
                    params_f = fit_gaussian2d(fullmap_sm_cut,fit_offset=True)
                    tradec_f = fullframe['T'].xy_to_angle(tempcenter[1]-npix_cutout_half+params_f[1],tempcenter[0]-npix_cutout_half+params_f[2])
                    ra_dict[obsid][thumbframe['Id']]['full'] = tradec_f[0]
                    ra_dict[obsid][thumbframe['Id']]['delta'] = tradec_f[0] - tradec_t[0]
                    dec_dict[obsid][thumbframe['Id']]['full'] = tradec_f[1]
                    dec_dict[obsid][thumbframe['Id']]['delta'] = tradec_f[1] - tradec_t[1]
                except:
                    pass

d = {}
d['ra'] = ra_dict
d['dec'] = dec_dict

pickle.dump(d,open(' /spt/user/tcrawfor/compare_ptsrc_pos_fullmap_offline_thumbs_05sep19.pkl'))
