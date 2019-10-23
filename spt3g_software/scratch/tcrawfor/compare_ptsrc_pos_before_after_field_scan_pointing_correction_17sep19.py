import numpy as np
from spt3g import core, std_processing, pointing, mapmaker, coordinateutils
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.pointing.check_astrometry_at20 import check_astrometry_at20
import sys, glob
from scipy import ndimage
from astropy.time import Time, TimeDelta

allfiles = []
offset_dict = {}

dir1 = '/spt/user/axf295/Condor_Return/Calibration_Sources/Offline_Pointing/'
dir2 = '/spt/user/axf295/Condor_Return/Calibration_Sources/FieldScan_Pointing/'

files1 = glob.glob(dir1+'*/*.g3*')
files2 = glob.glob(dir2+'*/*.g3*')

for thisfile in files1:
    allfiles.append(thisfile)
for thisfile in files2:
    allfiles.append(thisfile)
allfiles.sort()

for thisfile in allfiles:
    print(thisfile)
    maplist = []
    f1 = core.G3File(thisfile)
    for thisframe in f1:
        if thisframe.type is core.G3FrameType.Map:
            if 'GHz' in thisframe['Id']:
                frame = thisframe
                mm.RemoveWeightModule(frame)
                maplist.append(frame)
#    pdict = check_astrometry_at20(maplist,check_beam_size=False,close_cut_arcsec=360.,nsigma=20.,plot=True,gauss_fit_switch=True)
    pdict = check_astrometry_at20(maplist,check_beam_size=False,close_cut_arcsec=360.,nsigma=10.,plot=True,gauss_fit_switch=True)
    obsid = thisfile.split('/')[-1].split('_')[0]
    if obsid not in offset_dict.keys():
        offset_dict[obsid] = {}
    if 'Offline' in thisfile:
        offset_dict[obsid]['dra_orig'] = pdict['dra']
        offset_dict[obsid]['ddec_orig'] = pdict['ddec']
    else:
        offset_dict[obsid]['dra_corr'] = pdict['dra']
        offset_dict[obsid]['ddec_corr'] = pdict['ddec']
