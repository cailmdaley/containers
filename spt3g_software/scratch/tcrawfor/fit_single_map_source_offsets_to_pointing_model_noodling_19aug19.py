import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp, coordinateutils, pointing
import os
import glob
from spt3g.mapmaker import mapmakerutils as mm

# sandbox for developing code to fit source offsets from single-obs
# maps to obtain offline pointing corrections

dir1 = '/spt/data/onlinemaps/'

fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
obsids_all = []
for field1 in fields:
    files = glob.glob(dir1 + field1 + '/*150GHz_tonly.g3.gz')
    for file1 in files:
        obsid1 = file1.split('/')[-1].split('_')[0]
        if np.int(obsid1) > 70000000:
            obsids_all.append(obsid1)
obsids_all.sort()
stride = np.int(np.floor(len(obsids_all)/20))
#stride = np.int(np.floor(len(obsids_all)/3))
obsids2do = obsids_all[::stride]

bands = ['90','150','220']
kernel_dict = {}
kernel_dict['90'] = 1.8
kernel_dict['150'] = 1.5
kernel_dict['220'] = 1.2

offset_dict = {}

for band1 in bands:
    for field1 in fields:
        files = glob.glob(dir1 + field1 + '/*' + band1 + 'GHz_tonly.g3.gz')
        files.sort()

        for file1 in files:

            obsid1 = file1.split('/')[-1].split('_')[0]
            if obsid1 not in obsids2do:
                continue

            print(file1)
            if obsid1 not in offset_dict.keys():
                offset_dict[obsid1] = {}

            f1 = core.G3File(file1)
            for frame in f1:
                if frame.type is core.G3FrameType.Map:
                    if 'GHz' in frame['Id']:
                        break

            fwhm_pix = kernel_dict[band1]/(frame['T'].res/core.G3Units.arcmin)
            map_sm = ndimage.gaussian_filter(np.asarray(frame['T']),fwhm_pix*0.42)
            map_new = coordinateutils.FlatSkyMap(map_sm,res=frame['T'].res,is_weighted=False,alpha_center=frame['T'].alpha_center,delta_center=frame['T'].delta_center,proj=frame['T'].proj,pol_type=core.MapPolType.T)
            frame_new = core.G3Frame(core.G3FrameType.Map)
            frame_new['Wunpol'] = frame['Wunpol']
            frame_new['T'] = map_new

            wtemp = np.asarray(frame_new['Wunpol'].TT)
            wtemp_sm = ndimage.gaussian_filter(wtemp,16.)
            medwt = np.median(wtemp_sm[np.where(wtemp_sm > np.max(wtemp_sm)/2.)])
            whok = np.where(wtemp_sm > 0.8*medwt)
            pixel_mask = np.zeros(wtemp_sm.shape)
            pixel_mask[whok] = 1

            pdict = pointing.check_astrometry_at20.check_astrometry_at20(frame_new,check_beam_size=False,close_cut_arcsec=360.,pixel_mask=pixel_mask,nsigma=20.,plot=False)

            if pdict is not None:
                offset_dict[obsid1][band1] = pdict

