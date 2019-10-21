import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob
from spt3g.mapmaker import mapmakerutils as mm
from spt3g.mapmaker import summingmaps as sm
from spt3g.util import stats

outfile = '/spt/user/tcrawfor/weight_rms_stats_08jun19.pkl'
d1 = pickle.load(open(outfile,'rb'))
wrms_dict = d1['wrms']
obsid_dict = d1['obsid']
median_dict = d1['median']
std_dict = d1['std']
badobsdict = {}
nddict = {}

#badthresh1 = 10.
#badthresh2 = 20.
badthresh1 = 3.
badthresh2 = 5.
nbad = 0

coadd_all_dict = {}
coadd_cut1_dict = {}
coadd_cut2_dict = {}
weights_all_dict = {}
weights_cut1_dict = {}
weights_cut2_dict = {}
sign_all_dict = {}
sign_cut1_dict = {}
sign_cut2_dict = {}

for field in median_dict.keys():
    sign_all_dict[field] = {}
    sign_cut1_dict[field] = {}
    sign_cut2_dict[field] = {}
    for band in median_dict[field].keys():
        sign_all_dict[field][band] = 1.
        sign_cut1_dict[field][band] = 1.
        sign_cut2_dict[field][band] = 1.

for field in median_dict.keys():
    coadd_all_dict[field] = {}
    coadd_cut1_dict[field] = {}
    coadd_cut2_dict[field] = {}
    weights_all_dict[field] = {}
    weights_cut1_dict[field] = {}
    weights_cut2_dict[field] = {}
    for band in median_dict[field].keys():
        thismed = median_dict[field][band]
        thisstd = std_dict[field][band]
        coadd_all_dict[field][band] = np.zeros([1000,1000])
        coadd_cut1_dict[field][band] = np.zeros([1000,1000])
        coadd_cut2_dict[field][band] = np.zeros([1000,1000])
        weights_all_dict[field][band] = np.zeros([1000,1000])
        weights_cut1_dict[field][band] = np.zeros([1000,1000])
        weights_cut2_dict[field][band] = np.zeros([1000,1000])
        for wrms, obsid in zip(wrms_dict[field][band],obsid_dict[field][band]):
#        for wrms, obsid in zip(wrms_dict[field][band][0:3],obsid_dict[field][band][0:3]):
            file1 = '/spt/data/onlinemaps/'+field+'/'+obsid+'_'+band+'GHz_tonly.g3.gz'
            print()
            print(file1)
            # what pixels to use for weight and rms
            if field == 'ra0hdec-44.75':
                xmin = 8500
                xmax = 9500
                ymin = 8500
                ymax = 9500
            if field == 'ra0hdec-52.25':
                xmin = 8500
                xmax = 9500
                ymin = 6700
                ymax = 7700
            if field == 'ra0hdec-59.75':
                xmin = 8500
                xmax = 9500
                ymin = 4900
                ymax = 5900
            if field == 'ra0hdec-67.25':
                xmin = 8500
                xmax = 9500
                ymin = 3100
                ymax = 4100
            f1 = core.G3File(file1)
            for frame in f1:
                if frame.type is core.G3FrameType.Map:
                    if 'GHz' in frame['Id']:
                        wtemp = np.asarray(frame['Wunpol'].TT)[ymin:ymax,xmin:xmax]
                        mtemp = np.asarray(frame['T'])[ymin:ymax,xmin:xmax]
            coadd_all_dict[field][band] += mtemp*sign_all_dict[field][band]
            weights_all_dict[field][band] += wtemp
            sign_all_dict[field][band] *= -1.
            normdev = np.abs(wrms - thismed)/thisstd
            print('weight-rms^2 product is '+str(normdev)+' sigma high.')
            if normdev < badthresh1:
                coadd_cut1_dict[field][band] += mtemp*sign_cut1_dict[field][band]
                weights_cut1_dict[field][band] += wtemp
                sign_cut1_dict[field][band] *= -1.
            if normdev < badthresh2:
                coadd_cut2_dict[field][band] += mtemp*sign_cut2_dict[field][band]
                weights_cut2_dict[field][band] += wtemp
                sign_cut2_dict[field][band] *= -1.

for field in weights_all_dict.keys():
    for band in weights_all_dict[field].keys():
        whn0 = np.where(weights_all_dict[field][band] > 0.)
        coadd_all_dict[field][band][whn0] /= weights_all_dict[field][band][whn0]
        whn0 = np.where(weights_cut1_dict[field][band] > 0.)
        coadd_cut1_dict[field][band][whn0] /= weights_cut1_dict[field][band][whn0]
        whn0 = np.where(weights_cut2_dict[field][band] > 0.)
        coadd_cut2_dict[field][band][whn0] /= weights_cut2_dict[field][band][whn0]
        
