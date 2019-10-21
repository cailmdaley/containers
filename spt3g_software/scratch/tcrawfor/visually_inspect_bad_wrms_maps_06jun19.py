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

#outfile = '/spt/user/tcrawfor/weight_rms_stats_06jun19.pkl'
outfile = '/spt/user/tcrawfor/weight_rms_stats_08jun19.pkl'
d1 = pickle.load(open(outfile,'rb'))
wrms_dict = d1['wrms']
obsid_dict = d1['obsid']
median_dict = d1['median']
std_dict = d1['std']
badobsdict = {}
nddict = {}

badthresh = 5.
nbad = 0

for field in median_dict.keys():
    badobsdict[field] = {}
    nddict[field] = {}
    for band in median_dict[field].keys():
        badobsdict[field][band] = []
        nddict[field][band] = {}
        thismed = median_dict[field][band]
        thisstd = std_dict[field][band]
        for wrms, obsid in zip(wrms_dict[field][band],obsid_dict[field][band]):
            normdev = np.abs(wrms - thismed)/thisstd
            nddict[field][band][obsid] = normdev
            if normdev > badthresh:
                nbad += 1
                badobsdict[field][band].append(obsid)

badwts = np.zeros([nbad,1000,1000])
badmaps = np.zeros([nbad,1000,1000])
badfield = []
badband = []
badobs = []

qmap = 0
for field in median_dict.keys():
    for band in median_dict[field].keys():
        for obsid in badobsdict[field][band]:
            badfield.append(field)
            badband.append(band)
            badobs.append(obsid)
            file1 = '/spt/data/onlinemaps/'+field+'/'+obsid+'_'+band+'GHz_tonly.g3.gz'
            print()
            print(file1)
            print('weight-rms^2 product is '+str(nddict[field][band][obsid])+' sigma high.')
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
                        mm.RemoveWeightModule(frame)
                        mtemp = np.asarray(frame['T'])[ymin:ymax,xmin:xmax]
                        badwts[qmap,:,:] = wtemp
                        badmaps[qmap,:,:] = mtemp
                        qmap += 1
                            


