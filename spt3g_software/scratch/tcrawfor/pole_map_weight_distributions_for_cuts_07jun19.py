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

# get files
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
file_dict = {}
for field in fields:
    files = glob.glob('/spt/data/onlinemaps/'+field+'/7*.g3.gz')
    files.sort()
    file_dict[field] = files

# initialize output
weights_dict = {}
rms_dict = {}
obsid_dict = {}
bands = ['90','150','220']
    
# loop over files, get weight and map rms for each
for field in fields:
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
    weights_dict[field] = {}
    rms_dict[field] = {}
    obsid_dict[field] = {}
    for band in bands:
        weights_dict[field][band] = []
        rms_dict[field][band] = []
        obsid_dict[field][band] = []
    for file1 in file_dict[field]:
#    for file1 in file_dict[field][0:3]:
        print(file1)
        f1 = core.G3File(file1)
        obsid = file1.split('/')[-1].split('_')[0]
        for frame in f1:
            if frame.type is core.G3FrameType.Map:
                if 'GHz' in frame['Id']:
                    thiskey = frame['Id'].split('GHz')[0]
                    wtemp = np.asarray(frame['Wunpol'].TT)[ymin:ymax,xmin:xmax]
                    whgood = np.where(wtemp > 0.)
                    if len(whgood[0]) > 100:
                        weights_dict[field][thiskey].append(np.nanmean(wtemp[whgood]))
                        mm.RemoveWeightModule(frame)
                        mtemp = np.asarray(frame['T'])[ymin:ymax,xmin:xmax]
                        rms_dict[field][thiskey].append(np.nanstd(mtemp[whgood]))
                    else:
                        weights_dict[field][thiskey].append(0)
                        rms_dict[field][thiskey].append(np.nan)
                    obsid_dict[field][thiskey].append(obsid)
                    
wrms_dict = {}
median_dict = {}
std_dict = {}
for field in fields:
    wrms_dict[field] = {}
    median_dict[field] = {}
    std_dict[field] = {}
    for band in bands:
        wrms_dict[field][band] = np.asarray(weights_dict[field][band])*(np.asarray(rms_dict[field][band])**2)
        median_dict[field][band] = np.median(wrms_dict[field][band])
        std_dict[field][band] = stats.robust_sigma(wrms_dict[field][band])

outdict = {}
outdict['weights'] = weights_dict
outdict['rms'] = rms_dict
outdict['obsid'] = obsid_dict
outdict['wrms'] = wrms_dict
outdict['median'] = median_dict
outdict['std'] = std_dict
outfile = '/spt/user/tcrawfor/weight_rms_stats_07jun19.pkl'
pickle.dump(outdict,open(outfile,'wb'))





