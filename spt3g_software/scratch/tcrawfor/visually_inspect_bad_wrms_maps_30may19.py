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

outfile = '/spt/user/tcrawfor/weight_rms_stats_30may19.pkl'
d1 = pickle.load(open(outfile,'rb'))
median_dict = d1['median']
std_dict = d1['std']

# get files
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
file_dict = {}
for field in fields:
    files = glob.glob('/spt/data/onlinemaps/'+field+'/7*.g3.gz')
    files.sort()
    file_dict[field] = files
    
# initialize output
wmap_dict = {}
map_dict = {}
bands = ['90','150','220']
badobs = []    

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
    wmap_dict[field] = {}
    map_dict[field] = {}
    for band in bands:
        wmap_dict[field][band] = {}
        map_dict[field][band] = {}
    for file1 in file_dict[field]:
#    for file1 in file_dict[field][0:3]:
        print(file1)
        f1 = core.G3File(file1)
        obsid = file1.split('/')[-1].split('_')[0]
        for frame in f1:
            if frame.type is core.G3FrameType.Map:
                if 'GHz' in frame['Id']:
                    thiskey = frame['Id'].split('GHz')[0]
#                    wmap_dict[field][thiskey][obsid] = np.asarray(frame['Wunpol'].TT)[ymin:ymax,xmin:xmax]
#                    mm.RemoveWeightModule(frame)
#                    map_dict[field][thiskey][obsid] = np.asarray(frame['T'])[ymin:ymax,xmin:xmax]
#                    thiswrms = np.nanmean(wmap_dict[field][thiskey][obsid])*(np.nanstd(map_dict[field][thiskey][obsid]))**2
                    wtemp = np.asarray(frame['Wunpol'].TT)[ymin:ymax,xmin:xmax]
                    mm.RemoveWeightModule(frame)
                    mtemp = np.asarray(frame['T'])[ymin:ymax,xmin:xmax]
                    thiswrms = np.nanmean(wtemp)*(np.nanstd(mtemp))**2
                    print()
                    print(thiswrms)
                    print(median_dict[field][thiskey])
                    print(thiswrms-median_dict[field][thiskey])
                    print(std_dict[field][thiskey])
                    print(np.abs(thiswrms-median_dict[field][thiskey])/std_dict[field][thiskey])
                    print()
                    if np.abs(thiswrms-median_dict[field][thiskey]) > 3.*std_dict[field][thiskey]:
                        badobs.append(obsid)
                    
