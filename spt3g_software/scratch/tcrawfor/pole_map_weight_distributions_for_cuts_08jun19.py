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
wrms_dict = {}
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
    wrms_dict[field] = {}
    obsid_dict[field] = {}
    for band in bands:
        weights_dict[field][band] = []
        rms_dict[field][band] = []
        wrms_dict[field][band] = []
        obsid_dict[field][band] = []
    for file1 in file_dict[field]:
#    for file1 in file_dict[field][0:3]:
        print(file1)
        f1 = core.G3File(file1)
        obsid = file1.split('/')[-1].split('_')[0]
        for frame in f1:
            if frame.type is core.G3FrameType.PipelineInfo:
                for key in frame.keys():
                    if hasattr(frame[key],'vcs_githash'):
                        git_hash = frame[key].vcs_githash
                    if key == 'config_file':
                        config_file = frame[key]
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
                        wrms_dict[field][thiskey].append(np.nanstd(mtemp[whgood]*np.sqrt(wtemp[whgood]))**2)
                    else:
                        weights_dict[field][thiskey].append(0.)
                        rms_dict[field][thiskey].append(np.nan)
                        wrms_dict[field][thiskey].append(np.nan)
                    obsid_dict[field][thiskey].append(obsid)
                    
median_dict = {}
std_dict = {}
for field in fields:
    median_dict[field] = {}
    std_dict[field] = {}
    for band in bands:
        median_dict[field][band] = np.median(wrms_dict[field][band])
        std_dict[field][band] = stats.robust_sigma(wrms_dict[field][band])

outdict = {}
outdict['weights'] = weights_dict
outdict['rms'] = rms_dict
outdict['obsid'] = obsid_dict
outdict['wrms'] = wrms_dict
outdict['median'] = median_dict
outdict['std'] = std_dict
outfile = '/spt/user/tcrawfor/weight_rms_stats_08jun19.pkl'
pickle.dump(outdict,open(outfile,'wb'))

# also write text files (will be put in 3g repo)
outfile_indobs = '/home/tcrawfor/code/spt3g_software/std_processing/mapmakers/weight_rms_squared_individual_observation_stats_08jun19.txt'
outfile_dist = '/home/tcrawfor/code/spt3g_software/std_processing/mapmakers/weight_rms_squared_distribution_moments_08jun19.txt'

f = open(outfile_indobs,'w')
f.write('#Written by pole_map_weight_distributions_for_cuts_08jun19.py. \n')
f.write('#Maps created using configuration file \n')
f.write('#   ' + config_file + ',\n')
f.write('# with code at git hash ' + git_hash + '.\n')
f.write('#\n')
f.write('# Field    Band [GHz]   ObsId    Weight*RMS^2 \n')
for field in fields:
    for band in bands:
        for obsid, wrms in zip(obsid_dict[field][band],wrms_dict[field][band]):
            f.write(field+'   '+ '%3s' % band+'   '+obsid+'   '+'{:8.5f}'.format(wrms)+'\n')
f.close()

f = open(outfile_dist,'w')
f.write('#Written by pole_map_weight_distributions_for_cuts_08jun19.py. \n')
f.write('#Maps created using configuration file \n')
f.write('#   ' + config_file + ',\n')
f.write('# with code at git hash ' + git_hash + '.\n')
f.write('#\n')
f.write('# Field    Band [GHz]  Median(Weight*RMS^2)  RobustSigma(Weight*RMS^2) \n')
for field in fields:
    for band in bands:
        f.write(field+'   '+ '%3s' % band+'    '+'{:8.5f}'.format(median_dict[field][band])+'     '+'{:8.5f}'.format(std_dict[field][band])+'\n')
f.close()






