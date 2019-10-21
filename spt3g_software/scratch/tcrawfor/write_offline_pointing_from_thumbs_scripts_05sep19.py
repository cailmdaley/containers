from spt3g import core, std_processing
import numpy as np
import glob, os

# calls should look like this: 
#python offline_pointing_from_field_obs.py /spt/user/axf295/Condor_Return/Calibration_Sources/ra0hdec-44.75/82199407_calibration_thumbnails.g3

data_dir = '/spt/user/axf295/Condor_Return/Calibration_Sources/'
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
tfiles = glob.glob(data_dir+'*/*_calibration_thumbnails.g3')

# loop over sources, find obsids for actual source scans, then find associated elnod and calibrator files
scrfile = 'offline_pointing_from_thumbs_05sep19.scr'
f1 = open(scrfile,'w')
for tfile in tfiles:
    thisstr = 'python offline_pointing_from_field_obs.py ' + tfile
    f1.write(thisstr + '\n')

f1.close()


