import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

# get files
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
file_dict = {}
for field in fields:
    files = glob.glob('/spt/data/bolodata/downsampled/'+field+'/7*/0000.g3')
    files.sort()
    file_dict[field] = files
    
# declination vector
dec = np.arange(3800)/3800.*38.-76.
ddec = dec[1] - dec[0]
dec_mean = -56.
wvec0 = np.cos(dec_mean*np.pi/180.)/np.cos(dec*np.pi/180.)

# define field weights
wvec_dict = {}
nhann = np.int(np.round(2./ddec))
htemp = np.hanning(2*nhann)
hann_left = htemp[0:nhann]
hann_right = htemp[nhann:]
for field in fields:
    wvec_dict[field] = wvec0.copy()
    dec_mid = -np.float(field.split('-')[-1])
    dec_min = dec_mid - 4.75
    dec_max = dec_mid + 4.75
    ind1 = np.max(np.where(dec <= dec_min)[0])
    ind2 = np.min(np.where(dec >= dec_max)[0])
    wvec_dict[field][0:ind1] = 0.
    wvec_dict[field][ind1:ind1+nhann] *= hann_left
    wvec_dict[field][ind2:] = 0.
    wvec_dict[field][ind2-nhann:ind2] *= hann_right

# now make slice by adding one copy of field weights for every obs
slice = np.zeros(len(dec))
for field in fields:
    for file1 in file_dict[field]:
        slice += wvec_dict[field]

