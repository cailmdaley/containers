import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
from spt3g.util import tctools

obsid = '7016644'
File1 = '/spt/data/bolodata/downsampled/calibrator/'+obsid+'/0000.g3'
file2 = '/spt/user/production/calibration/calibrator/'+obsid+'.g3'

f2 = core.G3File(file2)
crframe = f2.next()
names = []
for name in crframe['CalibratorResponseSN'].keys():
    if crframe['CalibratorResponseSN'][name] > 50.:
        names.append(name)

f1 = core.G3File(file1)
tframe = f1.next()
wframe = f1.next()
dframe = f1.next()

psd_dict = {}
npts_psd = 7630/2
for name in names:
    try:
        thisbdata = dframe['RawTimestreams_I'][name]
        thisbdata = thisbdata[2000:12000]
        thisbdata -= np.mean(thisbdata)
        psdtemp = tctools.quick_pspec(thisbdata,npts_psd=npts_psd,rate=76.3)
        psd_dict[name] = psdtemp['psd']
    except:
        pass

psd_tot = np.zeros(npts_psd)
freqs = psdtemp['freqs']
for key in psd_dict.keys():
    psd_tot += psd_dict[key]**2

psd_tot = np.sqrt(psd_tot)
