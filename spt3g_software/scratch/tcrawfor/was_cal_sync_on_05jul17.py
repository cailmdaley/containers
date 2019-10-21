import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
from spt3g.util import tctools
import glob

dir1 = '/spt/data/bolodata/downsampled/calibrator/'
cfiles = glob.glob(dir1+'*/0000.g3')

iobsids = []
for cfile in cfiles:             
    str1=filter(None,cfile.split('/'))
    str2=filter(None,str1[5].split('.'))
    iobsids.append(np.int(str2[0]))
iobsids = np.asarray(iobsids)
aso = np.argsort(iobsids)
siobsids = iobsids[aso]
scfiles = (np.asarray(cfiles))[aso]

# !!!
scfiles = scfiles[0:100]
# !!!

was_on_at_6hz = []

for cfile2 in scfiles:
    f1 = core.G3File(cfile2)       
    try:
        oframe = f1.next()
        wframe = f1.next()
        frame = f1.next()
        sync = np.asarray(frame['CalibratorOn']) - 0.5
        sync = sync[0:256]
        npts = len(sync)
        freqs = np.arange(npts/2)/np.float(npts/2)*76.3/2.
        asyncf = np.abs(np.fft.fft(sync))
        maxfreq = freqs[np.argmax(asyncf[0:len(freqs)])]
        if np.abs(maxfreq-6.) < 0.1:
            was_on_at_6hz.append(1)
        else:
            was_on_at_6hz.append(0)
    except:        
        was_on_at_6hz.append(0)

##file1 = '/spt/data/bolodata/downsampled/calibrator/2815699/0000.g3'
##file1 = '/spt/data/bolodata/downsampled/calibrator/15918320/0000.g3'
#file1 = '/spt/data/bolodata/downsampled/calibrator/2976719/0000.g3'
#f1 = core.G3File(file1)       
#oframe = f1.next()
#wframe = f1.next()
#frame = f1.next()
#bdict = frame['RawTimestreams_I']
#names = bdict.keys()
#names = names[0:100]

