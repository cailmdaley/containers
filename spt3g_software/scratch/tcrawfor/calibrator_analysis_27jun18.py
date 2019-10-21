import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
from spt3g.util import tctools

#obsid = '47635164'
obsid = '49022539'

#file1 = '/spt/data/bolodata/downsampled/calibrator/2815699/0000.g3'
#file1 = '/spt/data/bolodata/downsampled/calibrator/15918320/0000.g3'
#file1 = '/spt/data/bolodata/downsampled/calibrator/46795228/0000.g3'
#file1 = '/spt/data/bolodata/downsampled/calibrator/36047770/0000.g3'
file1 = '/spt/data/bolodata/downsampled/calibrator/'+obsid+'/0000.g3'
f1 = core.G3File(file1)       
oframe = f1.next()
wframe = f1.next()
frame = f1.next()
bdict = frame['RawTimestreams_I']
names = bdict.keys()

#calfreq = 6.
calfreq = 4.
samplerate = 76.3
npts = len(frame['CalibratorOn'])

nchunks = 5
chunklen = 2*np.int(np.floor(npts/nchunks)/2.)
freqs = np.arange(chunklen/2)/np.float(chunklen/2)*samplerate/2.
dcf = np.abs(freqs-calfreq)
cfind = np.argmin(dcf)

cdict = {}
dcdict = {}
sndict = {}

for name in names:
    if np.max(np.abs(bdict[name])) > 0.:
        pchunk = np.zeros(nchunks)
        for i in np.arange(nchunks):
            dtemp = bdict[name][chunklen*i:chunklen*(i+1)]
            dtemp -= np.mean(dtemp)
            dtempf = np.abs(np.fft.fft(dtemp))
            pchunk[i] = dtempf[cfind]
        cdict[name] = np.mean(pchunk)
        dcdict[name] = np.std(pchunk)
        sndict[name] = np.mean(pchunk)/np.std(pchunk)

