import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
from spt3g.scratch.tcrawfor import tctools

def grab_data(frame, cmap_dict, data1 = [], data2 = [], bdict = {}):

    if frame.type != core.G3FrameType.Timepoint:
        return
#    try:
    data1.append(frame['CalibratorOn'])
    data2.append(frame['EventHeader'])
    for key in cmap_dict:
        inds = cmap_dict[key]
        bdict[key].append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
#    except:
#        pass

file1 = '/buffer/bolodata/20170213_200044.g3'

calfreq = 9.

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
bnames = np.asarray(wmap.keys())
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

# !!! this gets all bolos
bolos2get = bnames
#bolos2get = bnames[67:72]

cmap_dict = {}
for key in bolos2get:
    wmk = wmap[key]
    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
    cmap_dict[key] = ilist
names = cmap_dict.keys()

nbolo = len(names)

cdict = {}
dcdict = {}
for name in names:
    dcdict[name] = 1e12

nchunks = 5
chunklen = 1526
freqs = np.arange(1526/2)/10.
cfind = np.int(calfreq/10.)

data1 = []
data2 = []
bdict = {}
for name in names:
    bdict[name] = []
    
for frame in core.G3File(file1):
    grab_data(frame, cmap_dict, data1 = data1, data2 = data2, bdict=bdict)
        
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

