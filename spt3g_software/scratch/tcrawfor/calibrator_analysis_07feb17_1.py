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
    try:
#        data1.append(frame['CalibratorOn'])
        data2.append(frame['EventHeader'])
        for key in cmap_dict:
            inds = cmap_dict[key]
            bdict[key].append(frame['DfMux'][inds[0]][inds[1]][inds[2]])
    except:
        pass

files = ['/buffer/bolodata/20170206_194933.g3','/buffer/bolodata/20170206_195254.g3']
#files = ['/buffer/bolodata/20170206_194933.g3']

calfreqs = [6.,15.]

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
bnames = np.asarray(wmap.keys())
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

# !!! this gets all bolos
#bolos2get = bnames
bolos2get = bnames[67:72]

cmap_dict = {}
for key in bolos2get:
    wmk = wmap[key]
    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
    cmap_dict[key] = ilist
names = cmap_dict.keys()

nbolo = len(names)

cdict_6hz_6hz = {}
cdict_15hz_6hz = {}
cdict_6hz_15hz = {}
cdict_15hz_15hz = {}
dcdict_6hz_6hz = {}
dcdict_15hz_6hz = {}
dcdict_6hz_15hz = {}
dcdict_15hz_15hz = {}

nchunks = 5
chunklen = 1526

for file1 in files:

    data1 = []
    data2 = []
    bdict = {}
    for name in names:
        bdict[name] = []

    for frame in core.G3File(file1):
        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, bdict=bdict)

    bdict2 = {}
    for name in names:
        bdict2[name] = tctools.list_to_array(bdict[name])
    
