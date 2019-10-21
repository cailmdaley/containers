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

converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')
dir1 = '/spt_data/bolodata/fullrate/calibrator/3265029/'
file2 = dir1 + 'nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']
file3 = dir1 + '0000.g3'
f3 = core.G3File(file3)
oframe = f3.next()
wframe = f3.next()
converter(wframe)
dframe = f3.next()
converter(dframe)

psd_dict = {}
npts_psd = 1024
for name in names:
    try:
        thisbdata = dframe['CalTimestreams'][name]
        thisbdata -= np.mean(thisbdata)
        psdtemp = tctools.quick_pspec(thisbdata,npts_psd=npts_psd)
        psd_dict[name] = psdtemp['psd']
    except:
        pass

start_time = core.G3Time('20170207_184400')
stop_time = core.G3Time('20170207_184700')


files = ['/buffer/bolodata/20170207_184343.g3',
         '/buffer/bolodata/20170207_184433.g3',
         '/buffer/bolodata/20170207_184524.g3',
         '/buffer/bolodata/20170207_184614.g3']

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
#bolos2get = bnames[67:72]
bolos2get = pickle.load(open('elnod_gnames_08feb17.pkl'))

cmap_dict = {}
for key in bolos2get:
    wmk = wmap[key]
    ilist = [wmk.board_serial, wmk.module, wmk.channel*2]
    cmap_dict[key] = ilist
names = cmap_dict.keys()

nbolo = len(names)

data1 = []
data2 = []
bdict = {}
for name in names:
    bdict[name] = []
    
for file1 in files:
    for frame in core.G3File(file1):
        grab_data(frame, cmap_dict, data1 = data1, data2 = data2, bdict=bdict)

#bdict2 = {}
#for name in names:
#    bdict2[name] = tctools.list_to_array(bdict[name])
bdict2 = bdict

mjds = []
for d2 in data2:
    mjds.append(d2.mjd)
mjds = np.asarray(mjds)
whg = np.where(np.logical_and(mjds >= start_time.mjd, mjds < stop_time.mjd))
minwh = np.min(whg[0])
maxwh = np.max(whg[0])

tstot = np.zeros(len(bdict2[name]))
tstot = np.zeros(maxwh - minwh + 1)
for name in names:
    bdt = bdict2[name][minwh:maxwh+1]
    tstot += bdt - np.mean(bdt)

npts = np.int(np.floor(np.float(len(tstot))/2.)*2.)
pstot = np.abs(np.fft.fft(tstot[0:npts]))
pstot = pstot[0:npts/2]
freqs = np.arange(npts/2)/np.float(npts/2)*152.6/2.

boards = np.asarray([(wmap[name]).board_serial for name in names])
uboards = np.unique(boards)
boarddata_dict = {}
nb_dict = {}
for board in uboards:
    boarddata_dict[board] = np.zeros(maxwh - minwh + 1)
    nb_dict[board] = 0.
bpsd_dict = {}
bpt_ampl1 = {}
bpt_ampl2 = {}
for name in names:
    bdt = bdict2[name][minwh:maxwh+1]
    boarddata_dict[(wmap[name]).board_serial] += bdt - np.mean(bdt)
    nb_dict[(wmap[name]).board_serial] += 1.
for board in uboards:
    boarddata_dict[board] /= nb_dict[board]
    bpsd_dict[board] = np.abs(np.fft.fft(boarddata_dict[board][0:npts]))
    bpt_ampl1 = bpsd_dict[board][216]
    bpt_ampl2 = bpsd_dict[board][255]
