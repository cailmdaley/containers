from spt3g import core, dfmux, mapmaker
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle
import glob

elnod_dir = pickle.load(open('elnod_dir_20170203_183206.pkl'))

names = []
for name in elnod_dir.keys():
    if elnod_dir[name]['elnod_sn'] > 800:
#    if elnod_dir[name]['elnod_sn'] > 1200:
        names.append(name)
nbolo = len(names)

reso_arcmin = 0.5

f1 = core.G3File('/home/nwhitehorn/saturnbpm2.g3')
bpframe = f1.next()
bp = bpframe['BolometerProperties']

fname0 = '/poleanalysis/sptdaq/20170129_rawdata/20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
cframe = f1.next()
bp2 = cframe['NominalBolometerProperties']
wafnames = np.asarray([bp2[name].wafer_id for name in names])
uwafers = np.unique(wafnames)
uwafers.sort()
nwafers = len(uwafers)

sfiles = glob.glob('/poleanalysis/cal_results/saturn/*.g3')
sfiles.sort()
sfiles = sfiles[8:15]
nfiles = len(sfiles)

focmap_dir = {}
for wafer in uwafers:
    focmap_dir[wafer] = np.zeros([nfiles,360,360])

for i in np.arange(nfiles):
#for i in np.arange(1) + 4:
    mapdir = {}
    wtsdir = {}
    for wafer in uwafers:
        mapdir[wafer] = np.zeros([360,360])
        wtsdir[wafer] = np.zeros([360,360])
    f1 = core.G3File(sfiles[i])
    for frame in f1:
        if frame.type is core.G3FrameType.Map:
            if 'Wunpol' in frame:
                if np.max(np.asarray(frame['Wunpol'].TT)) > 0.:
                    twt = np.asarray(frame['Wunpol'].TT)
            if frame['Id'] in names:
                if bp2[frame['Id']].band/10. == 150.:
                    try:
                        xoff = np.int(np.round(bp[frame['Id']].x_offset/core.G3Units.arcmin/reso_arcmin))
                        yoff = np.int(np.round(bp[frame['Id']].y_offset/core.G3Units.arcmin/reso_arcmin))
                        thiswid = bp2[frame['Id']].wafer_id
                        map_weighted = np.asarray(frame['T'])
                        mapdir[thiswid] += tctools.shift(map_weighted,[-yoff,-xoff])
                        wtsdir[thiswid] += tctools.shift(twt,[-yoff,-xoff])
                    except KeyError:
                        pass
    for wafer in uwafers:
        wh0 = np.where(wtsdir[wafer] == 0.)
        if np.size(wh0) > 1:
            wtsdir[wafer][wh0] = 1.
        mapdir[wafer] /= wtsdir[wafer]
        focmap_dir[wafer][i,:,:] = mapdir[wafer]
