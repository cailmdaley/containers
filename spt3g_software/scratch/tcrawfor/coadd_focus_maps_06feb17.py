from spt3g import core, dfmux, mapmaker
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle
import glob

f4 = core.G3File('/poleanalysis/cal_results/elnod/3091165.g3')
enframe = f4.next()

gnames = []
for name in enframe['ElnodSlopes'].keys():
    ent = enframe['ElnodSigmaSlopes'][name]
    ert = enframe['ElnodSlopes'][name]
    if ent > 0.:
        if -ert/ent > 400.:
            gnames.append(name)
nbolo = len(gnames)

file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp2 = f2.next()['NominalBolometerProperties']
f1 = core.G3File('/home/nwhitehorn/saturnbpm2.g3')
bpframe = f1.next()
bp = bpframe['BolometerProperties']

reso_arcmin = 0.5

#file1 = '/poleanalysis/cal_results/0537-441/3091301.bolomaps.g3.gz'
#files = [file1]
files = ['/poleanalysis/cal_results/0537-441/3091301.bolomaps.g3.gz',
         '/poleanalysis/cal_results/0537-441/3095973.bolomaps.g3.gz',
         '/poleanalysis/cal_results/0537-441/3100645.bolomaps.g3.gz',
         '/poleanalysis/cal_results/0537-441/3105317.bolomaps.g3.gz',
         '/poleanalysis/cal_results/0537-441/3109993.bolomaps.g3.gz']
nfiles = len(files)

print files

focmaps = np.zeros([nfiles,360,360])

for i in np.arange(nfiles):
    maptemp = np.zeros([360,360])
    wtemp = np.zeros([360,360])
    print files[i]
    f1 = core.G3File(files[i])
    for frame in f1:
        if frame.type is core.G3FrameType.Map:
            if 'Wunpol' in frame:
                if np.max(np.asarray(frame['Wunpol'].TT)) > 0.:
                    twt = np.asarray(frame['Wunpol'].TT)
            if frame['Id'] in gnames:
                try:
                    xoff = np.int(np.round(bp[frame['Id']].x_offset/core.G3Units.arcmin/reso_arcmin))
                    yoff = np.int(np.round(bp[frame['Id']].y_offset/core.G3Units.arcmin/reso_arcmin))
                    map_weighted = np.asarray(frame['T'])
                    if np.min(map_weighted) < 0. and np.min(map_weighted) > -1e-12:
                        maptemp += tctools.shift(map_weighted,[-yoff,-xoff])
                        wtemp += tctools.shift(twt,[-yoff,-xoff])
                except:
                    pass
    wh0 = np.where(wtemp == 0.)
    if np.size(wh0) > 1:
        wtemp[wh0] = 1.
    maptemp /= wtemp
    focmaps[i,:,:] = maptemp
