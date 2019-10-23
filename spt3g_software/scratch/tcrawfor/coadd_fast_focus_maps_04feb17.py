from spt3g import core, dfmux, mapmaker
from spt3g.scratch.tcrawfor import tctools
import numpy as np
import pickle
import glob

elnod_dir = pickle.load(open('elnod_dir_20170203_183206.pkl'))

names = []
for name in elnod_dir.keys():
#    if elnod_dir[name]['elnod_sn'] > 800:
    if elnod_dir[name]['elnod_sn'] > 1200:
        names.append(name)
nbolo = len(names)

#file2 = '/spt_data/bolodata/fullrate/RCW38/2653414/nominal_online_cal.g3'
#f2 = core.G3File(file2)
#bp = f2.next()['NominalBolometerProperties']
#psfac = 41.8/15.
#xoffs = np.asarray([bp[name].x_offset*psfac/core.G3Units.arcmin for name in names])
#yoffs = np.asarray([bp[name].y_offset*psfac/core.G3Units.arcmin for name in names])
reso_arcmin = 0.5
#xoffs_pix = xoffs/reso_arcmin
#yoffs_pix = yoffs/reso_arcmin
f1 = core.G3File('/home/nwhitehorn/saturnbpm2.g3')
bpframe = f1.next()
bp = bpframe['BolometerProperties']

sfiles = glob.glob('/poleanalysis/cal_results/saturn/*.g3')
sfiles.sort()
sfiles = sfiles[8:15]
nfiles = len(sfiles)

focmaps = np.zeros([nfiles,360,360])

for i in np.arange(nfiles):
    maptemp = np.zeros([360,360])
    wtemp = np.zeros([360,360])
    f1 = core.G3File(sfiles[i])
    for frame in f1:
        if frame.type is core.G3FrameType.Map:
            if 'Wunpol' in frame:
                if np.max(np.asarray(frame['Wunpol'].TT)) > 0.:
                    twt = np.asarray(frame['Wunpol'].TT)
            if frame['Id'] in names:
                try:
                    xoff = np.int(np.round(bp[frame['Id']].x_offset/core.G3Units.arcmin/reso_arcmin))
                    yoff = np.int(np.round(bp[frame['Id']].y_offset/core.G3Units.arcmin/reso_arcmin))
                    map_weighted = np.asarray(frame['T'])
                    maptemp += tctools.shift(map_weighted,[-yoff,-xoff])
                    wtemp += tctools.shift(twt,[-yoff,-xoff])
                except:
                    pass
    wh0 = np.where(wtemp == 0.)
    if np.size(wh0) > 1:
        wtemp[wh0] = 1.
    maptemp /= wtemp
    focmaps[i,:,:] = maptemp
