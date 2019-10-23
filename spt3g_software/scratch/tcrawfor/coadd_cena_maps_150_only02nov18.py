from spt3g import core, dfmux, mapmaker
from spt3g.util import tctools
import numpy as np
import pickle
import glob
from scipy import ndimage

obsid = 56472737

f1 = core.G3File('/spt/user/production/calibration/CenA-pixelraster/maps/'+str(obsid)+'.g3')
f2 = core.G3File('/spt/user/production/calibration/calframe/CenA-pixelraster/'+str(obsid)+'.g3')
calframe = f2.next()
bp = calframe['BolometerProperties']
csndict = calframe['CalibratorResponseSN']

names = []
for name in csndict.keys():
    if csndict[name] > 20:
        names.append(name)
nbolo = len(names)

reso_arcmin = 0.5

bigmap = np.zeros([360,360])
bigwts = np.zeros([360,360])

counter = 0
maxcount = 10
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if counter >= maxcount:
            continue
        if 'Wunpol' in frame and frame['Id'] != 'bsmap' and frame['Id'] != 'azmap' and frame['Id'] != 'elmap':
            twt = np.asarray(frame['Wunpol'].TT)
            print(frame['Id'])
        if frame['Id'] in names:
            bpf = bp[frame['Id']]
            if bpf.band/10. == 150.:
                print frame['Id']
                print bpf.band
                print bpf.x_offset
                print bpf.y_offset
                print counter
                try:
                    xoff = np.int(np.round(bpf.x_offset/core.G3Units.arcmin/reso_arcmin))
                    yoff = np.int(np.round(bpf.y_offset/core.G3Units.arcmin/reso_arcmin))
                    bigmap += tctools.shift(frame['T'],[-yoff,-xoff])
                    bigwts += tctools.shift(twt,[-yoff,-xoff])
                    counter += 1
                except:
                    pass
