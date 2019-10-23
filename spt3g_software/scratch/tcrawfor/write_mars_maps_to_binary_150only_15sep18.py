from spt3g import core, dfmux, mapmaker
from spt3g.util import tctools
import numpy as np
import pickle
import glob
from scipy import ndimage

obsid = 52279854

f1 = core.G3File('/spt/user/production/calibration/mars-pixelraster/singlebolomaps/'+str(obsid)+'.g3')
f2 = core.G3File('/spt/user/production/calibration/calframe/mars-pixelraster/'+str(obsid)+'.g3')
calframe = f2.next()
bp = calframe['BolometerProperties']
csndict = calframe['CalibratorResponseSN']

names = []
for name in csndict.keys():
    if csndict[name] > 20:
        names.append(name)
nbolo = len(names)

maps = np.zeros([nbolo,360,360])

qqq = 0
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bsmap' and frame['Id'] != 'azmap' and frame['Id'] != 'elmap':
            twt = np.asarray(frame['Wunpol'].TT)
        if frame['Id'] in names:
            if bp[frame['Id']].band/10. == 150.:
                print(frame['Id'])
                try:
                    maps[qqq,:,:] = np.asarray(frame['T'])
                    qqq += 1
                except KeyError:
                    pass

