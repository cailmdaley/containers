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

reso_arcmin = 0.5

wafnames = np.asarray([bp[name].wafer_id for name in names])
uwafers = np.unique(wafnames)
uwafers.sort()
nwafers = len(uwafers)

mapdict = {}
wtsdict = {}
for wafer in uwafers:
    mapdict[wafer] = np.zeros([360,360])
    wtsdict[wafer] = np.zeros([360,360])

for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bsmap' and frame['Id'] != 'azmap' and frame['Id'] != 'elmap':
            twt = np.asarray(frame['Wunpol'].TT)
            print(frame['Id'])
        if frame['Id'] in names:
            if bp[frame['Id']].band/10. == 150.:
                try:
                    amap = np.asarray(frame['T'])
                    mapshape = np.shape(amap)
                    amap_sm = ndimage.gaussian_filter(amap,4)
                    ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
                    xoff = xcenter - 180
                    yoff = ycenter - 180
                    thiswid = bp[frame['Id']].wafer_id
                    mapdict[thiswid] += tctools.shift(amap,[-yoff,-xoff])
                    wtsdict[thiswid] += tctools.shift(twt,[-yoff,-xoff])
                except KeyError:
                    pass

for wafer in uwafers:
    wh0 = np.where(wtsdict[wafer] == 0.)
    if np.size(wh0) > 1:
        wtsdict[wafer][wh0] = 1.
        mapdict[wafer] /= wtsdict[wafer]
