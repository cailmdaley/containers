from spt3g import core, dfmux, mapmaker
from spt3g.util import tctools
import numpy as np
import pickle
import scipy
from scipy import ndimage

obsid = '2923042'

file2 = '/spt/data/bolodata/fullrate/saturn/'+obsid+'/nominal_online_cal.g3'
f2 = core.G3File(file2)
bp = f2.next()['NominalBolometerProperties']

mapdict = {}
bands = {}
file1 = '/spt/analysis/production/calibration/saturn/maps/'+obsid+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt[np.where(twt == 0.)] = 1.
        if frame['Id'] != 'bs':
            map_weighted = np.asarray(frame['T'])
            if np.max(np.abs(map_weighted)) > 0.:
                map_unw = map_weighted.copy()
                map_unw /= twt
                mapdict[frame['Id']] = map_unw

gnames = mapdict.keys()
bands = np.asarray([np.int(bp[name].band/10.) for name in gnames])

mapshape = np.shape(mapdict[gnames[0]])
map90 = np.zeros(mapshape)
n90 = 0
map150 = np.zeros(mapshape)
n150 = 0
map220 = np.zeros(mapshape)
n220 = 0

for i in np.arange(len(gnames)):
    name = gnames[i]
    amap = np.zeros(mapshape)
    amap[80:350,80:350] = -(mapdict[name])[80:350,80:350]*1e12
    amap_sm = ndimage.gaussian_filter(amap,4)
    if np.max(amap_sm)/np.std(amap_sm[330:,330:]) > 20.:
        ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
        amap_shift = np.roll(amap,mapshape[0]/2-ycenter,0)
        amap_shift = np.roll(amap_shift,mapshape[1]/2-xcenter,1)
        if bands[i] == 150:
            map150 += amap_shift
            n150 += 1
        if bands[i] == 90:
            map90 += amap_shift
            n90 += 1
        if bands[i] == 220:
            map220 += amap_shift
            n220 += 1

map150 /= np.float(n150)
map90 /= np.float(n90)
map220 /= np.float(n220)






