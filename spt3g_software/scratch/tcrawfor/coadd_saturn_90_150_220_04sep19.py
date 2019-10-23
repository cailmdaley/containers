from spt3g import core, std_processing, pointing, mapmaker, coordinateutils
import numpy as np
import pickle
import scipy
from scipy import ndimage

#obsid = '83559260'
obsid = '85628116'

file2 = '/spt/user/production/calibration/boloproperties/60000000.g3'
f2 = core.G3File(file2)
bp = f2.next()['BolometerProperties']

mapshape = [360,360]
map90 = np.zeros(mapshape)
n90 = 0
map150 = np.zeros(mapshape)
n150 = 0
map220 = np.zeros(mapshape)
n220 = 0

file1 = '/spt/user/production/calibration/saturn-pixelraster/maps/'+obsid+'.g3'
f1 = core.G3File(file1)
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt[np.where(twt == 0.)] = 1.
        if frame['Id'] in bp.keys():
            map_weighted = np.asarray(frame['T'])
            if np.isfinite(bp[frame['Id']].band) and np.max(np.abs(map_weighted)) > 0.:
                map_unw = map_weighted.copy()
                map_unw /= twt
                name = frame['Id']
                amap = np.zeros(mapshape)
                amap[80:350,80:350] = -map_unw[80:350,80:350]*1e12
                amap_sm = ndimage.gaussian_filter(amap,4)
                if np.max(amap_sm)/np.std(amap_sm[330:,330:]) > 20.:
                    ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
                    amap_shift = np.roll(amap,np.int(mapshape[0]/2-ycenter),0)
                    amap_shift = np.roll(amap_shift,np.int(mapshape[1]/2-xcenter),1)
                    if np.int(bp[name].band/10.) == 150:
                        map150 += amap_shift
                        n150 += 1
                    if np.int(bp[name].band/10.) == 90:
                        map90 += amap_shift
                        n90 += 1
                    if np.int(bp[name].band/10.) == 220:
                        map220 += amap_shift
                        n220 += 1

map150 /= np.float(n150)
map90 /= np.float(n90)
map220 /= np.float(n220)






