import numpy as np
from scipy import ndimage
from spt3g import core, mapmaker
import centroid_maybe as cent
from matplotlib import pyplot as pl

def centroid_maps(mapfile):
    f = core.G3File(mapfile)
    fr = f.next()
    source = fr['Id']
    centers = {}
    for i, fr in enumerate(f):
        # if not i == 0:
        #     continue
        bolo = fr['Id']
        res = fr['T'].res / core.G3Units.arcmin
        map = np.array(fr['T'])
        weights = np.array(fr['Wunpol'].TT)
        if len(np.nonzero(weights)[0]) == 0:
            continue
        mfilt = ndimage.gaussian_filter(map, 5)
        xcen, ycen = np.unravel_index(np.nanargmax(mfilt), np.shape(map))
        dpix = 50
        ymin = ycen - dpix
        ymax = ycen + dpix
        xmin = xcen - dpix
        xmax = xcen + dpix
        if ymin < 0:
            ymin = 0
        if xmin < 0:
            xmin = 0
        if ymax >= np.shape(map)[0]:
            ymax = np.shape(map)[0] - 1
        if xmax >= np.shape(map)[1]:
            xmax = np.shape(map)[1] - 1
        map = map[xmin:xmax, ymin:ymax]
        map[np.isnan(map)] = 0.
        
        centers[bolo] = (cent.fit_map(map, res), (xcen, ycen), map)
    return centers

def denanify(map, weightmap):
    if len(np.nonzero(weightmap)[0]) == 0:
        raise RuntimeError('No good data')
    mapfilt = ndimage.gaussian_filter(map, 5)
    yp, xp = np.unravel_index(np.nanargmax(mapfilt), np.shape(weightmap))
    goodx = np.where(weightmap[yp] > 0)[0]
    goody = np.where(weightmap[:, xp] > 0)[0]
    return (goody[3], goody[-4]), (goodx[3], goody[-4])

if __name__ == '__main__':
    import sys
    offs = centroid_maps(sys.argv[1])
    np.savez(sys.argv[2], offs = offs)
