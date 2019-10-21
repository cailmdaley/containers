from spt3g import core, dfmux, mapmaker
import numpy as np
import pickle
import scipy
from scipy import ndimage

nx = 360
ny = 360
reso_arcmin = 0.5
mapdict = {}

file1 = '/spt/user/production/calibration/saturn-pixelraster/maps/35315694.g3'
f1 = core.G3File(file1)
maxframe = 20
qq = 0
for frame in f1:
    if frame.type is core.G3FrameType.Map:
        if 'Wunpol' in frame and frame['Id'] != 'bs':
            twt = np.asarray(frame['Wunpol'].TT)
            twt_orig = twt.copy()
            twt[np.where(twt == 0.)] = 1.
        map_weighted = np.asarray(frame['T'])
        if np.max(np.abs(map_weighted)) > 0.:
            map_unw = map_weighted.copy()
            map_unw /= twt
            mapdict[frame['Id']] = map_unw
            qq+=0
            if qq > maxframe:
                print(notavariable)

print(notavariable)

# now make some composite maps
#   hmmm, why isn't there a calframe for this file, and why doesn't the calibration frame have good bolo props data? ok, use an old calframe
cframe = (core.G3File('/spt/user/production/calibration/calframe/saturn-pixelraster/35052899.g3')).next()
bp = cframe['BolometerProperties']
# pick out high-S/N bolos
whw = np.where(twt > 2.)
bands = np.asarray([900.,1500.,2200.])
pixpos = np.asarray(['x','y'])
mapnames = ['90x','90y','150x','150y','220x','220y']
bigmap = {}
for mn in mapnames:
    bigmap[mn] = np.zeros([ny,nx])
for name in mapdict.keys():
    bpt = bp[name]
    thispixpos = bpt.physical_name[-1]
    thisband = bpt.band
    thismap = -mapdict[name]
    thissn = np.max(thismap[whw])/np.std(thismap[whw])
    if thissn > 10:
        whb = np.where(bands == thisband)[0][0]
        whp = np.where(pixpos == thispixpos)[0][0]
        kmap = 2*whb + whp
        ycenter, xcenter = np.unravel_index(np.argmax(np.abs(thismap)),[ny,nx])
        bigmap[mapnames[kmap]][ycenter-20:ycenter+20,xcenter-20:xcenter+20] += thismap[ycenter-20:ycenter+20,xcenter-20:xcenter+20]

for mn in mapnames:
    figure()
    imshow(bigmap[mn],vmin=-1e-12,vmax=1e-12,cmap='bone') 
    title(mn)
