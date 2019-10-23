import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 900
nx = 900

obsids = ['33148130','33152051','33155972','33159892']
nobs = len(obsids)
mdir = '/poleanalysis/tcrawfor/calresult/calibration/PMNJ0210-5101-pixelraster/'

maps = np.zeros([nobs,ny,nx])
for i in np.arange(nobs):
    thisfile = mdir + obsids[i] + '.g3'
    f1 = core.G3File(thisfile)
    frame = f1.next()
    frame = f1.next()
    mt1=np.asarray(frame['T'])
    wt1=np.asarray(frame['Wunpol'].TT)
    wh1=np.where(wt1 > 0)
    mt1[wh1]/=wt1[wh1]
    maps[i,:,:] = mt1

ymax = np.zeros(nobs)
xmax = np.zeros(nobs)
for i in np.arange(nobs):
    amap_sm = ndimage.gaussian_filter(maps[i,200:700,200:700],4.)
    ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[500,500])
    ymax[i] = ycenter + 200
    xmax[i] = xcenter + 200
