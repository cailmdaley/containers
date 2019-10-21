import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 900
nx = 900

obsids = ['34181029',  '34181570',  '34182109',  '34182648',  '34183187']
obsids.sort()
nobs = len(obsids)
mdir = '/poleanalysis/tcrawfor/calresult/calibration/saturn/'

maps = np.zeros([nobs,3,ny,nx])
benchoffs = np.zeros([nobs,6])
for i in np.arange(nobs):
    thistime = std_processing.obsid_to_g3time(obsids[i])
    
    thisfile = mdir + obsids[i] + '.g3'
    f1 = core.G3File(thisfile)
    for j in np.arange(3):
        frame = f1.next()
        mt1=np.asarray(frame['T'])
        wt1=np.asarray(frame['Wunpol'].TT)
        wh1=np.where(wt1 > 0)
        mt1[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt1

ymax = np.zeros([nobs,3])
xmax = np.zeros([nobs,3])
for i in np.arange(nobs):
    for j in np.arange(3):
        amap_sm = ndimage.gaussian_filter(maps[i,j,200:700,200:700],4.)
        ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[500,500])
        ymax[i,j] = ycenter + 200
        xmax[i,j] = xcenter + 200
