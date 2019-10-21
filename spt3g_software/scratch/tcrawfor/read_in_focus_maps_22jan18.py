import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 900
nx = 900

obsids = ['33395462','33397925','33400386','33405321','33410251','33415176','33420113','33402859','33407784','33412714','33417650','33422576']
obsids.sort()
#obsids = ['33244105','33248876','33253653','33258424','33263201','33267971']
#obsids = ['33248876','33253653','33258424','33263201','33267971']
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
