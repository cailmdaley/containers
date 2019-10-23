import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 180
nx = 180

obsids = ['33395462','33397925','33400386','33405321','33410251','33415176','33420113','33402859','33407784','33412714','33417650','33422576']
obsids.sort()
nobs = len(obsids)
mdir = '/poleanalysis/tcrawfor/calresult/calibration/PMNJ0210-5101-pixelraster/half_arcmin/'

maps = np.zeros([nobs,3,ny,nx])
for i in np.arange(nobs):
#for i in np.arange(1)+11:
    thisfile = mdir + obsids[i] + '.g3'
    f1 = core.G3File(thisfile)
    for j in np.arange(3):
        frame = f1.next()
        mt1=np.asarray(frame['T'])
        wt1=np.asarray(frame['Wunpol'].TT)
        wh1=np.where(wt1 > 0)
        mt2 = mt1.copy()
        mt2[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt2

