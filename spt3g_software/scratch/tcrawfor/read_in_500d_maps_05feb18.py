import numpy as np
import scipy
from scipy import ndimage
import pickle, glob
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 450
nx = 750

ndhdir = '/poleanalysis/ndhuang/ra0hdec-57.5_auto/'
#files = [ndhdir + '34049133.g3',ndhdir + '34062650.g3']
files = glob.glob(ndhdir+'/3*.g3')
files.sort()
#files = files[1:]
nobs = len(files)

maps = np.zeros([nobs,3,ny,nx])
wts = np.zeros([nobs,3,ny,nx])
bigmaps = np.zeros([3,ny,nx])
bigwtss = np.zeros([3,ny,nx])
for i in np.arange(nobs):
    f1 = core.G3File(files[i])
    psframe = f1.next()
    for j in np.arange(3):
        frame = f1.next()
        mt1=np.asarray(frame['T'])
        bigmaps[j,:,:] += mt1
        wt1=np.asarray(frame['Wpol'].TT)
        bigwtss[j,:,:] += wt1
        wh1=np.where(wt1 > 0)
        mt1[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt1
        wts[i,j,:,:] = wt1

bigmaps_wtd = bigmaps.copy()
for j in np.arange(3):
    wt2 = bigwtss[j,:,:]
    mt2 = bigmaps[j,:,:]
    wh2 = np.where(wt2 > 0.)
    mt2[wh2] /= wt2[wh2]
    bigmaps[j,:,:] = mt2
