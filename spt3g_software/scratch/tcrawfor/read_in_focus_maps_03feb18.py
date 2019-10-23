import numpy as np
import scipy
from scipy import ndimage
import pickle, glob
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 900
nx = 900

files = glob.glob('/big_scratch/tcrawfor/calibration/saturn-pixelraster/coadds/3*.g3')
files.sort()
nobs = len(files)

maps = np.zeros([nobs,3,ny,nx])
benchoffs = np.zeros([nobs,6])
for i in np.arange(nobs):
    f1 = core.G3File(files[i])
    for j in np.arange(3):
        frame = f1.next()
        mt1=np.asarray(frame['T'])
        wt1=np.asarray(frame['Wunpol'].TT)
        wh1=np.where(wt1 > 0)
        mt1[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt1

