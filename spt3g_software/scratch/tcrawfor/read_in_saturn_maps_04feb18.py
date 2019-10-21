import numpy as np
import scipy
from scipy import ndimage
import pickle, glob
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 900
nx = 900

reso_arcmin = 0.5

file1 = '/poleanalysis/tcrawfor/calresult/calibration/saturn/34354304.g3'

maps = np.zeros([3,ny,nx])
f1 = core.G3File(file1)
for j in np.arange(3):
    frame = f1.next()
    mt1=np.asarray(frame['T'])
    if j == 0: wm90 = mt1
    if j == 1: wm150 = mt1
    if j == 2: wm220 = mt1
    wt1=np.asarray(frame['Wunpol'].TT)
    if j == 1: wt150 = np.asarray(frame['Wunpol'].TT)
    wh1=np.where(np.abs(wt1) > 0)
    mt1[wh1]/=wt1[wh1]
    maps[j,:,:] = mt1

