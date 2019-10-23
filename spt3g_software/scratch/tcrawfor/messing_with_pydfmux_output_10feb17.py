import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

allfiles_0210 = glob.glob('/poleanalysis/pydfmux_output/20170210/*measure_noise*/data/*.pkl')
allfiles_0210.sort()
allfiles_0209 = glob.glob('/poleanalysis/pydfmux_output/20170209/*measure_noise*/data/*.pkl')
allfiles_0209.sort()
allfiles = allfiles_0209[10:]
for file1 in allfiles_0210:
    allfiles.append(file1)



fbdict = {}
noisedict = {}

for file1 in allfiles:
    ff1 = filter(None,file1.split('/'))
    ff2 = filter(None,ff1[3].split('_'))
    thistime = core.G3Time(ff2[0] + '_' + ff2[1])
    d1 = pickle.load(open(file1))
    fbdict[thistime] = np.asarray([d1[key]['frequency'] for key in d1.keys()])
    noisedict[thistime] = np.asarray([d1[key]['noise']['median_noise'] for key in d1.keys()])


