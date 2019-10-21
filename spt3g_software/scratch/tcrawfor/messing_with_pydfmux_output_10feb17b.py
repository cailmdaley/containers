import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

dirs = glob.glob('/poleanalysis/pydfmux_output/20170209/20170209_0*take_rawdump*/')
dirs.sort()
dirs = dirs[12:]
dirs.append('/poleanalysis/pydfmux_output/20170210/20170210_012318_take_rawdump/')
dirs.append('/poleanalysis/pydfmux_output/20170210/20170210_013215_take_rawdump/')
dirs.sort()

rdave = np.zeros([len(dirs),24999])
for i in np.arange(len(dirs)):
    tfiles = glob.glob(dirs[i]+'data/IceCrate*.pkl')
    for tfile in tfiles:
        d1temp = pickle.load(open(tfile))
        rdave[i,:] += d1temp['freq_domain']['y']/np.float(len(tfiles))
