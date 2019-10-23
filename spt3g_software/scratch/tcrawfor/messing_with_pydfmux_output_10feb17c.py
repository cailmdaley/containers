import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

dirs = []
dirs.append('/spt/data/pydfmux_output/20170212/20170212_103952_take_rawdump/')
dirs.append('/spt/data/pydfmux_output/20170212/20170212_110110_take_rawdump/')
dirs.sort()

rdave = np.zeros([len(dirs),24999])
for i in np.arange(len(dirs)):
    tfiles = glob.glob(dirs[i]+'data/IceCrate*.pkl')
    for tfile in tfiles:
        d1temp = pickle.load(open(tfile))
        rdave[i,:] += d1temp['freq_domain']['y']/np.float(len(tfiles))
