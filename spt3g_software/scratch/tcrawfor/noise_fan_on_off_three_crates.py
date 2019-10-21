import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core
import os
import glob

dir1 = '/scratch/pydfmux_output/'

files = {}
ckeys = ['005','015','016']
for key in ckeys:
    files[key] = {}

files['015']['fan_off'] = glob.glob(dir1 + '20170203/20170202_205811_measure_noise/data/*.pkl')
files['015']['fan_on'] = glob.glob(dir1 + '20170203/20170202_205953_measure_noise/data/*.pkl')
files['005']['fan_off'] = glob.glob(dir1 + '20170206/20170206_083106_measure_noise_pulse_tubes_off/data/*.pkl')
files['005']['fan_on'] = glob.glob(dir1 + '20170206/20170206_082749_measure_noise/data/*.pkl')
files['016']['fan_off'] = glob.glob(dir1 + '20170306/20170306_080758_measure_noise_test3e2_benchtop_config/data/*.pkl')
files['016']['fan_on'] = glob.glob(dir1 + '20170307/20170306_123605_measure_noise_test3d2_crate_cabin_door_fan_tray/data/*.pkl')

fbdict = {}
noisedict = {}

for key1 in files.keys():
    fbdict[key1] = {}
    noisedict[key1] = {}
    for key2 in files[key1].keys():
        file1 = files[key1][key2][0]
        d1 = pickle.load(open(file1))
        fbdict[key1][key2] = np.asarray([d1[key]['frequency'] for key in d1.keys()])
        noisedict[key1][key2] = np.asarray([d1[key]['noise']['median_noise'] for key in d1.keys()])


