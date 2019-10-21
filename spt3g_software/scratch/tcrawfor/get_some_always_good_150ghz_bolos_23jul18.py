from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

files1 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/4???????.g3'))
files1.sort()
#files1 = files1[0:len(files1)/3]
files1 = files1[-100:]
files1 = files1[::10]

ngoodthresh = 8

for i in np.arange(len(files1)):
    file1 = files1[i]
    print file1
    obsid_str = file1.split('/')[-1].split('.')[0]
    file1b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/nominal_online_cal.g3'
    f2 = core.G3File(file1b)
    bp = f2.next()['NominalBolometerProperties']
    cftemp1 = core.G3File(file1).next()
    stemp1 = cftemp1['CalibratorResponse']
    sntemp1 = cftemp1['CalibratorResponseSN']
    if i == 0:
        sndict = {}
        names = sntemp1.keys()
        for name in names:
            if np.int(bp[name].band/10.) == 150:
                sndict[name] = np.zeros(len(files1))
    for name in names:
        try:
            sndict[name][i] = sntemp1[name]
        except:
            pass

ngood={}
for name in sndict.keys(): 
    ngood[name]=len(np.where(sndict[name] > 200.)[0])
goodnames_150=np.asarray([name for name in sndict.keys() if ngood[name] > ngoodthresh])
