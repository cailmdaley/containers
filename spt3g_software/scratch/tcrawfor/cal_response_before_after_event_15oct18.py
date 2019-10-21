from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

files1 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/????????.g3'))
files1.sort()
#files1 = files1[0:len(files1)/3]
#files1 = files1[-100:]
files1 = files1[2000:]
files1 = files1[::50]
#files1 = files1[::10]
obsids = np.asarray([np.float(file1.split('/')[-1].split('.')[0]) for file1 in files1])
whbefore = np.where(obsids < 40000000.)
maxbefore = np.max(whbefore)

for i in np.arange(len(files1)):
    file1 = files1[i]
    print file1
    obsid_str = file1.split('/')[-1].split('.')[0]
    file1b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/nominal_online_cal.g3'
    try:
        f2 = core.G3File(file1b)
        bp = f2.next()['NominalBolometerProperties']
    except:
        pass
    cftemp1 = core.G3File(file1).next()
    stemp1 = cftemp1['CalibratorResponse']
    sntemp1 = cftemp1['CalibratorResponseSN']
    if i == 0:
        sndict = {}
        sdict = {}
        names = sntemp1.keys()
        for name in names:
            if np.int(bp[name].band/10.) == 150:
                sndict[name] = np.zeros(len(files1))
                sdict[name] = np.zeros(len(files1))
    for name in names:
        try:
            sndict[name][i] = sntemp1[name]
            sdict[name][i] = stemp1[name]
        except:
            pass

ratio_dict = {}
for name in sndict.keys():
    whg = np.where(sndict[name] >= 20.)
    if len(whg[0]) > 30:
        whg2 = np.where(sndict[name][0:maxbefore+1] >= 20.)
        if len(whg2[0]) > 5:
            vec1 = sdict[name][0:maxbefore+1]
            vec2 = sdict[name][maxbefore+1:]
            meanresp1 = np.median(vec1[whg2])
            whg3 = np.where(sndict[name][12:] >= 20.)
            meanresp2 = np.median(vec2[whg3])
            ratio_dict[name] = meanresp2/meanresp1
