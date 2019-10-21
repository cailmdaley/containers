from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

rf8dict = pickle.load(open('/home/tcrawfor/rf8_obsids.pkl'))
n8_before = len(rf8dict['obsids'])

files1 = []
for obsid in rf8dict['obsids']:
    files1.append('/spt/user/production/calibration/calibrator/'+obsid+'.g3')

files2 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/????????.g3'))
files2.sort()

obsids = np.asarray([np.float(file2.split('/')[-1].split('.')[0]) for file2 in files2])
whbefore_orig = np.where(obsids < 40000000.)
maxbefore_orig = np.max(whbefore_orig)
files3 = files2[maxbefore_orig+1:]
nafter = len(files3)
cadence = nafter/n8_before
for i in np.arange(n8_before):
    files1.append(files3[i*cadence])
files1.sort()
whbefore = range(n8_before)
maxbefore = n8_before-1

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
#            if np.int(bp[name].band/10.) == 150:
            sndict[name] = np.zeros(len(files1))
            sdict[name] = np.zeros(len(files1))
    if np.int(np.round(cftemp1['CalibratorResponseFrequency']/core.G3Units.Hz)) == 4:
        for name in names:
            try:
                sndict[name][i] = sntemp1[name]
                sdict[name][i] = stemp1[name]
            except:
                pass

#snthresh = 100.
snthresh = 40.
ratio_dict = {}
for name in sndict.keys():
    if name in bp.keys():
        whg = np.where(sndict[name] >= snthresh)
        if len(whg[0]) > 30:
            whg2 = np.where(sndict[name][0:maxbefore+1] >= snthresh)
            if len(whg2[0]) > 5:
                vec1 = sdict[name][0:maxbefore+1]
                vec2 = sdict[name][maxbefore+1:]
                meanresp1 = np.nanmedian(vec1[whg2])
                whg3 = np.where(sndict[name][maxbefore+1:] >= snthresh)
                meanresp2 = np.nanmedian(vec2[whg3])
                ratio_dict[name] = meanresp2/meanresp1

ratios150 = np.asarray([ratio_dict[name] for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 150])
xp150 = np.asarray([bp[name].x_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 150])
yp150 = np.asarray([bp[name].y_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 150])
ratios90 = np.asarray([ratio_dict[name] for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 90])
xp90 = np.asarray([bp[name].x_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 90])
yp90 = np.asarray([bp[name].y_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 90])
ratios220 = np.asarray([ratio_dict[name] for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 220])
xp220 = np.asarray([bp[name].x_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 220])
yp220 = np.asarray([bp[name].y_offset for name in ratio_dict.keys() if np.int(bp[name].band/10.) == 220])
