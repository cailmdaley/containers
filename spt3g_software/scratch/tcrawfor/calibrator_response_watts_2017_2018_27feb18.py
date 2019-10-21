from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

files1 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/19*.g3'))
files1.sort()
files1 = files1[::10]

files2 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/36*.g3'))
files2.sort()
files2 = files2[::10]

nhbins = 200
hrange = [0,20]
histarr1 = np.zeros([len(files1),3,nhbins])
histarr2 = np.zeros([len(files2),3,nhbins])

bands = [90,150,220]

for i in np.arange(len(files1)):
    file1 = files1[i]
    obsid_str = file1.split('/')[-1].split('.')[0]
    file1b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/nominal_online_cal.g3'
    f2 = core.G3File(file1b)
    bp = f2.next()['NominalBolometerProperties']
    cftemp1 = core.G3File(file1).next()
    stemp1 = cftemp1['CalibratorResponse']
    for j in np.arange(3):
        dtemp = np.asarray([stemp1[name] for name in stemp1.keys() if np.int(bp[name].band/10.) == bands[j]])
        htemp = np.histogram(dtemp*1e15,bins=nhbins,range=hrange)
        histarr1[i,j,:] = htemp[0]

for i in np.arange(len(files2)):
    file2 = files2[i]
    obsid_str = file2.split('/')[-1].split('.')[0]
    file2b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/nominal_online_cal.g3'
    f2 = core.G3File(file2b)
    bp = f2.next()['NominalBolometerProperties']
    cftemp2 = core.G3File(file2).next()
    stemp2 = cftemp2['CalibratorResponse']
    for j in np.arange(3):
        dtemp = np.asarray([stemp2[name] for name in stemp2.keys() if np.int(bp[name].band/10.) == bands[j]])
        htemp = np.histogram(dtemp*1e15,bins=nhbins,range=hrange)
        histarr2[i,j,:] = htemp[0]

#for i in np.arange(len(files2)):
#    file2 = files2[i]
##    obsid_str = file2.split('/')[-1].split('.')[0]
##    file1b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/nominal_online_cal.g3'
##    f2 = core.G3File(file1b)
##    bp = f2.next()['NominalBolometerProperties']
##    names = bp.keys()
##    nbolo = len(names)
#    cftemp1 = core.G3File(file2).next()
#    dtemp = np.asarray(cftemp1['CalibratorResponse'].values())
#    htemp = np.histogram(dtemp*1e15,bins=nhbins,range=hrange)
#    histarr2[i,:] = htemp[0]


