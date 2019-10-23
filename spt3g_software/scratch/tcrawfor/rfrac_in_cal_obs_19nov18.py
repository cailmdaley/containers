from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

files1 = np.asarray(glob.glob('/spt/user/production/calibration/calibrator/????????.g3'))
files1.sort()
#files1 = files1[0:len(files1)/3]
#files1 = files1[-100:]
files1 = files1[2000:]
obsids = np.asarray([np.float(file1.split('/')[-1].split('.')[0]) for file1 in files1])
whbefore = np.where(obsids < 40000000.)
files1 = files1[whbefore]

## !!!
#files1 = files1[0:40:10]
## !!!

#print(notavariable)

dtemp = pickle.load(open('/home/tcrawfor/gnames_19nov18.pkl'))
gnames = dtemp['goodnames']
gnames = gnames[0:100]
rfdict = {}

for file1 in files1:
    obsid_str = file1.split('/')[-1].split('.')[0]
    rfdict[obsid_str] = np.zeros(len(gnames))

    file1b = '/spt/data/bolodata/downsampled/calibrator/'+obsid_str+'/0000.g3'
    print file1b
    try:
        f1b = core.G3File(file1b)
        for frame in f1b:
#        print frame
            if frame.type is core.G3FrameType.Wiring:
                wiringmap = frame['WiringMap']
            if frame.type is core.G3FrameType.Scan:
                hkmap = frame['DfMuxHousekeeping']
                q = 0
                for name in gnames:
                    try:
                        status = dfmux.HousekeepingForBolo(hkmap,wiringmap,name)
                        rfdict[obsid_str][q] = status.rfrac_achieved
                        q += 1
                    except:
                        pass
    except:
        pass
