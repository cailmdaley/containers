from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

# run get_some_always... before this

npts = 3000

cfiles=glob.glob('/spt/data/bolodata/downsampled/calibrator/4???????/0000.g3')
cfiles.sort()
cts=np.zeros([len(cfiles)/10,npts])
ctimes=[]
for i in np.arange(len(cfiles)/10):
    f1=core.G3File(cfiles[10*i])
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            ctimes.append(frame['ObservationStart'])
        if frame.type is core.G3FrameType.Scan:
            for name in goodnames:
                try:
                    cts[i,:] += frame['RawTimestreams_I'][name][0:npts]
                except:
                    pass

efiles=glob.glob('/spt/data/bolodata/downsampled/elnod/4???????/0000.g3')
efiles.sort()
ets=np.zeros([len(efiles)/5,npts])
etimes=[]
for i in np.arange(len(efiles)/5):
    f1=core.G3File(efiles[5*i])
    ntot = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            etimes.append(frame['ObservationStart'])
        if frame.type is core.G3FrameType.Scan:
            print frame
            thisn = len(frame['OnlineBoresightAz'])
            minn = ntot
            maxn = np.min([ntot+thisn,npts])
            ntot += maxn - minn
            print minn
            print maxn
            print ntot
            for name in goodnames:
                try:
                    ets[i,minn:maxn] += frame['RawTimestreams_I'][name][0:maxn-minn]
                except:
                    pass
#                    print(notavariable)
#    print(notavariable)

nfiles=glob.glob('/spt/data/bolodata/downsampled/noise/4???????/0000.g3')
nfiles.sort()
nts=np.zeros([len(nfiles),npts])
ntimes=[]
for i in np.arange(len(nfiles)):
    f1=core.G3File(nfiles[i])
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            ntimes.append(frame['ObservationStart'])
        if frame.type is core.G3FrameType.Scan:
            for name in goodnames:
                try:
                    nts[i,:] += frame['RawTimestreams_I'][name][0:npts]
                except:
                    pass

