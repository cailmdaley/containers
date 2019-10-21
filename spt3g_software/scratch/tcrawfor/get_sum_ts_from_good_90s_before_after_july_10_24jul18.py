from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

# run get_some_always... before this

npts = 3000

fobsids_highel = [
    '47124531',
    '47470890',
    '47481109',
    '47491320',
    '47698388',
    '47708655',
    '47718888',
    '47908944',
    '47919033',
    '47929134',
    '48159813',
    '48168442',
    '48177226',
    '48443027',
    '48451642',
    '48460318',
    '48629218',
    '48637906',
    '48646593',
    '48835311',
    '48843915',
    '48878273',
    '48886888',
    '48895577']

eobsids_highel = []
efiles=glob.glob('/spt/data/bolodata/downsampled/elnod/4???????/0000.g3')
eobsids_all_i = np.asarray([np.int((thisfile.split('/'))[-2].split('.')[0]) for thisfile in efiles])
for obsid in fobsids_highel:
    index = np.argmin(np.abs(eobsids_all_i - np.int(obsid)))
    eobsids_highel.append(str(eobsids_all_i[index]))
#efiles=glob.glob('/spt/data/bolodata/downsampled/elnod/4???????/0000.g3')
efiles = ['/spt/data/bolodata/downsampled/elnod/'+eobsid+'/0000.g3' for eobsid in eobsids_highel]
efiles.sort()
ets=np.zeros([len(efiles),npts])
etimes=[]
for i in np.arange(len(efiles)):
    f1=core.G3File(efiles[i])
    ntot = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            etimes.append(frame['ObservationStart'])
            print frame
        if frame.type is core.G3FrameType.Scan:
#            print frame
            thisn = len(frame['OnlineBoresightAz'])
            minn = ntot
            maxn = np.min([ntot+thisn,npts])
            ntot += maxn - minn
#            print minn
#            print maxn
#            print ntot
#            for name in goodnames:
            for name in goodnames_150:
                try:
                    ets[i,minn:maxn] += frame['RawTimestreams_I'][name][0:maxn-minn]
                except:
                    pass
#                    print(notavariable)
#    print(notavariable)

cobsids_highel = []
cfiles=glob.glob('/spt/data/bolodata/downsampled/calibrator/4???????/0000.g3')
cobsids_all_i = np.asarray([np.int((thisfile.split('/'))[-2].split('.')[0]) for thisfile in cfiles])
for obsid in fobsids_highel:
    index = np.argmin(np.abs(cobsids_all_i - np.int(obsid)))
    cobsids_highel.append(str(cobsids_all_i[index]))

cfiles = ['/spt/data/bolodata/downsampled/calibrator/'+cobsid+'/0000.g3' for cobsid in cobsids_highel]
cfiles.sort()
cts=np.zeros([len(cfiles),npts])
ctimes=[]
for i in np.arange(len(cfiles)):
    f1=core.G3File(cfiles[i])
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            ctimes.append(frame['ObservationStart'])
        if frame.type is core.G3FrameType.Scan:
#            for name in goodnames:
            for name in goodnames_150:
                try:
                    cts[i,:] += frame['RawTimestreams_I'][name][0:npts]
                except:
                    pass


#nfiles=glob.glob('/spt/data/bolodata/downsampled/noise/4???????/0000.g3')
#nfiles.sort()
#nts=np.zeros([len(nfiles),npts])
#ntimes=[]
#for i in np.arange(len(nfiles)):
#    f1=core.G3File(nfiles[i])
#    for frame in f1:
#        if frame.type is core.G3FrameType.Observation:
#            ntimes.append(frame['ObservationStart'])
#        if frame.type is core.G3FrameType.Scan:
#            for name in goodnames:
#                try:
#                    nts[i,:] += frame['RawTimestreams_I'][name][0:npts]
#                except:
#                    pass

