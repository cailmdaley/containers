from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

# run get_some_always... before this

npts = 3000

files_highel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-67.25/*')
files_highel.sort()
fobsids_highel = np.asarray([np.int((thisfile.split('/'))[-1]) for thisfile in files_highel])
eobsids_highel = []
efiles=glob.glob('/spt/data/bolodata/downsampled/elnod/*/0000.g3')
eobsids_all_i = np.asarray([np.int((thisfile.split('/'))[-2].split('.')[0]) for thisfile in efiles])
for obsid in fobsids_highel:
    index = np.argmin(np.abs(eobsids_all_i - np.int(obsid)))
    eobsids_highel.append(str(eobsids_all_i[index]))
efiles_highel = ['/spt/data/bolodata/downsampled/elnod/'+eobsid+'/0000.g3' for eobsid in eobsids_highel]
efiles_highel.sort()
ets_highel=np.zeros([len(efiles_highel),npts])
etimes_highel=[]

for i in np.arange(len(efiles_highel)):
    f1=core.G3File(efiles_highel[i])
    ntot = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            etimes_highel.append(frame['ObservationStart'])
            print frame
        if frame.type is core.G3FrameType.Scan:
            thisn = len(frame['OnlineBoresightAz'])
            minn = ntot
            maxn = np.min([ntot+thisn,npts])
            ntot += maxn - minn
#            for name in goodnames_150:
            for name in goodnames:
                try:
                    ets_highel[i,minn:maxn] += frame['RawTimestreams_I'][name][0:maxn-minn]
                except:
                    pass

files_lowel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-44.75/*')
files_lowel.sort()
fobsids_lowel = np.asarray([np.int((thisfile.split('/'))[-1]) for thisfile in files_lowel])
eobsids_lowel = []
efiles=glob.glob('/spt/data/bolodata/downsampled/elnod/*/0000.g3')
eobsids_all_i = np.asarray([np.int((thisfile.split('/'))[-2].split('.')[0]) for thisfile in efiles])
for obsid in fobsids_lowel:
    index = np.argmin(np.abs(eobsids_all_i - np.int(obsid)))
    eobsids_lowel.append(str(eobsids_all_i[index]))
efiles_lowel = ['/spt/data/bolodata/downsampled/elnod/'+eobsid+'/0000.g3' for eobsid in eobsids_lowel]
efiles_lowel.sort()
ets_lowel=np.zeros([len(efiles_lowel),npts])
etimes_lowel=[]

for i in np.arange(len(efiles_lowel)):
    f1=core.G3File(efiles_lowel[i])
    ntot = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Observation:
            etimes_lowel.append(frame['ObservationStart'])
            print frame
        if frame.type is core.G3FrameType.Scan:
            thisn = len(frame['OnlineBoresightAz'])
            minn = ntot
            maxn = np.min([ntot+thisn,npts])
            ntot += maxn - minn
#            for name in goodnames_150:
            for name in goodnames:
                try:
                    ets_lowel[i,minn:maxn] += frame['RawTimestreams_I'][name][0:maxn-minn]
                except:
                    pass

