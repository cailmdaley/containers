from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

# run get_some_always... before this

npts = 3000

files_highel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-67.25/*/0000.g3')
files_highel.sort()
fts_highel=np.zeros([len(files_highel),npts])
ftimes_highel=[]
cts_highel=np.zeros([len(files_highel),npts])

nfmax = 1
for i in np.arange(len(files_highel)):
    print files_highel[i]
    f1=core.G3File(files_highel[i])
    nfdone = 0
    while nfdone < nfmax:
        frame = f1.next()
        if frame.type is core.G3FrameType.Observation:
            ftimes_highel.append(frame['ObservationStart'])
            print frame
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
                nfdone += 1
                try:
                    for name in goodnames_150:
                        fts_highel[i,:] += frame['RawTimestreams_I'][name][0:npts]
                except:
                    pass
    thiscobs = tctools.get_cal_obs_before_field_scan(files_highel[i])
    thiscfile = '/spt/data/bolodata/downsampled/calibrator/' + thiscobs + '/0000.g3'
    f2 = core.G3File(thiscfile)
    for frame in f2:
        if frame.type is core.G3FrameType.Scan:
            try:
                for name in goodnames_150:
                    cts_highel[i,:] += frame['RawTimestreams_I'][name][0:npts]
            except:
                pass
        
fts_highel /= np.float(len(goodnames_150))
cts_highel /= np.float(len(goodnames_150))


files_lowel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-44.75/*/0000.g3')
files_lowel.sort()
fts_lowel=np.zeros([len(files_lowel),npts])
ftimes_lowel=[]
cts_lowel=np.zeros([len(files_lowel),npts])

nfmax = 1
for i in np.arange(len(files_lowel)):
    print files_lowel[i]
    f1=core.G3File(files_lowel[i])
    nfdone = 0
    while nfdone < nfmax:
        frame = f1.next()
        if frame.type is core.G3FrameType.Observation:
            ftimes_lowel.append(frame['ObservationStart'])
            print frame
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
                nfdone += 1
                try:
                    for name in goodnames_150:
                        fts_lowel[i,:] += frame['RawTimestreams_I'][name][0:npts]
                except:
                    pass
    thiscobs = tctools.get_cal_obs_before_field_scan(files_lowel[i])
    thiscfile = '/spt/data/bolodata/downsampled/calibrator/' + thiscobs + '/0000.g3'
    f2 = core.G3File(thiscfile)
    for frame in f2:
        if frame.type is core.G3FrameType.Scan:
            try:
                for name in goodnames_150:
                    cts_lowel[i,:] += frame['RawTimestreams_I'][name][0:npts]
            except:
                pass
        
fts_lowel /= np.float(len(goodnames_150))
cts_lowel /= np.float(len(goodnames_150))
