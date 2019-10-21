from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle, glob

# run get_some_always... before this

#npts = 3000
#
#files_highel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-67.25/*/0000.g3')
#files_highel.sort()
#fts_highel=np.zeros([len(files_highel),npts])
#ftimes_highel=[]
#
#nfmax = 1
#for i in np.arange(len(files_highel)):
#    f1=core.G3File(files_highel[i])
#    nfdone = 0
#    while nfdone < nfmax:
#        frame = f1.next()
#        if frame.type is core.G3FrameType.Observation:
#            ftimes_highel.append(frame['ObservationStart'])
#            print frame
#        if frame.type is core.G3FrameType.Scan:
#            if 'Turnaround' not in frame:
#                nfdone += 1
#                try:
#                    for name in goodnames_150:
#                        fts_highel[i,:] += frame['RawTimestreams_I'][name][0:npts]
#                except:
#                    pass

files_lowel = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-44.75/*/0000.g3')
files_lowel.sort()
fts_lowel=np.zeros([len(files_lowel),npts])
ftimes_lowel=[]

nfmax = 1
for i in np.arange(len(files_lowel)):
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
