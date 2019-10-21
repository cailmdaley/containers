#import dfmuxtools                                                                                                                                                      
from spt3g import core, dfmux
import numpy as np
import cPickle as pickle

noisefile = '/data/sptdaq/iceboard_testing/dfmux/dfmux_20151202_075912_to_20151202_102912.g3'
noisescanfile = 'noisescan_20151202_075912_to_20151202_102912.g3'

PSDs = dict()

nframes = 0
scanframes = []
def savescans(frame):
    global nframes, PSDs, myframe
    if frame.type == core.G3FrameType.Scan:
        print frame
        nframes += 1
        scanframes.append(frame)

def skip_frames(frame):
    global nframes
    if nframes >= 100:
        return []

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=noisefile)
pipe.Add(skip_frames)
pipe.Add(dfmux.FixedLengthScans, N=2048)
pipe.Add(dfmux.DfMuxCollator)
pipe.Add(savescans)
pipe.Run()

f = open('scans_motors_on.pkl', 'w')
pickle.dump(scanframes, f)
file.close()
