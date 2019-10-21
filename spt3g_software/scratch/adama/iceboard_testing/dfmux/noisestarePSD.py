#import dfmuxtools
from spt3g import core, dfmux
import numpy as np

noisefile = '/data/sptdaq/iceboard_testing/dfmux/dfmux_20151202_075912_to_20151202_102912.g3'
noisescanfile = 'noisescan_20151202_075912_to_20151202_102912.g3'

PSDs = dict()

nframes = 0
scanframes = []
def averagePSD(frame):
    global nframes, PSDs, myframe
    if frame.type == core.G3FrameType.Scan:
        print frame
        for key in frame['RawTimestreams_I'].keys():
            bolo_psd = np.abs(np.fft.fft(frame['RawTimestreams_I'][key]))
            if key not in PSDs.keys():
                PSDs[key] = bolo_psd
            else:
                PSDs[key] = PSDs[key] + bolo_psd
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
pipe.Add(averagePSD)
pipe.Run()



