'''
Test script for playing around at the 3g software meeting
'''
import matplotlib.pyplot as plt
import numpy as np

from spt3g import core

def find_cal_one_bolo(timestream):

    # Get the PSD of this timestream
    timestream -= np.median(timestream) # remove median
    ft = np.abs(np.fft.fft(timestream))**2.

    # delta-time between samples
    nsamp = len(ft)
    dt = 1/(timestream.sample_rate / core.G3Units.Hz)
    f_vec = np.fft.rfftfreq(nsamp, d=dt)

    # Get power of 6 Hz line by searching for the max value between 1 and 10 Hz
    value = np.max(ft[np.where((f_vec > 1) & (f_vec < 10))])    

    return value


## Write a calibration module
def get_bolo_calibration(frame):

    # Only operate on scan frames
    if frame.type != core.G3FrameType.Scan:
        return

    tsmap = frame['RawTimestreams_I']  # timestreams

    cal_values = core.G3MapDouble() # Output object
    for name, ts in tsmap:
        cal_values[name] = find_cal_one_bolo(ts)

    # Add this to the frame
    frame['Cal'] = cal_values


## Build a pipe
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3')
pipe.Add(get_bolo_calibration)
pipe.Add(core.G3Writer, filename='/scratch/kstory/out_0403.g3')

# Run the pipe
pipe.Run()

# Debugging
'''
frames = list(core.G3File('/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3'))
timestream = frames[-1]['RawTimestreams_I'][ 'W136/2017.W136.5.4.7.912']
'''
