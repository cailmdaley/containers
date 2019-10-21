'''
Test script for playing around at the 3g software meeting
'''
import matplotlib.pyplot as plt
import numpy as np

from spt3g import core
from spt3g import dfmux

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


class GetCal(object):
    def __init__(self, input='Data'):
        self.wiringmap = None
        self.input = input

    def __call__(self, frame):

        # Get wiring map
        if frame.type == core.G3FrameType.Wiring:
            self.wiringmap = frame['WiringMap']

        # For Scan frames, calculate the Cal response
        if frame.type == core.G3FrameType.Scan:
            hkmap = frame['DfMuxHousekeeping']

            cal_values = core.G3MapDouble() # Output object
            for bolo, ts in frame['RawTimestreams_I']:
                status = dfmux.HousekeepingForBolo(hkmap, self.wiringmap , bolo)

                # identify only good bolos
                bolo_is_good =  ((status.state == 'tuned') and 
                                 (status.carrier_amplitude != 0))

                if bolo_is_good:
                    cal_values[bolo] = find_cal_one_bolo(ts)
                    print("Add bolo cal: %s, %f"%(bolo, cal_values[bolo]))
                else:
                    print("Bad bolo %s:"%bolo, status.state, status.carrier_amplitude)
                    cal_values[bolo] = 0.

            # Add Cal information to the frame
            frame['Cal'] = cal_values




## Build a pipe
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3')
pipe.Add(GetCal)
pipe.Add(core.G3Writer, filename='/scratch/kstory/out_0404.g3')

# Run the pipe
pipe.Run()


# Plot a histogram
data = core.G3File('/scratch/kstory/out_0404.g3')
tmp = data.next()
print(tmp.type)
while tmp.type != core.G3FrameType.Scan:
    print(tmp.type)
    tmp = data.next()

y = np.asarray(tmp['Cal'].values(), dtype='float')
y = y[y != 0]  # remove zeros

plt.hist(y)

#core.G3VectorString
# Make new G3TimestreamMap that only includes the good bolos

# Debugging
'''
frames = list(core.G3File('/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3'))
bolo = 'W136/2017.W136.5.4.7.912'
wiringmap = frames[1]['WiringMap']
hkmap = frames[-1]['DfMuxHousekeeping']
'''
