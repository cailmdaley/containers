'''
Test script for playing around at the 3g software meeting
'''
import matplotlib.pyplot as plt
import numpy as np
import sys

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




if __name__ == '__main__':
    g3file = sys.argv[1]
    
    ## Build a pipe
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=g3file)
    pipe.Add(GetCal)
    pipe.Add(core.G3Writer, filename='calibrator_grid_0404.g3')

    # Run the pipe
    pipe.Run()

