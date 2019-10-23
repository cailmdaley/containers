import scipy.fftpack
import cPickle as pickle
import numpy as np
import scipy.fftpack
from spt3g import core
import matplotlib.pyplot as plt

f = open('scans_motors_on.pkl')
scanframes = pickle.load(f)

ts_length = 2048

PSDs = dict()
for scan in scanframes:
    for key in scan['RawTimestreams_I'].keys():
        if key not in PSDs.keys():
            PSDs[key] = np.zeros(len(scan['RawTimestreams_I'][key]))
        f = scipy.fftpack.rfft(scan['RawTimestreams_I'][key]).real
        PSDs[key] += np.conj(f)*f/ len(scanframes)

        time = (scan['RawTimestreams_I'][key].stop.time - scan['RawTimestreams_I'][key].start.time) / core.G3Units.second
        freq = scipy.fftpack.rfftfreq(ts_length, time / ts_length)

for readoutmod in range(4):
    plt.figure(readoutmod)
    for key in PSDs.keys():
        if int(key[-1]) == readoutmod+1:
            plt.semilogy(freq, PSDs[key])
    plt.xlabel('frequency [Hz]')
    plt.xlim([np.min(freq), np.max(freq)])
    plt.title('readout module ' + str(readoutmod+1))        
plt.show()
            
