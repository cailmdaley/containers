import numpy as np
import scipy
from spt3g import core, std_processing, gcp, dfmux

class BufferThings(object):
    def __init__(self, input='Data'):
        self.wiringmap = None
        self.input = input
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wiringmap = frame['WiringMap']
        if frame.type == core.G3FrameType.Scan:
            hkmap = frame['DfMuxHousekeeping']
            names = frame['RawTimestreams_I'].keys()
            frame['IsGood'] = core.G3MapInt()
            print notavariable
            for bolo in names:
                status = dfmux.HousekeepingForBolo(hkmap,self.wiringmap,bolo)
                if status.state == 'None' and status.carrier_amplitude != 0.:
                    frame['IsGood'][bolo] = 1
                else:
                    frame['IsGood'][bolo] = 0

def make_bt1(frame,bt1 = None):
    bt1 = BufferThings()

def add_wiring(frame,bt1 = None):
    if frame.type == core.G3FrameType.Wiring:
        bt1(frame)
#        print notavariable

def check_status(frame,bt1 = None):
    if frame.type == core.G3FrameType.Scan:
        bt1(frame)
#        print notavariable

def calc_calfreq(frame):

    if frame.type != core.G3FrameType.Scan:
        return
#    try:
    calon = frame['CalibratorOn']
    srate_hz = calon.sample_rate/core.G3Units.Hz
    calon = np.asarray(calon)
    npts = len(calon)
    freqs = np.arange(npts/2)/(npts/2.)*srate_hz/2.
    calon_f = np.abs(np.fft.fft(calon-np.mean(calon)))
    calfreq = freqs[np.argmax(calon_f[0:npts/2])]
    frame['CalibratorFrequency'] = calfreq
#    except:
#        core.log_warn('No cal sync signal found!\n')
    print notavariable
        
def calc_ampls(frame):

    if frame.type != core.G3FrameType.Scan:
        return
    rti = frame['RawTimestreams_I']
    names = rti.keys()
    nbolos = len(names)
#    frame['CalibratorResponse'] = core.G3MapDouble()
    cr = core.G3MapDouble()
    for name in names:
        if frame['IsGood'][name]:
            bolots = rti[name]
            npts = len(bolots)
            bolots -= np.mean(bolots)
            bolots[0:100] = 0.
            bolots[npts-100:] = 0.
            srate_hz = bolots.sample_rate/core.G3Units.Hz
            freqs = np.arange(npts/2)/(npts/2.)*srate_hz/2.
            whcal = np.argmin(freqs-frame['CalibratorFrequency'])
            bolots_f = np.abs(np.fft.fft(bolots-np.mean(bolots)))
            cr[name] = bolots_f[whcal]
            frame['CalibratorResponse'] = cr

file1 = '/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3'
outfile1 = '/home/tcrawfor/cal_out_3009897_v2.g3'

pipe1 = core.G3Pipeline()
pipe1.Add(make_bt1, bt1=bt1)
pipe1.Add(core.G3Reader,filename=file1)
pipe1.Add(add_wiring, bt1=bt1)
pipe1.Add(check_status, bt1=bt1)
pipe1.Add(calc_calfreq)
pipe1.Add(calc_ampls)
pipe1.Add(core.G3Writer,filename=outfile1)
pipe1.Run()
