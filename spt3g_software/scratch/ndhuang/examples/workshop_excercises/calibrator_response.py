import numpy as np
from scipy import signal
from spt3g import core

def AddCalResponse(fr, cal_freq = 1.0 * core.G3Units.Hz,
                   ts_key = 'RawTimestreams_I',
                   response_key = 'CalibratorResponse'):
    '''
    Add a vector of G3Doubles that give the approximate calibrator response
    of each bolometer.
    `cal_freq` is the frequency of the calibrator chopper, in native units
    `ts_key` is the key that points to the detector data
    `response_key` is the name of the new object that is added to the frame
    
    Some additional notes on the arguments: cal_freq is in "native"
    units.  This is preferable, as it allows the user to choose the
    units.  For example, one could look for calibrator signal at 25
    MHz by running 
    pipe.Add(AddCalResponse, cal_freq = 25 * core.G3Units.MHz)

    ts_key and response_key are both provided for flexibility.  It is
    entirely reasonable (almost expected) that someone will want to
    run analysis on other timestreams (such as the Q timestream from
    the iceboards).  Allowing the output key to be changed adds
    additional flexibility for additional modules further down the
    pipe.    
    '''
    if fr.type != core.G3FrameType.Scan:
        # Skip frames that are not of the Scan type (such as Wiring frames)
        return
    cal_response = core.G3MapDouble()
    for boloname, timestream in fr[ts_key].iteritems():
        # Super lazy calibrator analysis: 
        # just find the response at the frequency nearest to the 
        # cal frequency
        psd = np.abs(np.fft.rfft(signal.detrend(timestream)))**2
        dt = 1 / timestream.sample_rate
        freq = np.fft.rfftfreq(len(timestream), dt)
        cal_freq_ind = np.argmin(np.abs(freq - cal_freq))
        cal_response[boloname] = psd[cal_freq_ind]
    fr[response_key] = cal_response

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = '/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3')
pipe.Add(core.Dump)
pipe.Add(AddCalResponse, cal_freq = 6.0 * core.G3Units.Hz)
pipe.Add(core.Dump, type = core.G3FrameType.Scan, 
         added_message = 'Note that a new key has been added')
pipe.Add(core.G3Writer, filename = '/scratch/ndhuang/cal_response_out.g3')
pipe.Run()
