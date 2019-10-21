#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.todfilter import timeconstant_filter as tc
from spt3g.todfilter.convolve_filter import ConvolveFilter

class DeltaFunctionFrames(object):
    '''
    Creates frames with a user-input sample rate to test the convolution - deconvolution algorithms.  This class produces timestreams 
    with a delta function, to ensure that we have constructed a zero-lag, zero-leakage filter.
    '''
    def __init__(self, nsamples = 1001, delta_location = [23], bolos = ['a'], sample_rate = 100.*core.G3Units.Hz):
        self.nsamples = nsamples
        self.delta_loc = delta_location
        self.sample_rate = sample_rate
        self.bolos = bolos
        self.maxframe = len(delta_location)
        self.nframes = 0

    def __call__(self, fr):
        assert(fr is None)
        if self.nframes >= self.maxframe:
            return []
        ts = np.zeros(self.nsamples)
        if self.delta_loc[self.nframes] is not None:
            ts[self.delta_loc[self.nframes]] = 1024

        ts = core.G3Timestream(ts)
        ts.start = core.G3Time.Now()
        ts.stop = ts.start + (self.nsamples - 1.)*(1./self.sample_rate)
        assert(np.isclose(ts.sample_rate, self.sample_rate))
        tm = core.G3TimestreamMap()
        for b in self.bolos:
            tm[b] = ts

        outfr = core.G3Frame(core.G3FrameType.Scan)
        outfr['tods'] = tm
        outfr['Has Delta'] = (self.delta_loc[self.nframes] is not None)
        self.nframes += 1
        return outfr

# This function generates an exponential filter to simulate the effect of a time constant. 
# The filter length is set by requiring the value in the last bin to be smaller than rel_err
# for the longest time constant in the data set (with a minimum length of 5).

def exponential_filter(taus, sample_rate, rel_err):
    filter = {}
    mult = -np.log(rel_err) 
    klen = int(max(taus.values())*sample_rate*mult)*2 + 1 
    klen = max(klen, 5)
    ktimes = (np.arange(klen)-klen//2)/sample_rate
    for bolo in taus:
        kernel = np.where(ktimes >= 0., np.exp(-ktimes/taus[bolo]), 0.)
        filter[bolo] = kernel/kernel.sum()
    return filter

# Inject a given time constant dictionary into the frame
def inject_taus(frame, taus):
    if frame.type == core.G3FrameType.Scan:
        tau_md = core.G3MapDouble()
        for bolo, tau in taus.items():
            tau_md[bolo] = tau
        frame['TimeConstants'] = tau_md
    return

def compare_peak_location(frame, ts1_key, ts2_key, filter_name):
    if ts1_key in frame and ts2_key in frame:
        bolos = frame[ts1_key].keys()
        for b in bolos:
            ts1 = frame[ts1_key][b]
            ts2 = frame[ts2_key][b]
            if not np.argmax(ts1) == np.argmax(ts2):
                raise RuntimeError("Filter %s has a lag of %d samples"%(filter_name, np.argmax(ts1) - np.argmax(ts2)))

def compare_integrated_power(frame, ts1_key, ts2_key, rel_err, filter_name):
    if ts1_key in frame and ts2_key in frame:
        bolos = frame[ts1_key].keys()
        for b in bolos:
            ts1 = frame[ts1_key][b]
            ts2 = frame[ts2_key][b]
            if not np.isclose(np.sum(ts1), np.sum(ts2), rtol = rel_err):
                raise RuntimeError("Filter %s fails integrated power comparison at the %.2E level: result %.2E"%(filter_name, rel_err, np.abs(np.sum(ts1) - np.sum(ts2))/np.sum(ts1)))

def compare_peak_amplitude(frame, ts1_key, ts2_key, rel_err, filter_name):
    if ts1_key in frame and ts2_key in frame:
        bolos = frame[ts1_key].keys()
        for b in bolos:
            ts1 = frame[ts1_key][b]
            ts2 = frame[ts2_key][b]
            if not np.isclose(np.max(ts1), np.max(ts2), rtol = rel_err):
                    raise RuntimeError("Filter %s fails peak amplitude comparison at the %.2E level: result %.2E"%(filter_name, rel_err, np.abs(np.max(ts1) - np.max(ts2))/np.sum(ts1)))

######################################################################
# The test:
#
# * defines some time constants for simulated bolometers
# * simulates the exponential decay effect of a time constant on a unit impulse signal
# * does the deconvolution with the inverse kernel
# * checks that there is no phase lag or leakage in any of the filters
#
######################################################################

if __name__ == '__main__':

    bolos = ['a','b','c','d']
    sample_rate = 134.*core.G3Units.Hz
    taus = {'a': 0.5*core.G3Units.ms, 'b': 3*core.G3Units.ms, 'c': 10*core.G3Units.ms, 'd': 30*core.G3Units.ms}
    exp_err = 1e-10
    rel_err = 1e-5

    # begin and end with an empty frame to ignore edge effects in the convolution code
    # (these are tested for in the convolve_filter tester)

    delta_loc = [None,523,None]

    pipe = core.G3Pipeline()
    pipe.Add(DeltaFunctionFrames, nsamples = 1001, bolos = bolos, delta_location = delta_loc, sample_rate = sample_rate)
    pipe.Add(inject_taus, taus=taus)
    pipe.Add(ConvolveFilter, filter = exponential_filter(taus, sample_rate, rel_err = exp_err), keys = ['tods'], key_suffix='_Convolved')
    pipe.Add(tc.TimeConstantFilter, keys = ['tods_Convolved'])
    pipe.Add(lambda fr: 'Has Delta' in fr and fr['Has Delta'])
    pipe.Add(compare_peak_location, ts1_key = 'tods', ts2_key = 'tods_Convolved', filter_name = 'exponential_filter')
    pipe.Add(compare_peak_location, ts1_key =  'tods_Convolved', ts2_key = 'tods_Convolved_Deconvolved', filter_name = 'TimeConstantFilter')
    pipe.Add(compare_integrated_power, ts1_key = 'tods', ts2_key = 'tods_Convolved', rel_err = rel_err, filter_name = 'exponential_filter')
    pipe.Add(compare_integrated_power, ts1_key = 'tods_Convolved', ts2_key = 'tods_Convolved_Deconvolved', rel_err = rel_err, filter_name = 'TimeConstantFilter')
    pipe.Add(compare_peak_amplitude, ts1_key = 'tods', ts2_key = 'tods_Convolved_Deconvolved', rel_err = rel_err, filter_name = 'TimeConstantFilter after exponential_filter')
    pipe.Run()



