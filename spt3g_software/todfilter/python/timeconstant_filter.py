import numpy as np
from spt3g import core
from .convolve_filter import ConvolveFilter

#####################################################################################
# This class grabs timeconstants from the tau_key in the Calibration frame,
# constructs a real-space deconvolution kernel for each bolo and sends it to 
# ConvolveFilter. Default behavior for missing timeconstant values is to convolve
# with a Delta function (i.e., do nothing).
#
# The convolution kernel of a timeconstant is modeled as np.exp(-t/tau) for t >= 0.
# Going to frequency domain, inverting, and going back - i.e. taking 
#
#   inverse Fourier Transform ( 1. / Fourier Transform(kernel) )
#
# gives the corresponding real-space inverse kernel and can be solved analytically as
#
#   Delta[0] - tau*Delta'[0]
#
# which can be discretized as the normalized 3-sample kernel [0 1 -np.exp(-dt/tau)] 
# where dt is the sample spacing
#####################################################################################

class TimeConstantFilter(ConvolveFilter):
    def __init__(self, keys=['RawTimestreams_I'], key_suffix='_Deconvolved', tau_key = 'TimeConstants'):
        super(TimeConstantFilter, self).__init__(filter=[1.], keys=keys, key_suffix=key_suffix)
        self.tau_key = tau_key
        self.taus = None

    def __call__(self, frame):
        if self.tau_key in frame:
            self.taus = frame[self.tau_key]
        if frame.type == core.G3FrameType.Scan:
            tsunits = np.array([frame[tskey].values()[0].units for tskey in self.keys])
            if (tsunits == core.G3TimestreamUnits.Counts).any():
                core.log_fatal("Using time constant deconvolution with timestreams in Counts will produce rounding errors. Convert to another unit first.")
            if self.taus is None:
                core.log_fatal("Could not find time constant key %s in this frame or any preceding frame"%self.tau_key)
            else:
                filter = {}
                dt = 1./frame[self.keys[0]].sample_rate
                for bolo in frame[self.keys[0]].keys():
                    kernel = np.zeros(3)
                    kernel[1] = 1.0
                    if bolo in self.taus:
                        kernel[2] = -np.exp(-dt/self.taus[bolo])
                    filter[bolo] = kernel/kernel.sum()
                self.filter = filter
                self.lenfilt = len(list(filter.values())[0])
                self.filter_per_bolo = True
        return super(TimeConstantFilter, self).__call__(frame)
