import numpy as np
import math
from scipy import signal
from spt3g import core

@core.indexmod
class ConvolveFilter(object):
    '''
    A module to convolve time-domain filters with bolometer data.  This module
    stitches together data from adjacent scan frames to ensure that the 
    convolution happens on continuous data.  Hence, this module should always
    appear before dropping turnarounds.
    The actual behavior of the module is controlled by initialization 
    parameters.  See __init__ documentation for details.
    

    Potentially complicated features one could add:
    1. Windowing.  Since we are using entire timestreams, this would require two
        passes through the data, which makes life rather difficult.
        One could imagine implementing a window function that only apodizes the
        a few samples at the beginning and end, which would be much easier, but
        less ideal.
    2. Support for filter kernels that are longer than the shortest two adjacent
        scans.  To be fully general, this would require choosing the number of
        scan frames to load dynamically, and handle repackaging the filtered
        timestreams.
    '''
    def __init__(self, filter = None, keys=['RawTimestreams_I'],
                 key_suffix = '_filt', temp_suffix = '_FILT_TEMP', 
                 debug = False):
        '''
        filter: array-like or dict-like
            The filter kernel.  If `filter` has a `keys` method, 
            bolometer `bolo` will be convolved with `filter[bolo]`.  In this
            all kernels must have the same length.
            NB: If you only want a few bolometers to use a different filter,
            you can do something clever with `collections.default_dict`.
        keys: list, str
            A list of keys to filter.  Each key must be a valid key
            in every ScanFrame.  This is usually bolomter data.
        key_suffix: str, '_filt'
            This string is added to the end of each key in `keys` and
            and the filtered data is stored in the new key.
        temp_suffix: str, '_FILT_TEMP'
            Internal temporary keys will be stored with this suffix.  
            You probably don't need to worry about this.
        debug: bool, False
            Include debugging output.  You probably don't want this.
        '''
        assert(filter is not None)
        self.filter = filter
        self.filter_per_bolo = hasattr(self.filter, 'keys')
        if self.filter_per_bolo:
            self.lenfilt = len(list(self.filter.values())[0])
            for k in self.filter.keys():
                if len(self.filter[k]) != self.lenfilt:
                    raise ValueError('Not all filters have the same length')
        else:
            self.lenfilt = len(self.filter)

        if isinstance(keys, str):
            keys = [keys]
        self.keys = keys
        self.key_suffix = key_suffix
        self.temp_suffix = temp_suffix
        self.i_start = 0
        self.nframes = 0
        self.last_fr = None
        self.debug = debug

    def _convolve_renorm_ends(self, kernel, a):
        # This is convolution that treats edges specially.
        # If the kernel extends beyond the end of the data,
        # those parts of the kernel are zeroed, and the entire kernel
        # is renormalized.  Then convolution continues as usual.
        if self.lenfilt > len(a):
            raise ValueError("`filt` should not be longer than `a`")
        out = np.empty_like(a)
        offset = self.lenfilt // 2 + (self.lenfilt % 2)
        out[:] = np.convolve(a, kernel, mode = 'same')
        for i in range(offset):
            out[i] = np.sum(a[:offset + i] * kernel[:offset + i])
            out[i] /= np.sum(kernel[:offset + i])
            out[-(i + 1)] = np.sum(a[-(offset + i):] * 
                                   kernel[-(offset + i):])
            out[-(i + 1)] /= np.sum(kernel[-(offset + i):])
        return out

    def _first_scan_frame(self, fr):
        for k in self.keys:
            new_k = k + self.key_suffix
            outmap = core.G3TimestreamMap()
            for subkey in fr[k].keys():
                if self.filter_per_bolo:
                    kernel = self.filter[subkey]
                else:
                    kernel = self.filter
                conv = self._convolve_renorm_ends(kernel, fr[k][subkey])
                outmap[subkey] = core.G3Timestream(fr[k][subkey])
                if self.in_counts:
                    conv = np.round(conv)
                outmap[subkey][:] = conv
            fr[new_k + self.temp_suffix] = outmap
        self.last_fr = fr
        return []

    def _intermediate_frames(self, fr):
        # The recipe here is:
        # 1. grab enough data from this frame to filter the last frame
        # 2. filter the last frame
        # 3. replace the last len(filter) / 2 samples of the data with properly
        #    filtered data.
        half1 = int(np.floor((self.lenfilt - 1) / 2.))
        half2 = int(np.ceil((self.lenfilt - 1) / 2.))
        for k in self.keys:
            new_k = k + self.key_suffix
            outmap = core.G3TimestreamMap()
            outmap_last = self.last_fr.pop(new_k + self.temp_suffix, None)
            for subkey in fr[k].keys():
                # Steps 1-3
                stitched = np.concatenate((self.last_fr[k][subkey][-self.lenfilt:],
                                           fr[k][subkey]))
                if self.filter_per_bolo:
                    kernel = self.filter[subkey]
                else:
                    kernel = self.filter
                filtered = np.convolve(stitched, kernel, mode = 'same')
                if self.in_counts:
                    filtered = np.round(filtered)
                outmap_last[subkey][-half2:] = filtered[half2 + self.lenfilt % 2:self.lenfilt]
                outmap[subkey] = core.G3Timestream(fr[k][subkey])
                outmap[subkey][:] = filtered[self.lenfilt:]
            self.last_fr[new_k] = outmap_last
            fr[new_k + self.temp_suffix] = outmap
        ret_frame = self.last_fr
        self.last_fr = fr
        return ret_frame
        
    def _last_frame(self, fr):
        for k in self.keys:
            new_k = k + self.key_suffix
            outmap_last = self.last_fr.pop(new_k + self.temp_suffix, None)
            for subkey in self.last_fr[k].keys():
                thisdata = outmap_last[subkey]
                if self.nframes > 1:
                    if self.filter_per_bolo:
                        kernel = self.filter[subkey]
                    else:
                        kernel = self.filter
                    conv = self._convolve_renorm_ends(kernel, thisdata)[-self.lenfilt:]
                    if self.in_counts:
                        conv = np.round(conv)
                    outmap_last[subkey][-self.lenfilt:] = conv
            self.last_fr[new_k] = outmap_last
        return [self.last_fr, fr]

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            if self.last_fr is None:
                return
            return self._last_frame(fr)
        self.nframes += 1
        if fr.type != core.G3FrameType.Scan:
            return
        # If the timestreams are in counts, we need to round to ints
        tsunits = fr[self.keys[0]].values()[0].units
        self.in_counts =  tsunits == core.G3TimestreamUnits.Counts
        if self.last_fr is None:
            return self._first_scan_frame(fr)
        else:
            return self._intermediate_frames(fr)

    def drop_short_scans(self, fr):
        '''
        Drop any scans shorter than the filter length.
        '''
        if fr.type != core.G3FrameType.Scan:
            return
        for k in self.keys:
            if fr[k].n_samples < self.lenfilt:
                return False
        return True


##########################################################################
# Filter Kernels

#######################################################################
#                Some notes on filtering kernels
#
# These are meant as guidelines.  If you want to really understand
# what you're doing, go read Bracewell (note that the author of
# this section has not read Bracewell).
#
# In order to avoid phase shifts from your filter, you want
#     1. An odd length kernel, and
#     2. The peak of the kernel to be in the center.
# The second condition _should_ be broken if it is appropriate to
# model the process your filter is representing.
#
# If you want your filter to have gain one (i.e. convserve power or 
# keep the same units), normalize it such that its sum is equal to 1.
#
# Determining the overall length of the filter requires some artistry.
# You need to decide, based on your application, what is acceptable.
# For the downsampler, this was set based on the amount of aliased
# power we thought was acceptable, and the flatness of the passband.
#######################################################################

def timedomain_sinc(nsamples, cutoff, window = True):
    '''
    Create a sinc function for time-domain convolution.
    This is used as an anti-aliasing filter for the downsampler.
    The sinc is `nsamples` long, and cuts frequencies above
    `cutoff` * Nyquist Frequency.
    INPUTS
    ------
    nsamples: (int)
        The number of samples in the filter.  Longer filters will
        have a sharper cutoff (at the cost of larger support).
    cutoff: (float, (0.0, 1.0)
        The cutoff frequency specified as a fraction of the Nyquist
        frequency.
    window: bool: True
        If True, apply a Hann window.
    '''
    if cutoff == 1.0:
        # One would never do this, but it's important for a test case
        return np.ones(nsamples) / nsamples
    if cutoff > 1.0 or cutoff <= 0.0:
        raise RuntimeError('Cutoff frequency must be in the range (0.0, 1.0]')
    t = np.arange(nsamples) - ((nsamples  - 1) / 2.)
    t /= 2 # set Nyquist f at 1.0
    filt = np.ones(nsamples)
    inds = np.where(t != 0)[0]
    arg = 2 * np.pi * t[inds] * cutoff
    filt[inds] = np.sin(arg) / arg
    if window:
        # Construct a Hann window that is symmetric for an even number of
        # samples.  scipy.signal's window functions are not symmetric.
        window = (1 - np.cos((2 * np.pi * np.arange(nsamples + 2)) / 
                             (nsamples + 1))) / 2
        window = window[1:-1]
    filt *= window
    if np.sum(filt) == 0:
        raise RuntimeError(str(filt))
    filt /= np.sum(filt)
    return filt
