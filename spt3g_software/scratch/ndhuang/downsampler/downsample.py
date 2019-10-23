import os
from collections import deque
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pl
from spt3g import core, gcp, std_processing
from spt3g.todfilter import downsampler

# def convolve_renorm_ends(a, filt):
#     if len(filt) > len(a):
#         raise ValueError("`filt` should be longer than `a`")
#     out = np.empty_like(a)
#     offset = len(filt) / 2 + (len(filt) % 2)
#     out = np.convolve(a, filt, mode = 'same')
#     for i in xrange(offset):
#         if len(filt) % 2 == 1:
#             renorm = np.sum(filt[:offset + i])
#         else:
#             renorm = np.sum(filt[:offset + i])
#         out[i] /= renorm
#         renorm = np.sum(filt[offset - i - 1:])
#         out[-(i + 1)] /= renorm
#     return out
        
        
# def timedomain_sinc(nsamples, cutoff, window = 'hanning'):
#     '''
#     `cutoff` is a fraction of the Nyquist frequency
#     '''
#     t = np.arange(nsamples) - ((nsamples  - 1) / 2.)
#     t /= 2 # set Nyquist f at 1.0
#     filt = np.ones(nsamples)
#     inds = np.where(t != 0)[0]
#     arg = 2 * np.pi * t[inds] * cutoff
#     filt[inds] = np.sin(arg) / arg
#     if window:
#         filt *= signal.get_window(window, len(filt))

#     # fft method
#     # filt_f = np.ones(nsamples)
#     # freq = np.fft.fftfreq(nsamples, .5) # puts Nyquist f at 1.0
#     # inds = np.where(freq > cutoff)
#     # filt_f[inds] = 0.
#     # filt = np.fft.ifftshift(np.fft.ifft(filt_f))
#     # filt = abs(filt)

#     filt /= np.sum(filt)
#     return filt

def truncate_filter(filt, nbefore, nafter = None):
    if nafter is None:
        nafter = nbefore
    fout = filt[nbefore:-nafter]
    fout /= np.sum(fout)
    return fout
    
def quality_check(filt, cutoff, plot = True, newfig = False, verbose = True):
    freq_all = np.fft.fftfreq(len(filt), .5)
    freq = freq_all[np.where(freq_all >= 0)]
    ft = np.fft.fft(filt, norm = 'ortho')
    psd = 2 * np.abs(ft)**2
    psd = psd[np.where(freq_all >= 0)]
    tot_pow = np.trapz(psd)
    pass_inds = np.where(freq <= cutoff)
    pass_pow = np.trapz(psd[pass_inds])
    alias_pow = np.trapz(psd[np.where(freq >= cutoff)])
    
    pass_perfect = np.trapz(abs(filt)**2) / len(filt)
    pass_mean = np.mean(psd[pass_inds])
    pass_std = np.std(psd[pass_inds])
    
    half_pow = min(freq[np.where(psd <= .5 * pass_perfect)])
    
    if newfig and plot:
        pl.figure()
    if plot:
        pl.plot(freq, psd * len(psd), marker = 'o')
    
    
    if verbose:
        print "Power in passband: {:.03e}".format(pass_pow / tot_pow)
        print "Aliased power: {:.03e}".format(alias_pow / tot_pow)
        print "-3 dB: {:.1f} Hz".format(half_pow)
    
    out = dict(pass_pow = pass_pow, alias_pow = alias_pow, tot_pow = tot_pow,
               half_pow = half_pow, psd = psd, freq = freq)
    return out


def test_length(n_start, n_cut, savefig = False):
    def _single(ns, nc, cutoff):
        filt = timedomain_sinc(ns, cutoff)
        if nc > 0:
            filt = truncate_filter(filt, nc)
        out = quality_check(filt, cutoff, verbose = False)
        pass_frac = out['pass_pow'] / out['tot_pow']
        txt = "{:d}, pass: {:.3e}".format(len(filt), pass_frac)
        pl.plot(out['freq'], out['psd'] * len(out['psd']), marker = 'o', label = txt)
        cutoffs = 1. / np.array([2, 3, 4])
    for c in cutoffs:
        pl.figure()
        if np.iterable(n_start):
            for ns in n_start:
                for nc in n_cut:
                    _single(ns, nc, c)
        else:
            for nc in n_cut:
                _single(n_start, nc, c)
        ylim = pl.ylim()
        pl.vlines(c, *ylim)
        pl.ylim(ylim)
        pl.legend(loc = 'Best', fontsize = 'small')
        if savefig:
            pl.savefig('/home/ndhuang/plots/downsample/cut_{:.2f}-offset_0.png'.format(c))
            pl.close()
            

def test_equal_length(n, cut_frac, cutoff, savefig = False, window = 'hanning'):
    def _single(ns, nc, cutoff, window):
        filt = timedomain_sinc(ns, cutoff, window = window)
        if nc > 0:
            filt = truncate_filter(filt, nc)
        out = quality_check(filt, cutoff, verbose = False, plot = False)
        pass_frac = out['pass_pow'] / out['tot_pow']
        txt = "{:.2f}, pass: {:.3e}".format(np.round(float(nc) / ns * 2, 2), pass_frac)
        if window:
            txt += ', win'
        pl.plot(out['freq'], out['psd'] * len(out['psd']), marker = 'o', label = txt)
    cut = (np.round(n / (1 - np.array(cut_frac))).astype(int) - n) / 2
    pl.figure()
    for c in cut:
        _single(n + 2 * c, c, cutoff, window)
        _single(n + 2 * c, c, cutoff, False)
    ylim = pl.ylim()
    pl.vlines(cutoff, *ylim)
    pl.ylim(ylim)
    pl.xlim(0, 1)
    pl.legend(loc = 'best', fontsize = 'small')
    pl.xlabel('Freq (fraction of Nyquist)')
    pl.ylabel('Transmitted Power')
    if savefig:
        pl.savefig('/home/ndhuang/plots/downsample/cut_{:.2f}-offset_0.png'.format(c))
        pl.close()

def fancy_downsample(fr, ds_factor = 2, debug = False,
                     key_suffix = '_downsampled',
                     keys = None, decimate_keys = None):
    if fr.type != core.G3FrameType.Scan:
        return
    filter = timedomain_sinc(32, 1. / ds_factor)
    for k in keys:
        try:
            outmap = fr[k].copy()
            if debug:
                outmap_debug = outmap.copy()
        except AttributeError:
            outmap = core.G3Timestream(fr[k])
            outmap_debg = fr[k].__new__(fr[k])
        for subkey in fr[k].keys():
            if debug:
                outmap_debug[subkey] = convolve_renorm_ends(fr[k][subkey], filter)
                outmap[subkey] = outmap_debug[subkey][::ds_factor]
            else:
                outmap[subkey] = convolve_renorm_ends(fr[k][subkey], filter)[::ds_factor]
        new_k = k + key_suffix
        fr[new_k] = outmap
        if debug:
            fr[new_k + '_debug'] = outmap_debug
    return fr

# class Downsampler(object):
#     def __init__(self, ds_factor, debug = False, key_suffix = '_downsampled', 
#                  temp_suffix = '_DS_TEMP' keys = None, decimate_keys = None):
#         self.ds_factor = int(ds_factor)
#         self.filter = timedomain_sinc(32, 1. / ds_factor)
#         if isinstance(keys, str):
#             keys = [keys]
#         self.keys = keys
#         if decimate_keys is None:
#             decimate_keys = []
#         if isinstance(decimate_keys, str):
#             decimate_keys = [decimate_keys]
#         self.decimate_keys = decimate_keys
#         self.key_suffix = key_suffix
#         self.temp_suffix = temp_suffix
#         self.last_fr = None
#         self.debug = debug

#     def _first_scan_frame(self, fr):
#         for k in self.keys:
#             new_k = k + self.key_suffix
#             outmap = fr[k].copy()
#             for subkey in fr[k].keys():
#                 outmap[subkey] = np.round(convolve_renorm_ends(fr[k][subkey], 
#                                                                self.filter))
#                 # fr[new_k][subkey] = np.convolve(fr[k][subkey], 
#                 #                                 self.filter,
#                 #                                 mode = 'same')
#             fr[new_k + self.temp_suffix] = outmap
#         self.last_fr = fr
#         return []

#     def _intermediate_frames(self, fr):
#         # The recipe here is:
#         # 1. grab enough data from this frame to filter the last frame
#         # 2. filter the last frame
#         # 3. replace the last len(filter) / 2 samples of the data with properly
#         #    filtered data.
#         # 4. pack the last frame up with downsampled data
#         # 5. decimate any additional keys

#         half1 = int(np.floor((len(self.filter) - 1) / 2.))
#         half2 = int(np.ceil((len(self.filter) - 1) / 2.))
#         for k in self.keys:
#             new_k = k + self.key_suffix
#             outmap = core.G3TimestreamMap()
#             outmap_last = self.last_fr.pop(new_k + self.temp_suffix)
#             for subkey in outmap.keys():
#                 # Steps 1-3
#                 stitched = np.concatenate((self.last_fr[k][subkey][-len(self.filter):],
#                                            fr[k][subkey]))
#                 filtered = np.round(np.convolve(stitched,
#                                                 self.filter,
#                                                 mode = 'same'))
#                 outmap_last[subkey][-half2:] = filtered[half2:len(self.filter)]
#                 outmap[subkey] = filtered[len(self.filter):]
#             # Step 4
#             outmap_ds = core.G3TimestreamMap()
#             for subkey, value in outmap_last.iteritems():
#                 outmap_ds[subkey] = value[::self.ds_factor]
#             self.last_fr[new_k] = outmap_ds
#             fr[new_k + self.temp_suffix] = outmap
#             if debug:
#                 self.last_fr[new_k + '_debug'] = outmap_last
#         for k in self.decimate_keys:
#             # Step 5
#             new_k = k + self.key_suffix
#             self.last_fr[new_k] = core.G3Timestream(self.last_fr[k][::self.ds_factor])
#             self.last_fr[new_k].start = self.last_fr[k].start
#             self.last_fr[new_k].stop = self.last_fr[k].stop
#         ret_frame = self.last_fr
#         self.last_fr = fr
#         return ret_frame
        
#     def _last_frame(self, fr):
#         for k in self.keys:
#             new_k = k + self.key_suffix
#             outmap = core.G3TimestreamMap()
#             outmap_ds = core.G3TimestreamMap()
#             for subkey in self.last_fr[k].keys():
#                 thisdata = self.last_fr[new_k + self.temp_suffix][subkey]
#                 outmap[subkey][-len(self.filter):] = \
#                     np.round(convolve_renorm_ends(this_data,
#                                                   self.filter)[-len(self.filter):])
#                 outmap_ds[subkey] = outmap[subkey][::self.ds_factor]
#             self.last_fr[new_k] = outmap_ds
#             if debug:
#                 self.last_fr[new_k + '_debug'] = outmap
#             self.last_fr.pop(new_k + self.temp_suffix)
#         for k in self.decimate_keys:
#             new_k = k + self.key_suffix
#             self.last_fr[new_k] = core.G3Timestream(self.last_fr[k][::self.ds_factor])
#             self.last_fr[new_k].start = self.last_fr[k].start
#             self.last_fr[new_k].stop = self.last_fr[k].stop
#         return [self.last_fr, fr]

#     def __call__(self, fr):
#         if fr.type == core.G3FrameType.EndProcessing:
#             return self._last_frame(fr)
#         if fr.type != core.G3FrameType.Scan:
#             return
#         if self.last_fr is None:
#             return self._first_scan_frame(fr)
#         else:
#             return self._intermediate_frames(fr)

def make_plots(length, savedir, saveprefix = 'downsampler_'):
    for ds in [2., 3., 4.]:
        test_equal_length(length, [0, .5, .25], 1./ds, window = 'hanning')
        pl.title('Transfer Func, {:d}x'.format(int(ds)))
        fname = '{}{:d}equal_{:d}x.png'.format(saveprefix,
                                               length,
                                               int(ds))
        pl.savefig(os.path.join(savedir, fname))
    

class NoiseFrames(object):
    def __init__(self):
        self.i = 0

    def __call__(self, fr):
        assert(fr is None)
        fr = core.G3Frame(core.G3FrameType.Scan)
        ts = core.G3Timestream(np.random.randn(1024))
        tm = core.G3TimestreamMap()
        tm['noise'] = ts
        fr['noise'] = tm
        if self.i % 2 == 0:
            fr['Turnaround'] = True
        self.i += 1
        if self.i >= 5:
            return []
        return fr

class LinearFrames(object):
    def __init__(self):
        self.i = 0
    
    def __call__(self, fr):
        assert(fr is None)
        fr = core.G3Frame(core.G3FrameType.Scan)
        ts = core.G3Timestream(np.arange(self.i * 1024, (self.i + 1) * 1024))
        tm = core.G3TimestreamMap()
        tm['noise'] = ts
        fr['noise'] = tm
        if self.i % 2 == 0:
            fr['Turnaround'] = True
        self.i += 1
        if self.i >= 5:
            return []
        return fr
        

class Accumulator():
    def __init__(self):
        self.d = deque()
        self.breaks = -1
    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        self.d.append(fr['noise_downsampled']['noise'])
        return


def test(linear = False):
    accum = Accumulator()
    if linear:
        frames = LinearFrames()
    else:
        frames = NoiseFrames()
    ds = Downsampler(2, 'noise')
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(core.Dump)
    pipe.Add(ds)
    pipe.Add(accum)
    pipe.Run()
    out = np.array([np.array(ts) for ts in accum.d])
    return out.flatten()

if __name__ == '__main__':
    # A happy little tester for the official downsampler
    frames = LinearFrames()
    ds = downsampler.Downsampler(2, 'noise', compress = False)
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(ds)
    # pipe.Add(core.Dump)
    pipe.Run()
