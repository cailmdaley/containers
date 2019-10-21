#!/usr/bin/env python

from spt3g import core
from spt3g import todfilter

import unittest
import numpy as np

np.random.seed(4)
def make_noise_scan_frame(ts_len, noise, offset, sample_rate, scale_poly = 3):
    scan_frame = core.G3Frame(core.G3FrameType.Scan)
    ts = core.G3TimestreamMap()
    for d in range(3):
        ts['test_noise_%d' % d] = core.G3Timestream(np.random.normal(loc=offset, scale=noise,size=ts_len) + np.arange(ts_len)**scale_poly)

    scan_frame['SampleRate'] = float(sample_rate)
    scan_frame['TsMap'] = ts
    return scan_frame

def is_mostly_zero(arr, slop  = 10):
    return np.max(np.abs(arr)) < slop


def shitty_psd( a, sample_rate):
    return np.fft.fftfreq(len(a), 1.0/sample_rate)[:len(a)//2], np.abs(np.fft.fft(np.hanning(len(a)) * a)[:len(a)//2])


def test_poly():
    poly_order = 3
    masked_hpf_freq = -1
    sample_rate = 100
    sf = make_noise_scan_frame(ts_len = 1001, noise = 1.0, offset = 12.3, sample_rate = sample_rate, scale_poly = 3)
    mhpf_poly3 = todfilter.MaskedPolyHpf('TsMap', 'PolyMhpfTsMap3', 3, masked_hpf_freq, False, '')
    mhpf_poly2 = todfilter.MaskedPolyHpf('TsMap', 'PolyMhpfTsMap2', 2, masked_hpf_freq, False, '')
    mhpf_poly3(sf)
    mhpf_poly2(sf)
    
    '''
    import pylab as pl
    
    pl.plot(sf['TsMap']['test_noise_0'])
    pl.plot(sf['PolyMhpfTsMap2']['test_noise_0'])
    pl.show()
    '''
    for d in range(3):
        assert(not is_mostly_zero(sf['PolyMhpfTsMap2']['test_noise_%d'%d], slop=10))
        '''
        import pylab as pl
        pl.plot(sf['PolyMhpfTsMap3']['test_noise_%d'%d])
        pl.show()
        '''
        assert(is_mostly_zero(sf['PolyMhpfTsMap3']['test_noise_%d'%d], slop=10))
def test_lowpass():
    sample_rate = 100
    sf = make_noise_scan_frame(ts_len = 1000, noise = 1.0, offset = 12.3, sample_rate = sample_rate, scale_poly = 0)
    sf['padding'] = 2048
    
    todfilter.dftutils.LowPassFilterSpecifier(sf,
                                              cutoff_freq = 30.0,        
                                              sample_rate_override_key = 'SampleRate',
                                              ts_map_key = 'TsMap', 
                                              input_filter_field = None, 
                                              output_filter_field = 'LowPassFilter', 
                                              already_specified_key = 'youdontmaptter', 
                                              padding_key = 'padding')
    todfilter.dftutils.FftFilter(sf, 'TsMap', 'LowPassFilter', 'OutFFTTs')
    
    for ts_name in ['test_noise_%d' % d for d in range(3)]:
        freq, befpsd = shitty_psd( sf['TsMap'][ts_name], sample_rate)
        freq, aftpsd = shitty_psd( sf['OutFFTTs'][ts_name], sample_rate)
        expected = todfilter.dftutils.lowpass_func(freq, cutoff_freq = 30)

        if 0:
            import pylab as pl
            pl.plot(freq, aftpsd/befpsd)
            pl.plot(freq, expected)
            pl.show()

        assert(is_mostly_zero(aftpsd/befpsd - expected, slop = 0.05))

def pypoly(ts, deg):
    x = np.arange(len(ts))
    p = np.polyfit(x, ts, deg)
    fit_poly = np.polyval(p, x)
    return ts-fit_poly, fit_poly

def test_mhpf():
    poly_order = 10
    masked_hpf_freq = 20
    sample_rate = 100
    
    sf   = make_noise_scan_frame(ts_len = 1000, noise = 1.0, offset = 12.3, 
                                 sample_rate = sample_rate, scale_poly = 1)
    #mhpf = todfilter.MaskedPolyHpf('TsMap', 'PolyMhpfTsMap', 7, poly_order, False, '', 
    #                               sample_rate_override_key = 'SampleRate')
    

    poly = todfilter.MaskedPolyHpf('TsMap', 'inter', poly_order, -1,
                                   False, '', 
                                   sample_rate_override_key = 'SampleRate')
    mhpf = todfilter.MaskedPolyHpf('inter', 'PolyMhpfTsMap', 0, masked_hpf_freq, 
                                   False, '', 
                                   sample_rate_override_key = 'SampleRate')
    poly(sf)
    mhpf(sf)
    
    for ts_name in ['test_noise_%d' % d for d in range(3)]:
        freq, befpsd = shitty_psd( sf['TsMap'][ts_name], sample_rate)
        freq, aftpsd = shitty_psd( sf['PolyMhpfTsMap'][ts_name], sample_rate)    

        if 0:
            import pylab as pl
            pl.plot(sf['TsMap'][ts_name])
            pl.plot(sf['PolyMhpfTsMap'][ts_name])
            pl.plot(sf['TsMap'][ts_name]-sf['PolyMhpfTsMap'][ts_name])

            pyts, pyp = pypoly(sf['TsMap'][ts_name], poly_order)
            pl.plot(pyp)
            print(np.std(pyp - (sf['TsMap'][ts_name]-sf['PolyMhpfTsMap'][ts_name])))
            pl.show()

        if 0:
            import pylab as pl

            pl.plot(freq, aftpsd/befpsd)
            pyts, pyp = pypoly(sf['TsMap'][ts_name], poly_order)
            freq, aftpsdpol = shitty_psd( pyts, sample_rate)
            pl.plot(freq, aftpsdpol/befpsd)

            #pyts, pyp = pypoly(sf['TsMap'][ts_name], poly_order)
            freq, aftpsdpol = shitty_psd( sf['inter'][ts_name], sample_rate)
            pl.plot(freq, aftpsdpol/befpsd)
            pl.show()


        delta = (aftpsd/befpsd) - (freq > 20).astype(int)
        delta = delta[np.logical_or( freq < 19, freq > 21)]
        if 0:
            import pylab as pl
            pl.plot(delta)
            pl.show()
        assert(is_mostly_zero(delta, slop = 0.05))



test_mhpf()
test_poly()
test_lowpass()
