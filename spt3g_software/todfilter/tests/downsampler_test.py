#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.todfilter import downsampler as ds

class TestFrames(object):
    '''
    Creates frames to test the downsampler.  To test the transfer function,
    we want something that is flat in frequency space (at least up to 
    Nyquist), so this uses the sinc function from the downsampler.
    '''
    def __init__(self, nsamples):
        self.n_frames = 0
        self.nsamples = nsamples
        
    def __call__(self, fr):
        assert(fr is None)
        if self.n_frames >= 5:
            return []
        start = core.G3Time.Now()
        nsamples = self.nsamples + self.n_frames % 5
        ts = core.G3Timestream(ds.timedomain_sinc(nsamples, 1, False)) * 1e5
        ts.start = start
        stop = core.G3Time.Now()
        ts.stop = stop
        tm = core.G3TimestreamMap()
        tm['ts'] = ts
        fr = core.G3Frame(core.G3FrameType.Scan)
        fr['tods'] = tm
        dec_ts = core.G3Timestream(np.arange(nsamples, dtype='float64'))
        dec_ts.start = start
        dec_ts.stop = stop
        fr['dec'] = dec_ts
        self.n_frames += 1
        return fr

class SinusoidFrames(object):
    '''
    Creates frames to test the downsampler.  These are designed to test whether
    the downsampler is maintaining an even samplerate at frame boundaries.
    '''
    
    def __init__(self, nsamples = 100, nframes = 10, ds_factor = 2):
        self.nsamples = nsamples
        self.ds_factor = ds_factor
        self.maxframe = nframes
        self.nframes = 0
        self.last = 0
        self.last_x = 0
        self.length_offsets = [0, 1, 1, 0]
        self.start = core.G3Time.Now()

    def __call__(self, fr):
        assert(fr is None)
        if self.nframes >= self.maxframe:
            return []
        # test both odd and even length frames, but not always back to back
        ts_length = self.nsamples + self.length_offsets[self.nframes % 4]
        x = np.arange(self.last_x, self.last_x + ts_length)
        self.last_x = x[-1] + 1
        ts = np.sin(x * np.pi / self.ds_factor / 2)
        ts = core.G3Timestream(ts)
        # ts.start = self.start
        # ts.stop = self.start + ts_length 
        # self.start = ts.stop + 1
        tm = core.G3TimestreamMap()
        tm['ts'] = ts
        ts_decimate = core.G3Timestream(np.zeros(ts_length))
        fr = core.G3Frame(core.G3FrameType.Scan)
        fr['tm'] = tm
        fr['length_canary'] = ts_decimate
        self.nframes += 1
        return fr

class TestFramesWithUnits(object):
    '''
    If timestream units are in counts, they must be integers.  Make sure
    this is preserved
    '''
    def __init__(self, nsamples = 100):
        self.n_frames = 0
        self.nsamples = nsamples
        
    def __call__(self, fr):
        assert(fr is None)
        if self.n_frames >= 1:
            return []
        ts = core.G3Timestream(np.random.randint(10, size = self.nsamples).astype(float),
                               units = core.G3TimestreamUnits.Counts)
        tm = core.G3TimestreamMap()
        tm['ts'] = ts
        fr = core.G3Frame(core.G3FrameType.Scan)
        fr['tods'] = tm
        self.n_frames += 1
        return fr

class DeltaFunctionFrames(object):
    '''
    Creates frames to test the downsampler.  This class produces timestreams 
    with a delta function, to ensure that we have constructed a zero-lag filter.
    '''
    def __init__(self, nsamples = 1001, delta_location = [23]):
        self.nsamples = nsamples
        self.delta_loc = delta_location
        self.maxframe = len(delta_location)
        self.nframes = 0

    def __call__(self, fr):
        assert(fr is None)
        if self.nframes >= self.maxframe:
            return []
        ts = np.zeros(self.nsamples)
        ts[self.delta_loc[self.nframes]] = 1024
        ts = core.G3Timestream(ts)
        tm = core.G3TimestreamMap()
        tm['ts'] = ts
        fr = core.G3Frame(core.G3FrameType.Scan)
        fr['tm'] = tm
        self.nframes += 1
        return fr

class TestFrameBoundarySampling(object):
    '''
    Ensure that we're not injecting or skipping samples at frame boundaries.
    If everything is setup right, downsampled timestreams should give
    [-1., 0., 1., 0.] (i.e. a sin wave at Nyquist) in sequence, after the edge 
    effects are gone.
    '''
    def __init__(self):
        self.nframes = 0
        self.last = None
        self.template = np.array([-1., 0., 1., 0.])

    def __call__(self, fr):
        if fr.type != core.G3FrameType.Scan:
            return 
        if self.last is None:
            self.last = fr
            return fr
        for ds in (2,):
            tm_key = 'tm_{:d}x'.format(ds)
            dec_key = 'length_canary_{:d}x'.format(ds)
            joint = np.concatenate([self.last[tm_key]['ts'][-2:],
                                    fr[tm_key]['ts'][:2]])
            good = False
            for i in range(4):
                if all(self.template == np.roll(joint, i)):
                    good = True
                    break
            if not good:
                raise RuntimeError('Frames are not joining up correctly!')
            if len(fr[dec_key]) != len(fr[tm_key]['ts']):
                raise RuntimeError('Decimated and downsampled keys have different length!')
        self.last = fr
        return fr

def test_ds_transferfn(fr):
    '''
    Check the transfer function.  Require that at least 95% of the
    power be below the new Nyquist.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    for dsfactor in (2, 3):
        new_key = 'tods' + str(dsfactor) + 'x_debug'
        ft = np.fft.fft(fr[new_key]['ts'])
        psd = np.abs(ft[:len(ft) // 2])**2
        if sum(psd[:len(psd) // dsfactor]) / sum(psd) < .95:
            raise RuntimeError('Too much aliasing')
            
def test_ds_decimate(fr):
    '''
    Check the keys that are only decimated (without filtering)
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    for dsfactor in (2, 3):
        new_key = 'tods' + str(dsfactor) + 'x'
        # print(len(fr[new_key]['ts']), len(fr['dec' + str(dsfactor) + 'x']))
        if len(fr[new_key]['ts']) != len(fr['dec' + str(dsfactor) + 'x']):
            raise RuntimeError('Decimated timestream is not the same length as downsampled timestream.')

def test_ds_metadata(fr):
    '''
    Check the metadata that should be associated with the timestreams.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    start = fr['tods']['ts'].start
    stop = fr['tods']['ts'].stop
    for dsfactor in (2, 3):
        ds_key = 'tods' + str(dsfactor) + 'x'
        dec_key = 'dec' + str(dsfactor) + 'x'
        dt = 1. / fr[ds_key].sample_rate
        if abs(fr[ds_key]['ts'].start.time - start.time) > dt:
            raise RuntimeError('Start time is wrong')
        if abs(fr[dec_key].start.time - start.time) > dt:
            raise RuntimeError('Decimated start time is wrong')
        if abs(fr[ds_key]['ts'].stop.time - stop.time) > dt:
            raise RuntimeError('Stop time is wrong')
        if abs(fr[dec_key].stop.time - stop.time) > dt:
            raise RuntimeError('Decimated stop time is wrong')

def test_ds_length(fr):
    '''
    Check that the resulting time streams are the appropriate length.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    for dsfactor in (2, 3):
        new_key = 'tods' + str(dsfactor) + 'x'
        expected_len = np.floor(len(fr['tods']['ts']) / dsfactor + 1)
        if abs(len(fr[new_key]['ts']) - expected_len) > 1:
            raise RuntimeError('Timestream wrong length')

class TestLag(object):
    '''
    Check that the filtering has zero lag
    '''
    def __init__(self, delta_location = []):
        self.nframes = 0
        self.delta_loc = delta_location

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            assert(self.nframes == len(delta_location))
        if fr.type != core.G3FrameType.Scan:
            return
        this_delta = self.delta_loc[self.nframes]
        errmsf = 'x ds has lag with delta location: {:d}'.format(this_delta)
        if np.argmax(fr['tm_2x_debug']['ts']) != this_delta:
            raise RuntimeError('Lag in 2x ds with delta location = {:d}'.format(this_delta))
        if np.argmax(fr['tm_3x_debug']['ts']) != this_delta:
            raise RuntimeError('Lag in 3x ds with delta location = {:d}'.format(this_delta))
        self.nframes += 1

def test_counts(fr):
    '''
    If units are in `core.G3TimestreamUnits.Counts`, they must come out as 
    whole numbers.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    tsvals = np.array(fr['tods']['ts'])
    if not all(tsvals == tsvals.astype(int)):
        raise RuntimeError('Counts went in and floats came out!')
        
if __name__ == '__main__':
    frames = TestFrames(1001)
    pipe = core.G3Pipeline()
    pipe.Add(frames)
    pipe.Add(ds.Downsampler, ds_factor = 2, keys = 'tods', decimate_keys = 'dec',
             key_suffix = '2x', debug = True, compress = False)
    pipe.Add(ds.Downsampler, ds_factor = 3, keys = 'tods', decimate_keys = 'dec',
             key_suffix = '3x', debug = True, compress = False)
    pipe.Add(test_ds_transferfn)
    pipe.Add(test_ds_decimate)
    pipe.Add(test_ds_length)
    pipe.Add(test_ds_metadata)
    pipe.Run()
    

    pipe = core.G3Pipeline()
    pipe.Add(SinusoidFrames)
    pipe.Add(ds.Downsampler, ds_factor = 2, keys = 'tm', key_suffix = '_2x', compress = False,
             decimate_keys = 'length_canary')
    pipe.Add(ds.Downsampler, ds_factor = 3, keys = 'tm', key_suffix = '_3x', compress = False,
             decimate_keys = 'length_canary')
    pipe.Add(TestFrameBoundarySampling)
    pipe.Run()

    pipe = core.G3Pipeline()
    delta_location = [323, 0, 1, 1000, 350, 336]
    pipe.Add(DeltaFunctionFrames, delta_location = delta_location)
    pipe.Add(ds.Downsampler, ds_factor = 2, keys = 'tm', decimate_keys = 'dec',
             key_suffix = '_2x', compress = False, debug = True)
    pipe.Add(ds.Downsampler, ds_factor = 3, keys = 'tm', decimate_keys = 'dec',
             key_suffix = '_3x', compress = False, debug = True)
    pipe.Add(TestLag, delta_location = delta_location)
    pipe.Run()

    pipe = core.G3Pipeline()
    pipe.Add(TestFramesWithUnits, nsamples = 518)
    pipe.Add(ds.Downsampler, ds_factor = 2, keys = 'tods', decimate_keys = 'dec',
             key_suffix = '_2x', compress = False)
    pipe.Add(test_counts)
    pipe.Run()
