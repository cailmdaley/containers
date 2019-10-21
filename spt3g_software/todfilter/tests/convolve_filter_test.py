import numpy as np
from spt3g import core
from spt3g.todfilter.convolve_filter import ConvolveFilter

class MultipleTimestreamTestFrames(object):
    '''
    Produce frames with multiples timestreams in a Map to test the per-bolometer 
    filtering.
    '''
    def __init__(self, nsamples = 625, nts = 2):
        self.nframes = 0
        self.nsamples = nsamples
        self.nts = nts
    
    def __call__(self, fr):
        assert fr is None
        if self.nframes >= 3:
            return []
        tm = core.G3TimestreamMap()
        # put something in the timestreams...
        ts = core.G3Timestream(.01 * np.arange(self.nsamples)**2)
        for i in range(self.nts):
            tm[str(i)] = core.G3Timestream(ts)
        fr = core.G3Frame(core.G3FrameType.Scan)
        fr['tods'] = tm
        self.nframes += 1
        return fr

def test_multiple_filters(fr):
    '''
    If we provided different filters, just make sure that they create
    different outputs.
    '''
    if fr.type != core.G3FrameType.Scan:
        return
    ts0 = fr['tods_filt']['0']
    ts1 = fr['tods_filt']['1']
    if all(ts0 == ts1):
        raise RuntimeError('Per-bolometer filtering produces identical output!')

if __name__ == '__main__':
    # Test that multiple filters at least do different things
    kern0 = np.zeros(5)
    kern0[2] = 1
    kern1 = np.ones(5) / 5.
    filters = {'0': kern0, '1': kern1}
    pipe = core.G3Pipeline()
    pipe.Add(MultipleTimestreamTestFrames)
    pipe.Add(ConvolveFilter, filter = filters, keys = 'tods', key_suffix = '_filt')
    pipe.Run()
