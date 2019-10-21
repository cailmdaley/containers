from spt3g import core
import numpy as np
from copy import copy

@core.indexmod
class CutTimestreamsWithoutProperties(object):
    '''
    Cut detectors in the timestream map <input>, writing the results to
    <output>, not contained in the given <map>. By default, removes
    detectors without a calibration in the bolometer properties map.
    '''
    def __init__(self, input='CalTimestreams', output='TrimmedTimestreams',
      map='BolometerProperties'):
        self.input = input
        self.output = output
        self.map = map
        self.stored_map = None

    def __call__(self, frame):
        if self.map in frame:
            self.stored_map = frame[self.map]
        if self.input not in frame:
            return
        out = core.G3TimestreamMap()
        for b,ts in frame[self.input].iteritems():
            if b in self.stored_map:
                out[b] = ts
        frame[self.output] = out

@core.usefulfunc
def mask_timestream(ts_in, masked_locs, mask_with = None):
    '''
    Masks point sources at the timestream level for the input timestream
    ts_in.  By default, masks with inpainting (draws a line from the sample
    directly before the masking to the sample directly after), though you
    can mask with zeros or nans or whatever you want by setting the argument
    mask_with.  Mask locations must be provided by masked_locs (0 for
    unmasked, 1 for masked).
    '''
    ts_masked = copy(ts_in)
    if mask_with is not None:
        ts_masked[np.asarray(masked_locs) == 1] = mask_with

    else:
        i = 0
        mask_size = 0
        for masked in masked_locs:

            if masked:
                mask_size += 1
                if i+1 == len(ts_masked): #final sample is masked
                    inpaint = ts_masked[i-mask_size] * np.ones(mask_size)
                    ts_masked[(i-mask_size)+1:] = inpaint

            else:
                if mask_size > 0:
                    if i == mask_size: #first sample is masked
                        inpaint = ts_masked[i] * np.ones(mask_size)
                        ts_masked[0:mask_size] = inpaint
                    else:
                        inpaint = np.linspace(ts_masked[i-(mask_size+1)],ts_masked[i],num=mask_size+2)
                        ts_masked[i-mask_size:i] = inpaint[1:-1]
                    mask_size = 0

            i += 1

    return ts_masked

@core.usefulfunc
def mask_ts_map(ts_map_in, pntsrc_map, mask_with = None):
    '''
    Mask the point sources for an entire group of timestreams.  By
    default, masks with inpainting (draws a line from the sample directly
    before the masking to the sample directly after), though you can
    mask with zeros or nans or whatever you want by setting the argument
    mask_with.  Mask locations must be provided by pntsrc_map.
    '''
    ts_map_out = core.G3TimestreamMap()
    for k in ts_map_in.keys():
        if k in pntsrc_map:
            ts_map_out[k] = mask_timestream(ts_map_in[k], pntsrc_map[k],
                                            mask_with = mask_with)
        else:
            core.log_fatal('You have not provided the samples to be '
                           'masked for bolo %s'%k)
    return ts_map_out



