import numpy as np
from spt3g import core
from .convolve_filter import timedomain_sinc, ConvolveFilter

class LPFAndDecimate(object):
    def __init__(self, keys = ['RawTimestreams_I'], ds_factor = 2, 
                 decimate_keys = None, compress = True,
                 key_suffix = '_downsampled', temp_suffix = '_DS_TEMP',
                 debug = False):
        '''
        keys: list, str
            A list of keys to downsample.  Each key must be a valid key
            in every ScanFrame.
        ds_factor: int
            The factor by which to downsample.  `ds_factor` = 2 will give
            timestreams that are half their original length.
        decimate_keys: list, None
            A list of keys to decimate.  Decimation is down without an
            anti-aliasing filter.  Usually things like pointing go here.
            If `None`, the list will be populated automatically from all the 
            keys that are G3Timestreams in the first Scan frame.
            If you don't want to decimate any keys, set `decimate_keys` = [].
        compress: bool, True
            If true, then turn on flac compression for `keys`.  This
            only works on timestreams in units of counts.
        key_suffix: str, '_downsampled'
            This string is added to the end of each key in `keys` and
            `decimate_keys`, and the downsampled/decimated data is stored
            in the new key.
        temp_suffix: str, '_DS_TEMP'
            Internal temporary keys will be stored with this suffix.
            You probably don't need to worry about this.
        debug: bool, False
            Include debugging output.  You probably don't want this.
        '''
        self.ds_factor = int(ds_factor)
        if isinstance(keys, str):
            keys = [keys]
        self.keys = keys
        if isinstance(decimate_keys, str):
            decimate_keys = [decimate_keys]
        self.decimate_keys = decimate_keys
        self.compress = compress
        self.key_suffix = key_suffix
        self.temp_suffix = temp_suffix
        self.debug = debug
        self.i_start = 0
        
    def __call__(self, fr):
        # Do the actual downsampling and decimating
        if fr.type != core.G3FrameType.Scan:
            return
        if self.decimate_keys is None:
            # Find keys to decimate
            self._find_decimate_keys(fr)
        for k in self.keys:
            new_k = k + self.key_suffix
            outmap_last = fr.pop(new_k + self.temp_suffix, None)
            outmap_ds = core.G3TimestreamMap()
            for subkey, value in outmap_last.iteritems():
                outmap_ds[subkey] = value[self.i_start::self.ds_factor]
                outmap_ds[subkey].SetFLACCompression(self.compress)
            fr[new_k] = outmap_ds
            # fr[new_k + self.temp_suffix] = outmap
            if self.debug:
                fr[new_k + '_debug'] = outmap_last
        for k in self.decimate_keys:
            if k not in fr:
                continue
            # Step 5
            new_k = k + self.key_suffix
            if isinstance(fr[k], core.G3Timestream):
                fr[new_k] = core.G3Timestream(fr[k][self.i_start::self.ds_factor])
                fr[new_k].start = fr[k].start
                fr[new_k].stop = fr[k].stop
            else:
                if k != 'DetectorSampleTimes':
                    core.log_warn(('Decimating frame[{}], which is not a timestream.  ' + \
                                   'This is not normal').format(k))
                fr[new_k] = type(fr[k])(np.array(fr[k])[self.i_start::self.ds_factor])
        self.i_start = (self.ds_factor - (len(value) - self.i_start) % self.ds_factor) % self.ds_factor
        return fr

    def _find_decimate_keys(self, frame):
        """
        Find and cache all auxilliary timestreams that require
        simple decimation to match the bolometer timestream sample rate.
        In addition, add the special G3VectorTime "DetectorSampleTimes",
        if it exists.
        """
        self.decimate_keys = []
        for k in frame.keys():
            if k in self.decimate_keys:
                continue

            # avoid errors for unrecognized object types
            try:
                v = frame[k]
            except:
                pass
            else:
                # field must be a timestream instance
                if isinstance(v, core.G3Timestream):
                    self.decimate_keys += [k]
        if ('DetectorSampleTimes' in frame and 
            isinstance(frame['DetectorSampleTimes'], core.G3VectorTime)):
            # This is a special instance of a G3VectorTime that contains
            # the actual timestamps of samples, as reported by the iceboards.
            # Hence, we want to decimate it with the appropriate phase.
            self.decimate_keys += ['DetectorSampleTimes']

    def cleanup(self, frame):
        """
        Move derived downsampled fields to their fullrate names,
        so that the downsampled frame looks otherwise identical to
        the original fullrate one.

        Also remove Q data for any downsampled I fields, if the Q data
        are not also downsampled.
        """
        if frame.type != core.G3FrameType.Scan:
            return

        for key in self.keys:
            dkey = key + self.key_suffix
            data = frame.pop(dkey, None)
            if data is None:
                raise ValueError('Missing downsampled data!')
            del frame[key]
            if key.endswith('_I'):
                qkey = '{}_Q'.format(key[:-2])
                if qkey not in self.keys and qkey in frame:
                    del frame[qkey]
            frame[key] = data

        for key in self.decimate_keys:
            if key in frame:
                dkey = key + self.key_suffix
                val = frame.pop(dkey, None)
                if val is None:
                    raise ValueError(
                        'Missing downsampled {} data'.format(dkey))
                del frame[key]
                frame[key] = val

@core.pipesegment
def Downsampler(pipe, ds_factor = 2, keys = ['RawTimestreams_I'], 
                decimate_keys = None, compress = True, keep_input_ts = False,
                key_suffix = '_downsampled', temp_suffix = '_DS_TEMP', 
                debug = False):
    '''
    This pipesegment applies an anti-aliasing filter and downsamples data.
    The only work the pipesegment does is to build the kernel.
    It takes one argument that is not passed on to the sub-modules:

    keep_input_ts: bool, False
        If `True`, the original timestreams will be kept, and the 
        downsampled outputs will be stored in <key>_`key_suffix`.
        Otherwise, the downsampled output will overwrite the original
        timestreams.
    '''
    kernel_len = 16 * 2**int(np.ceil(np.log(ds_factor)/np.log(2))) + 1
    kernel = timedomain_sinc(kernel_len, 1. / ds_factor)
    filtered_fullrate_suffix = key_suffix + temp_suffix + 'internal'
    cf = ConvolveFilter(filter = kernel, keys = keys, 
                        key_suffix = key_suffix + temp_suffix, 
                        temp_suffix = filtered_fullrate_suffix,
                        debug = debug)
    pipe.Add(cf.drop_short_scans)
    pipe.Add(cf)
    lad = LPFAndDecimate(keys = keys, ds_factor = ds_factor,
                                decimate_keys = decimate_keys, 
                                compress = compress, key_suffix = key_suffix, 
                                temp_suffix = temp_suffix, debug = debug)
    pipe.Add(lad)
    if not (debug or keep_input_ts):
        pipe.Add(lad.cleanup)

