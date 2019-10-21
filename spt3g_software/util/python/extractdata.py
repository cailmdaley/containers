"""
This file contains classes that make interacting with data in G3 files 
easier.
"""

from collections import deque
import numpy as np
from spt3g import core, gcp

@core.indexmod
class ExtractDataFromPipeline(object):
    def __init__(self, keys):
        '''
        This module just stores the last values of the frame data mapped to 
        by keys in the pipeline.  The keys are set to members of this object.
        This is meant to get data out of the pipeline for interactive use and 
        not for production code.

        Example::
            data = ExtractDataFromPipeline(keys=["BolometerProperties","WiringMap"])
            pipe.Add(data)
            pipe.Run()

            data.BolometerProperties

        '''
        assert('keys' not in keys) 
        self.keys = keys
    def __call__(self, frame):
        for k in self.keys:
            if k in frame:
                self.__setattr__(k,frame[k])

@core.indexmod
class ExtractCombinedTimestreams(object):
    '''
    This module concatenates G3Timestreams or G3TimestreamMaps in the Scan 
    frames of the pipe.  keys is a list of the keys you want.  These values 
    are stored as the members of the object.

    This is meant to get data out of the pipeline for interactive use and not 
    for production code.

    Example::
        data = ExtractCombinedTimestreams(keys = ['CalTimestreams'])
        pipe.Add(data)
        pipe.Run()

        data.CalTimestreams

    '''

    def __init__(self, keys):
        assert('keys' not in keys) 
        self.keys = keys
        
    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return
        for k in self.keys:
            if not k in self.__dict__:
                self.__setattr__(k, frame[k])
            else:
                self.__setattr__(k, core.concatenate_timestreams([self.__getattribute__(k), 
                                                                  frame[k]]))

@core.indexmod
class Accumulator(object):
    '''
    A handy module for grabbing a key (and possibly subkeys) from frames.
    Extracted data is available in Accumulator.values
    '''
    def __init__(self, keys, callback = None, frametype = None):
        '''
        keys: str or list of str
            The keys (and subkeys) you want to extract.  For example,
            keys = ['key', 'subkey1', 'subkey2']
            will result in Accumulator.values containing the value of
            fr['key']['subkey1']['subkey2'] for each frame.
        callback: callable
            `callback` is run on the object extracted from each frame.
            For example, callback = lambda az: az / core.G3Units.deg
            would convert azimuth to degrees before storing it in
            Accumulator.values.
        frametype: core.G3FrameType
            Only extrac objecst from frames of this type.        
        '''
        # For REASONS, the spt3g arc file reader doesn't know how to parse
        # registers that contain strings.  Since those registers aren't hard
        # coded in ARCFileReader, we have to do it here.  Ditto for bools
        # Obviously, this list is incomplete.
        self.string_keys = ['scan_name']
        self.bool_keys = ['off_source']
        self.keys = keys
        if isinstance(keys, str):
            self.keys = [keys]
        else:
            self.keys = keys
        self.values = deque()
        if callback is None:
            self.callback = lambda v: v
        else:
            self.callback = callback
        self.frametype = frametype

    def __call__(self, fr):
        def fix_arc_string(val):
            if isinstance(val, core.G3String):
                return val
            return core.G3String(''.join([chr(v) for v in val]).split(chr(0))[0])
        def fix_arc_bool(val):
            if isinstance(val, core.G3String):
                # FYI, I'm not 100% sure about this condition
                return len(val) > 0
            if isinstance(val, core.G3VectorInt):
                return val[0] > 0

        if fr.type == core.G3FrameType.EndProcessing:
            return
        if self.frametype is not None and fr.type != self.frametype:
            return
        try:
            thingy = fr[self.keys[0]]
        except KeyError:
            return fr
        for k in self.keys[1:]:
            try:
                thingy = thingy[k]
            except KeyError:
                return fr
        if self.keys[-1] in self.string_keys:
            thingy = fix_arc_string(thingy)
        self.values.append(self.callback(thingy))
        return fr

    def extract_values(self):
        if len(self.values) == 0:
            core.log_warn('{} is empty, and probably missing'.format(self.keys))
            return np.array([])
        if (isinstance(self.values[0], core.G3Timestream) or 
            isinstance(self.values[0], core.G3TimestreamMap)):
            self.values = list(self.values)
            arr = self.values[0].concatenate(self.values[1:])
        else:
            arr = np.array(self.values).flatten()
            if arr.dtype == np.dtype('O'):
                for i, v in enumerate(arr):
                    arr[i] = convert_G3_types(v)
            shape = np.shape(self.values[0])
            if len(shape) > 0 and shape[0] != 1:
                if len(shape) == 1:
                    if shape[0] == 100:
                        return arr
                    else:
                        arr = np.concatenate(arr)
                else:
                    arr = arr.reshape((len(arr) // (shape[0] * shape[1]),
                                       shape[1] * shape[0])).T
        return arr

@core.indexmod
class MultiAccumulator(object):
    '''
    Like Accumulator, but extracts many keys.
    '''
    def __init__(self, keys, frametype = None, **kwargs):
        '''
        Arguments to this init are all dictionaries (or None).  Each key will 
        generate one Accumulator, and all the optional arguments to Accumulator
        will be filled in from optional arguments to this init.  The individual
        Accumulators are available in MultiAccumulator.accums

        keys: dict
            Each key should map to a frame key or list of keys that will be 
            passed to Accumulator instances.
        kwargs:
            All further arguments are passed to Accumulator.  If a dict, must
            have the same keys as `keys`.  If anything else, the argument will 
            be passed as-is to every Accumulator instance.
        '''
        self.accums = {}
        for k in keys.keys():
            this_kwargs = {}
            for kw in kwargs.keys():
                try:
                    this_kwargs[kw] = kwargs[kw][k]
                except TypeError:
                    this_kwargs[kw] = kwargs[kw]
            self.accums[k] = Accumulator(keys[k], **this_kwargs)
        
    def __call__(self, frame):
        for acc in self.accums.values():
            acc(frame)

    def extract_values(self):
        '''
        For each Accumulator, process the collected data into a single 
        array-like object.
        '''
        out = {}
        for k, acc in self.accums.items():
            out[k] = acc.extract_values()
        return out
        

@core.usefulfunc
def extract_keys(file, keys, frametype = None, verbose = False, 
                 arc_extract = True):
    '''
    This is a wrapper around Accumulator to extract one or more keys
    from a list of files.

    INPUTS
    ------
    file: str or list
        The filename or list of filenames from which to extract information
    keys: str, list or dict
        If `keys` is a string or list, only that key (or the specified subkey)
        are extracted from the frames.  The output will be the concatenated
        data from the specified key.
        If `keys` is a dict, each key in `keys` should map to a string
        or list of strings that specify the key or subkey(s) to be extracted.
        The return from `exract_keys` will be a dictionary, with the 
        concatenated data from each key.
    frametype: core.G3FrameType
        Only extrac objecst from frames of this type.
    verbose: bool
        If True, print every frame
    arc_extract: bool
        If the input files are archive files, and `arc_extract` is True,
        run extra processing (gcp.ARCExtractor) to make the output slightly 
        nicer.

    Example::
        >>> extract_keys(filename, keys={'az': 'OnlineBoresightAz', 'bolo': ['RawTimestreams_I', <boloname>]}
        {'az': <the azimuth for all frames>, 'bolo': <the timestream for <boloname> >}

    TODO: respect G3 types instead of dumping into np.array
    '''
    if isinstance(file, str):
        if file.endswith('.g3'):
            arcfile = False
            reader = core.G3Reader(file)
        elif file.endswith('.dat') or file.endswith('.dat.gz'):
            arcfile = True
            reader = gcp.ARCFileReader(file)
        else:
            raise ValueError('Extension not recognized on {}'.format(file))
    elif np.iterable(file):
        if file[0].endswith('.g3'):
            arcfile = False
            reader = core.G3Reader(file)
        elif file[0].endswith('.dat') or file[0].endswith('.dat.gz'):
            arcfile = True
            reader = gcp.ARCFileReader(file)
        else:
            raise ValueError('Extension not recognized on {}'.format(file[0]))
    else:
        arcfile = False
        reader = file
    if isinstance(keys, dict):
        accums = {}
        for k in keys:
            accums[k] = Accumulator(keys[k], frametype = frametype)
    else:
        accums = {'value': Accumulator(keys, frametype = frametype)}
    pipe = core.G3Pipeline()
    pipe.Add(reader)
    if arcfile and arc_extract:
        pipe.Add(gcp.ARCExtract)
    if verbose:
        pipe.Add(core.Dump, added_message = 'Before accumulate')
    for a in accums.values():
        pipe.Add(a)
    if verbose:
        pipe.Add(core.Dump, added_message = 'After accumulate')
    pipe.Run()
    out = {}
    for k, acc in accums.items():
        out[k] = acc.extract_values()
    if 'value' in out.keys() and len(out.keys()) == 1:
        out = out['value']
    return out

def convert_G3_types(a):
    '''
    Helper function that attempts to convert G3 data types
    '''
    try:
        if np.iterable(a):
            a = np.array([v.value for v in a])
        else:
            a = a.value
    except AttributeError:
        pass
    return a        
