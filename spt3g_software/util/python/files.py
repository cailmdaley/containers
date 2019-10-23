"""
This contains functions for reading/writing specific kinds of files
"""
import warnings
import numpy as np
import os
from glob import glob
try:
    import cPickle as pickle
except ImportError:
    import pickle
from spt3g import core, gcp
from .extractdata import extract_keys, MultiAccumulator

@core.usefulfunc
def save_pickle(data, filename):
    """
    Save a backward-compatible pickle file with the given filename.
    Ensures that the file is open in mode 'wb' (required for python3), and
    the most backward-compatible protocol is used (so that new files can still
    be opened in python2).
    """
    if hasattr(filename, 'write'):
        if filename.mode == 'wb':
            pickle.dump(data, filename, protocol=0)
            return
        warnings.warn("Reopening file {} in mode 'wb' for pickling".format(
            filename.name))
        filename.close()
        filename = filename.name
    with open(filename, 'wb') as f:
        pickle.dump(data, f, protocol=0)

@core.usefulfunc
def load_pickle(filename):
    """
    Load a pickle file from the given filename.
    Ensure that the file is open in mode 'rb' (required for python3), and
    that the encoding is set to 'bytes' in python3.
    """
    if hasattr(filename, 'read'):
        if filename.mode == 'rb':
            try:
                return pickle.load(f, encoding='latin1')
            except TypeError:
                return pickle.load(f)
        warnings.warn("Reopening file {} in mode 'rb' for unpickling".format(
            filename.name))
        filename.close()
        filename = filename.name
    with open(filename, 'rb') as f:
        try:
            return pickle.load(f, encoding='latin1')
        except TypeError:
            return pickle.load(f)

        
@core.usefulfunc
def get_fridge_log(start_time = None, stop_time = None, obsid = None,
                   source = '*', basedir = '/spt/data/arc'):
    """
    Given a start and stop time (as core.G3Time), or an observation id,
    will find the relevant fridge log and return it as a dictionary.

    Warning: this can take anywhere from 1 to 10 minutes.
    """
    assert any((obsid != None , (start_time != None and
                                 stop_time != None)))
    if obsid is not None:
        # Grab start and end times of obs
        if not isinstance(obsid, str):
            obsid = str(obsid)
        data_dir = '/spt/data/bolodata/downsampled/'
        obs_datafiles = sorted(glob(os.path.join(data_dir, source,
                                                 obsid, '0*.g3')))
        first = core.G3File(obs_datafiles[0])
        for fr in first:
            if 'ObservationStart' in fr:
                start_time = fr['ObservationStart']
                stop_time = fr['ObservationStop']
                break

    # Find the relevant arc files
    from spt3g.std_processing.FindARCFiles import ARCForTimerange
    arcfiles = ARCForTimerange(start_time, stop_time, basedir)

    # Set up the He10 fridge cryoboard register indices
    registers = {0:'UCHEAD', 1:'ICHEAD', 2:'HE4HEAD', 3:'HE4FB',
                4:'HE4PUMP', 5:'ICPUMP', 6:'UCPUMP', 7:'HE4SW',
                8:'ICSW', 9:'UCSW', 10:'UC_STAGE', 11:'LC_TOWER',
                12:'IC_STAGE', 13:'4K_HEAD', 14:'4K_SQUID_STRAP',
                15:'50K_HEAD'}
    keys = {}
    for reg, therm in registers.items():
        keys[therm] = ['array','cryo', 'temperature', 0, reg]
        keys['Time'] = ['array','cryo','utc']

    # Initialize the MultiAccumulator object
    fridge_temps = MultiAccumulator(keys = keys)

    # Start the pipeline
    pipe = core.G3Pipeline()
    pipe.Add(gcp.ARCFileReader, filename = arcfiles)
    pipe.Add(fridge_temps)
    pipe.Run()

    # Return dictionary
    return fridge_temps.extract_values()
