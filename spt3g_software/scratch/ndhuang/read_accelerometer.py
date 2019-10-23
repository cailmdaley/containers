'''
Read data from accelerometers mounted on the receiver in early 2018.
Absolute time is not guaranteed to better than ~.1s, however the sampling
rate should be quite even.
'''

import os
from datetime import datetime
from pytz import UTC
import struct
import numpy as np

headerfmt = 'sII'
headersize = struct.calcsize(headerfmt) + struct.calcsize('d')

def _read_header(file):
    headerstr = file.read(headersize)
    encoding, nchan, srate = struct.unpack_from(headerfmt, headerstr)
    start = struct.unpack_from('d', headerstr, headersize - 8)[0]
    print(encoding, nchan, srate, start)
    return encoding, nchan, srate, start

def _read_data(file, encoding, nchan):
    datstr = file.read()
    datstr = datstr[headersize:]
    sampsize = struct.calcsize(encoding)
    nsamps = (len(datstr) - headersize) // sampsize
    if isinstance(encoding, str):
        # py2 treats bytes and strings as the same data type
        fmt = str(nsamps) + encoding
    else:
        # py3 doesn't
        fmt = str(nsamps) + str(encoding)[2]
    data = struct.unpack_from(fmt, datstr, headersize)
    # drop partial samples
    if len(data) % nchan != 0:
        data = data[:-(len(data) % nchan)]
    return np.array(data).reshape(nsamps // nchan, nchan)
    

def read_file(fname):
    '''
    Read one file's worth of accelerometer data.
    Returns data (nsamples x nchannels), sample rate and start time (ctime).
    '''
    with open(fname, 'rb') as f:
        encoding, nchan, srate, start = _read_header(f)
        data = _read_data(f, encoding, nchan)
    return data, srate, start

def get_times(srate, data, start, return_datetime = True):
    '''
    Generate sample times from accelerometer data.
    
    INPUTS
    ------
    srate: float
        The sample rate, in Hz
    data: int, array-like
        Either the data array, or the number of samples in the data array.
    start: float,
        Start time (ctime)
    return_datetime: bool
        If True, return times as datetimes, otherwise as ctimes
    '''
    try:
        nsamp = len(data)
    except TypeError:
        nsamp = data
    dt = 1. / srate
    times = start + np.arange(nsamp) * dt
    if return_datetime:
        times = np.array([datetime.utcfromtimestamp(t) for t in times])
    return times

def read_for_time(start, stop, 
                  datadir = '/spt/data/rsync/accelerometer_data/'):
    '''
    Read acclerometer data from `start` to `stop`, where `start` and
    `stop` are datetimes.

    Returns a tuple of `(acclerometer_data, time, samplerate)`, where
    `accelerometer_data` is an array of 3 timestreams in volts.
    `time` is a list of datetime objects.
    `samplerate` is in Hz.

    NB: Relative timing (i.e. sampling rate) should be very accurate and consistent.  
    However, the absolute time may be off by ~.2s globally.
    '''
    if isinstance(start, datetime):
        epoch = datetime(1970, 1, 1, 0, tzinfo = UTC)
        start = (start.replace(tzinfo = UTC) - epoch).total_seconds()
        stop = (stop.replace(tzinfo = UTC) - epoch).total_seconds()
    files = sorted(os.listdir(datadir))
    # find the files we need
    for i, fname in enumerate(files):
        with open(os.path.join(datadir, fname), 'rb') as f:
            enc, n, srate, f_start = _read_header(f)
        if f_start > stop and i == 0:
            print(datetime.fromtimestamp(f_start))
            print(datetime.fromtimestamp(stop))
            raise RuntimeError('No data in your requested interval')
        if f_start < start:
            i_start = i
            data_start = f_start
        if f_start < stop:
            i_stop = i + 1
    files = files[i_start:i_stop]
    data = []
    srate_last = None
    for fname in files:
        _data, srate, junk = read_file(os.path.join(datadir, fname))
        if srate_last is not None:
            assert(srate == srate_last)
        srate_last = srate
        data.append(_data)
    data = np.concatenate(data)
    times = get_times(srate, len(data), data_start, False)
    inds = np.where((times <= stop) & (times >= start))[0]
    data = data[inds]
    times = get_times(srate, len(data), times[inds][0])
    return data, times, srate

def lpf(data, srate, cutoff):
    '''
    Low pass filter a data stream.  Uses an unwindowed sinc filter in the 
    Fourier domain.
    
    INPUTS
    ------
    data: array-like
        data stream to be filtered
    srate: float
        sample rate
    cutoff: float
        cutoff frequency.  All frequencies above this are zeroed.
    '''
    n = len(data)
    # window = .54 - .46 * np.cos(2 * np.pi * np.arange(n) / (n - 1)) # Hamming
    window = np.ones_like(data)
    ft = np.fft.rfft(window * data)
    freq = np.fft.rfftfreq(n, 1. / srate)
    ft[freq > cutoff] = 0
    return np.fft.irfft(ft)
