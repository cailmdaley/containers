import numpy as np
from spt3g import core, todfilter
from spt3g.timestreamflagging import add_flag
from copy import copy

U = core.G3Units

def get_unphysical_variance(tod_units, sample_rate):
    '''
    Returns an unphysically low variance value for use with flagging
    '''
    if tod_units == core.G3TimestreamUnits.Current:
        return ( 5e-12 * U.amp )**2.0 / U.hz * sample_rate
    elif tod_units == core.G3TimestreamUnits.Tcmb:
        return ( 200 * U.uK )**2.0 / U.hz * sample_rate #i don't even know any more
    elif tod_units == core.G3TimestreamUnits.Power:
        return ( 30 * U.attowatt )**2.0 / U.hz * sample_rate
    else:
        raise RuntimeError("The timestream units " + str(tod_units) +
                           " have yet to be included in unphysical variance list")

@core.indexmod
def FlagUnphysicallyLowVariance(frame, ts_key, prefilter_poly_order = -1,
                                flag_key = 'Flags', plot_vals = False):
    '''
    Flags detectors with unphysically low variances as determined by the function
    get_unphysical_variance.  If  prefilter_poly_order > 0, applies a poly filter
    to the timestream data before estimating the variance.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    ts_map = frame[ts_key]
    if prefilter_poly_order >= 0:
        ts_map = todfilter.polyutils.poly_filter_g3_timestream_map(
            ts_map, prefilter_poly_order)
    bad_bolos = []

    for k, v in ts_map.iteritems():
        var = np.var(v)
        if var < get_unphysical_variance(v.units, v.sample_rate):
            bad_bolos.append(k)
    add_flag(frame, flag_key, 'UnphysicalLowVariance', bad_bolos)

@core.indexmod
def FlagUnphysicallyLowCustomVariance(frame, ts_key, var_key, flag_key = 'Flags'):
    '''
    Flags detectors with unphysically low variance as determined by the function
    get_unphysical_variance. The variance has to be precomputed by another module
    and stored in var_key.
    '''
    if var_key not in frame:
        return
    var_map = frame[var_key]
    ts_map = frame[ts_key]
    bad_bolos = []
    for k, var in var_map.iteritems():
        if var < get_unphysical_variance(ts_map[k].units, ts_map[k].sample_rate):
            bad_bolos.append(k)
    add_flag(frame, flag_key, 'UnphysicalLowVariance', bad_bolos)

@core.indexmod
def FlagOscillatingChannels(frame, ts_key, threshold = 5, flag_key = 'Flags'):
    '''
    Flags detectors with excess power in the 6-10 Hz or 8-12 Hz bins
    compared to the 2-6 Hz bin, symptomatic of oscillating bolometers.
    '''
    if ts_key not in frame:
        return
    psd_map, freqs = todfilter.dftutils.get_psd_of_ts_map(frame[ts_key])
    low_bin = np.where((freqs/core.G3Units.Hz > 2) &
                       (freqs/core.G3Units.Hz < 6))[0]
    high_bin1 = np.where((freqs/core.G3Units.Hz > 6) &
                         (freqs/core.G3Units.Hz < 10))[0]
    high_bin2 = np.where((freqs/core.G3Units.Hz > 8) &
                         (freqs/core.G3Units.Hz < 12))[0]
    bad_bolos = []
    for bolo in psd_map.keys():
        psd = np.array(psd_map[bolo])
        if (np.sum(psd[high_bin1])/np.sum(psd[low_bin]) > threshold):
            bad_bolos.append(bolo)
        elif (np.sum(psd[high_bin2])/np.sum(psd[low_bin]) > threshold):
            bad_bolos.append(bolo)
    add_flag(frame, flag_key, 'Oscillating',bad_bolos)