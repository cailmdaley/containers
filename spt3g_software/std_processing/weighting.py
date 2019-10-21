from spt3g import core, dfmux, todfilter
import numpy as np 
import scipy.stats

@core.indexmod
def AddPSDWeights(frame, input='PolyFilteredTimestreams', output='TodWeights',
                  low_f = 3.*core.G3Units.Hz, high_f = 5.*core.G3Units.Hz,
                  store_psd = False, use_stored_psd = False,
                  store_freq_key = 'PsdFreqs', store_psd_key = 'Psd',
                  store_power = False, store_power_key = 'PsdPower'):
    '''
    Adds PSD based weights for timestreams. Integrates timestream power between
    3-5 Hz (default, can set this range by flagging low_f and high_f), and 
    inversely weights based on that.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    if use_stored_psd and store_psd_key in frame and store_freq_key in frame:
        psds, freqs = frame[store_psd_key], frame[store_freq_key]
    else:
        psds, freqs = todfilter.dftutils.get_psd_of_ts_map(frame[input])
    w = core.G3MapDouble()
    power = core.G3MapDouble()
    low_idx = np.where(np.asarray(freqs) > low_f)[0][0]
    hi_idx = np.where(np.asarray(freqs) < high_f)[0][-1]
    for k,ts in psds.items():
         #integrates timestream power between low_f and high_f
        noise_power = np.trapz(np.asarray(ts)[low_idx:hi_idx],
                               x = np.asarray(freqs)[low_idx:hi_idx])
        power[k] = noise_power
        w[k] = 1./noise_power
        if not np.isfinite(w[k]):
            w[k] = 0
    frame[output] = w
    if store_psd:
        frame[store_freq_key] = core.G3VectorDouble(freqs)
        frame[store_psd_key] = core.G3MapVectorDouble(psds)
    if store_power:
        frame[store_power_key] = power

@core.indexmod
def AddVarWeight(frame,input='PolyFilteredTimestreams', output='TodWeights',
                 store_var = False, store_var_key = 'TodVariances'):
    '''
    Add weights based on timestream variance. 
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    tsm = frame[input]
    weights = core.G3MapDouble()
    variances = core.G3MapDouble()
    for k in tsm.keys():
        v = np.var(tsm[k])
        variances[k] = v
        if np.isfinite(v) and v!=0:
            weights[k] = 1.0/v
        else:
            weights[k] = 0
    frame[output] = weights
    if store_var:
        frame[store_var_key] = variances

@core.indexmod
def AddMadVarWeight(frame,input='PolyFilteredTimestreams', output='TodWeights',
                    store_var = False, store_var_key = 'TodVariances'):
    '''
    Add weights based on the MAD (median absolute deviation) estimator. This is
    more robust to outliers (bright point sources) than PSD-based estimation.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    ts = frame[input]
    w = core.G3MapDouble()
    variances = core.G3MapDouble()
    for k in ts.keys():
        v = todfilter.get_mad_std(ts[k])**2.0
        variances[k] = v
        if v != 0:
            w[k] = 1.0/v
        else:
            w[k] = 0
        if not np.isfinite(w[k]):
            w[k] = 0
    frame[output] = w
    if store_var:
        frame[store_var_key] = variances

@core.indexmod
def AddMaskedVarWeight(frame, input='PolyFilteredTimestreams',
                       output='TodWeights', mask_key = 'FilterMask',
                       store_var = False, store_var_key = 'TodVariances'):
    '''
    Add weights based on masked variance. Pixels that see the point source in a
    given scan are masked (excluded) from the variance estimator. Requires the
    FilterMask key added by mapmaker.todfilter.
    '''
    if input not in frame or mask_key not in frame:
        return
    tsm = frame[input]
    has_mask = frame[mask_key].has_masked_pixels
    mask = frame[mask_key].pixel_mask
    weights = core.G3MapDouble()
    variances = core.G3MapDouble()
    for k in tsm.keys():
        if has_mask[k]:
            v = np.var(tsm[k][np.logical_not(np.array(mask[k]).astype(bool))])
        else:
            v = np.var(tsm[k])
        variances[k] = v
        if np.isfinite(v) and v != 0:
            weights[k] = 1.0/v
        else:
            weights[k] = 0.
    frame[output] = weights
    if store_var:
        frame[store_var_key] = variances

@core.indexmod
def AddMaskedMadVarWeight(frame, input='PolyFilteredTimestreams',
                          output='TodWeights', mask_key = 'FilterMask',
                          store_var = False, store_var_key = 'TodVariances'):
    '''
    Add weights based on masked MAD variance. Pixels that see the point source
    in a given scan are masked (excluded) from the MAD variance estimator.
    Requires the FilterMask key added by mapmaker.todfilter.
    '''
    if input not in frame or mask_key not in frame:
        return
    tsm = frame[input]
    has_mask = frame[mask_key].has_masked_pixels
    mask = frame[mask_key].pixel_mask
    weights = core.G3MapDouble()
    variances = core.G3MapDouble()
    for k in tsm.keys():
        if has_mask[k]:
            v = todfilter.get_mad_std(core.G3Timestream(
                tsm[k][np.logical_not(np.array(mask[k]).astype(bool))]) )
        else:
            v = todfilter.get_mad_std(tsm[k])
        variances[k] = v**2.0
        if np.isfinite(v) and v != 0:
            weights[k] = 1.0/v**2.0
        else:
            weights[k] = 0.
    frame[output] = weights
    if store_var:
        frame[store_var_key] = variances

@core.indexmod
def AddSigmaClippedWeight(frame, input='PolyFilteredTimestreams',
  output='TodWeights', threshold=2.5):
    '''
    Add weights based on a sigma-clipped RMS. This is more robust to point
    sources. Weights on computed based on total variance after removing
    points more than <threshold> sigma from the mean.
    ''' 

    if input not in frame:
        return
    w = core.G3MapDouble()
    for k,ts in frame[input].iteritems():
        if not np.isfinite(ts).all():
            w[k] = 0
        elif (ts == 0).all():
            w[k] = 0
        else:
            w[k] = 1./np.var(scipy.stats.sigmaclip(ts, threshold, threshold).clipped)
            if not np.isfinite(w[k]):
                w[k] = 0
    frame[output] = w

