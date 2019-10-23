import numpy as np
from spt3g import core
from spt3g.todfilter import rolling_mean_filter_ts
from spt3g.todfilter.util import mask_timestream
from spt3g.timestreamflagging import add_flag

from copy import copy

@core.usefulfunc
def get_num_glitches(ts_unmasked, thresholds, ts_masked = None,
                     kernel_size = 11, also_get_hist=False, hist_bins=None):
    '''
    Flag detectors which have glitches exceeding the threshold.  For a map should be 
    called on the individual scans, so that a single glitched scan doesn't result
    in a flag for the entire bolo.  Also, glitches in turnarounds and el steps won't be
    considered.
    
    INPUTS
    ------
    kernel_size : int
        The window width for the tophat filter.
    thresholds : array-like
        An array of sigma values, data points above these sigma values
        are glitches

    Returns
    -------
    number of glitches above a threshold.  If thresholds is a list/array returns 
    this in a list
      
    History:
        created S. Patil
        24 May 2016: changed output lines: CR
        Modified by nlharr to fix format
        July 2019: Modified by kferguson to allow for point source masking.
    '''
    if np.shape(thresholds) == ():
        thresholds = [thresholds]

    if ts_masked is not None:
        ts = copy(ts_masked)
    else:
        ts = copy(ts_unmasked)

    if kernel_size > 0:
        outliers = rolling_mean_filter_ts(ts, kernel_size)
    else:
        outliers = copy(ts)

    rms = np.nanstd(outliers)

    if np.any(np.isnan(ts_unmasked)) or rms == 0:
        out_nglitches = []
        for thresh in thresholds:
            out_nglitches.append(len(outliers))
    else:
        outliers = np.abs(outliers/rms)
        out_nglitches = []
        for thresh in thresholds:
            out_nglitches.append(int(np.sum(outliers > thresh)))
    if also_get_hist:
        if hist_bins is None:
            hist_bins = range(0,30)
        h = np.histogram(outliers[~np.isnan(outliers)],bins=hist_bins)
        above = np.size(np.where(outliers > max(hist_bins)))
        return out_nglitches, h, above
    else:
        return out_nglitches


@core.indexmod
def AddNumGlitches(frame, thresholds, 
                   input_ts_key, 
                   output_glitch_num_key = 'GlitchesNumberOf',
                   output_thresholds_key = 'GlitchesThresholds',
                   kernel_size = 11, 
                   plot_hist_normed_ts=False, save_dir=None,
                   plot_scan_tag_key = 'ScanId',
                   filter_mask_key = None):
    '''
    Adds the number of "glitches" for the timestreams over the scan.

    These found glitches are stored in a vector corresponding to the thresholds.
    So if you provide 3 thresholds [5, 10, 20], the output_glitch_num_key will store a vector
    with 3 elements, [the number of glitches above 5 sigma, 10 sigma, 20 sigma]

    Glitches are found by heavily highpass filtering the data and then searching
    for deviations N sigma away where sigma is estimated from the standard deviation
    of the source assuming a Normal distribution.


    thresholds: An iterable of the "N sigma" thresholds we are using as cutoffs
    input_ts_key [->G3TimestreamMap] The timestreams we care about
    output_glitch_num_key [->G3MapVectorInt]: where we store the glitches the index maps to 
                 the thresholds list
    output_thresholds_key [->G3VectorDouble]: A copy of the thresholds are stored in the frame
    
    kernel_size: Right now the highpass filter is a rolling mean filter with width specified
                 by the kernel size
    '''
    if (plot_hist_normed_ts):
        assert(not save_dir is None)

    if frame.type != core.G3FrameType.Scan:
        return 

    ts_map = frame[input_ts_key]
    out_map = core.G3MapVectorInt()
    thresholds = core.G3VectorDouble(thresholds)

    if plot_hist_normed_ts:
        h_bins = None
        h_sum = 0
        n_above_sum = 0

    if filter_mask_key is not None:
        ps_map = frame[filter_mask_key].pixel_mask
        #is_masked = frame[filter_mask_key].has_masked_pixels
    else:
        ps_map = None

    for k in ts_map.keys():

        if ps_map is not None:
            ts_masked = mask_timestream(ts_map[k], ps_map[k], mask_with = np.nan)
        else:
            ts_masked = None

        if plot_hist_normed_ts:
            nglitches,hist,n_above = get_num_glitches(ts_map[k], thresholds,
                                                      ts_masked = ts_masked,
                                                      kernel_size = kernel_size,
                                                      also_get_hist=True)
            n_above_sum += n_above
            h_bins = hist[1]
            h_sum = h_sum+hist[0]
        else:
            nglitches = get_num_glitches(ts_map[k], thresholds,
                                         ts_masked = ts_masked,
                                         kernel_size = kernel_size)
        out_map[k] = nglitches
    frame[output_glitch_num_key] = out_map
    frame[output_thresholds_key] = thresholds

    if plot_hist_normed_ts:
        import matplotlib.pyplot as plt
        width = 0.7 * (h_bins[1] - h_bins[0])
        center = (h_bins[:-1] + h_bins[1:]) / 2
        plt.clf()
        plt.yscale('log')
        plt.bar(center, h_sum, align='center', width=width,bottom=0.1)
        plt.title("Histogram of abs(TOD/std(TOD)) for all detectors in a scan")
        plt.xlabel("abs(TOD/std(TOD)) where TOD is high passed")
        plt.savefig(save_dir+'/GlitchNormedTod_'+frame[plot_scan_tag_key]+'.png')

@core.indexmod
def FlagGlitches(frame, min_num_glitch_map,
                 glitch_num_key = 'GlitchesNumberOf',
                 glitch_thresholds_key = 'GlitchesThresholds',
                 flag_key = 'Flags', flag_reason = 'Glitchy'):
    '''
    n_glitch_map has the form { float(threshold) : min number of glitches to flag}

    so if n_glitch_map = {5.0: 3} it will flag any detector with 3 glitches above 5.0 sigma
    per scan.  The glitch thresholds need to have been precomputed with AddNumGlitches.
    '''


    if frame.type != core.G3FrameType.Scan:
        return
    n_glitch_map = copy(min_num_glitch_map)

    num_glitches = frame[glitch_num_key]
    glitch_thresholds = frame[glitch_thresholds_key]
    
    n_glitch_allowed = core.G3VectorInt(np.zeros(len(glitch_thresholds), dtype = 'int'))
    glitch_mask = core.G3VectorInt(np.zeros(len(glitch_thresholds), dtype = 'int'))

    for i, gt in enumerate(glitch_thresholds):
        if gt in n_glitch_map:
            n_glitch_allowed[i] = n_glitch_map[gt]
            glitch_mask[i] = 1
            n_glitch_map.pop(gt)
        else:
            n_glitch_allowed[i] = 0
            glitch_mask[i] = 0
    if (len(n_glitch_map) != 0):
        core.log_fatal("not all glitches you want to flag on were precalculated")
    bad_set = set()
    for k, v in num_glitches.iteritems():
        v = np.asarray(v)
        if np.any(((v - n_glitch_allowed + 1) * glitch_mask) > 0):
            bad_set.add(k)
    add_flag(frame, flag_key, flag_reason, bad_set)



