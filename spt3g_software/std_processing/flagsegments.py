from spt3g import core, timestreamflagging, todfilter
from spt3g.timestreamflagging import FlagBadG3MapValue, RemoveFlaggedTimestreams
from spt3g.timestreamflagging import FlagIncompletePixelPairs, FlagTimestreamsWithoutProperties
from spt3g.timestreamflagging import FlagBadHousekeeping, FlagBadRfrac
from spt3g.timestreamflagging import FlagMissing



'''
An assortment of pipe segments for flagging the data for a variety of reasons.
Some first pass pipe segments are:

FieldFlaggingPreKcmbConversion
FieldFlaggingPostKcmbConversion

The former flags any invalid data or data we don't have calibration info for
The latter function contains more general flags.

Both of these pipesegments only set the flag values.  You will need to add a:
timestreamflagging.RemoveFlaggedTimestreams pipesegment after setting the flags

These are not set in stone and probably need a lot more work improving them.

All of the other pipesegments are trying to flag data that fits some broad description.  The
two listed above are composed of them.

The file spt3g_software/timestreamflagging/python/ contains the function:
get_unphysical_variance  These values are extremely conservative.  You are encouraged to
set those values to be less conservative.

'''

def TagScan(frame, ts_map_key):
    '''
    Adds the key ScanId to the frame with the time string of the start of a scan
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    frame['ScanId'] = frame[ts_map_key].start.isoformat().split('.')[0]

@core.pipesegment
def FlagInvalidData(pipe, flag_key, ts_key, 
                          bolo_props="BolometerProperties",
                          wiring_map = 'WiringMap', hk='DfMuxHousekeeping',
                          flag_nan_offsets = False):
    '''
    Removes any timstreams with invalid data.  The definition of invalid data is 
    more strict with this module.

    This checks for data with missing bolometer properties.
    Also it looks at iceboard housekeeping to see the data may 
    not be an accurate representation of the TES state.
    '''
    pipe.Add(timestreamflagging.FlagNaNs,ts_key=ts_key,flag_key=flag_key)
    pipe.Add(timestreamflagging.FlagNegativeDANChannels,ts_key=ts_key,flag_key=flag_key)
    pipe.Add(FlagTimestreamsWithoutProperties, ts_key=ts_key, bolo_props=bolo_props,
             flag_nan_offsets=flag_nan_offsets)
    pipe.Add(FlagBadHousekeeping, ts_key=ts_key, flag_key=flag_key,
             wiring_map=wiring_map, hk=hk)

FlagInvalidDataStrict = FlagInvalidData #for backwards compatibility

@core.pipesegment
def FlagUncalibratable(pipe, ts_key, flag_key='Flags'):
    '''
    Checks that we have calibration entries for the the detectors in ts_key.

    The variable cal_sources lists any sources that could potentially be used
    to calibrate these data.
    '''
    needed_cal_frame_data = ['CalibratorResponse']
    for needed_data in needed_cal_frame_data:
        pipe.Add(FlagTimestreamsWithoutProperties,
                 ts_key=ts_key, flag_key=flag_key,
                 flag_id='Missing'+needed_data,
                 bolo_props=needed_data)
    pipe.Add(timestreamflagging.FlagMissingFluxCalibration, ts_key = ts_key,
            flag_key = flag_key, flag_id = 'MissingFluxCalibration')


@core.pipesegment
def FlagNonResponsive(pipe, flag_key, min_cal_sn=20, min_elnod_sn = 20):
    '''
    Checks elnod and calibrator SN to remove non-responsive detectors.
    '''
    pipe.Add(FlagBadG3MapValue, m_key='CalibratorResponseSN', min_val=min_cal_sn,
             flag_reason="BadCalSn",flag_key=flag_key)
    pipe.Add(FlagBadG3MapValue, m_key='ElnodSNSlopes', min_val=min_elnod_sn,
             flag_reason="BadElnodSn",flag_key=flag_key)


@core.pipesegment
def FlagNoisyNonGaussianData(
    pipe, flag_key, ts_key,

    variance_prefilter_poly_scale = -1,

    glitch_thresholds = [10, 8],
    glitch_num_above = [1, 5],
    glitch_filter_scale = 11,

    max_derivative_val= None,
    deriv_presmooth_scale = 10, 

    glitch_num_key = 'GlitchesNumberOf', glitch_thresholds_key = 'GlitchesThresholds',
    max_derivative_key='MaxDerivative', filter_mask_key = None,
    ps_mask_id = None, ps_pointing_key = None,

    plot_statistics = False, plot_directory = None,
    plot_by_scan = False, plot_n_bins = 100):

    '''
    Performs 3 checks for invalid data from the properties of the timestream.
    1) Are the variances unphyiscally low
    2) Are there an excess of glitches
    3) If we smooth the data do we find any extremely large first derivatives.
    
    The unphysically low variance values are stored in the timestreamflagging function:
    get_unphysical_variance.
    
    Glitch finding is determined by:
    glitch_thresholds, glitch_num_above, glitch_filter_scale

    The TOD is smoothed by a rolling mean filter of sample width glitch_filter_scale.  It
    is then normalized to have a std of 1.  The absolute value is then taken.  For every
    pair of values in glitch_num_above and glitch_thresholds, it counts the number of samples 
    above the threshold.  If that number is larger than the glitch_num_above entry it flags the 
    data.

    For the max derivative, it smooths the data and then finds the 
    abs(max(first derivative of TOD)).  If there a value above the threshold it flags the data
    
    if plot_statistics is true it generates plots of the distribution of max first derivatives
    and the values used by the glitch thresholds and saves them to plot_directory.
    If plot_by_scan is true, generates the maximum derivative values for each scan.

    '''

    if plot_statistics:
        pipe.Add(TagScan,ts_map_key=ts_key)
    #unphysically low variance
    pipe.Add(timestreamflagging.noiseflagging.FlagUnphysicallyLowVariance,
             ts_key = ts_key, prefilter_poly_order = variance_prefilter_poly_scale)
    #glitchy

    if filter_mask_key is not None and ps_mask_id is not None:
        pipe.Add(todfilter.FilterMaskInjector,
                 point_src_mask_id = ps_mask_id,
                 filter_mask_key = filter_mask_key,
                 pointing_key = ps_pointing_key)

    pipe.Add(timestreamflagging.glitchfinding.AddNumGlitches,
             thresholds = glitch_thresholds,
             input_ts_key = ts_key,
             kernel_size = glitch_filter_scale,
             output_glitch_num_key = glitch_num_key,
             output_thresholds_key = glitch_thresholds_key,
             plot_hist_normed_ts=plot_statistics,
             save_dir=plot_directory,
             plot_scan_tag_key='ScanId',
             filter_mask_key = filter_mask_key)
    min_num_glitch_map = {}
    for gt, gn in zip(glitch_thresholds, glitch_num_above):
        min_num_glitch_map[gt]=gn
    pipe.Add(timestreamflagging.glitchfinding.FlagGlitches,
             min_num_glitch_map = min_num_glitch_map,
             glitch_num_key = glitch_num_key,
             glitch_thresholds_key = glitch_thresholds_key,
             flag_key = flag_key)

    #large first derivative
    if not max_derivative_val is None:
        pipe.Add(todfilter.dftutils.AddMaxDerivative, ts_map_key=ts_key,
                 output_key=max_derivative_key, pre_smooth_scale=deriv_presmooth_scale)    
        pipe.Add(FlagBadG3MapValue, m_key=max_derivative_key, max_val=max_derivative_val, 
                 flag_reason="LargeFirstDeriv")

##########################
##########################
# Combined Pipesegments ##
##########################
##########################




@core.pipesegment
def FieldFlaggingPreKcmbConversion( pipe, ts_key, flag_key='Flags',
                                    cal_sources=['RCW38', 'MAT5A'],
                                    flag_on_delta_rfrac=False,
                                    flag_nan_offsets=False):
    pipe.Add(FlagInvalidData, ts_key = ts_key, flag_key = flag_key,
             flag_nan_offsets = flag_nan_offsets)
    pipe.Add(FlagUncalibratable, ts_key=ts_key, flag_key=flag_key)
    pipe.Add(FlagBadRfrac, ts_key=ts_key, latch_co = 0.3, overbias_co = 0.98,
             delta_rfrac_co = 0.035, flag_key = flag_key,
             flag_on_delta_rfrac = flag_on_delta_rfrac)

@core.pipesegment
def FieldFlaggingPostKcmbConversion( pipe, ts_key, flag_key='Flags', 
                                     variance_prefilter_poly_scale = -1,
                                     filter_mask_key = None,
                                     ps_mask_id = None, ps_pointing_key = None):

    pipe.Add(timestreamflagging.FlagNaNs,ts_key=ts_key,flag_key=flag_key,
             flag_id='PostCalibrationNaNs')
    pipe.Add(FlagBadG3MapValue, m_key='CalibratorResponseSN', min_val=20,
             flag_reason="BadCalSn")
    pipe.Add(FlagNoisyNonGaussianData, 
             ts_key=ts_key, flag_key=flag_key,
             variance_prefilter_poly_scale = variance_prefilter_poly_scale,
             glitch_thresholds = [20, 7], glitch_num_above = [1, 5],
             glitch_filter_scale = 11, filter_mask_key = filter_mask_key,
             ps_mask_id = ps_mask_id, ps_pointing_key = ps_pointing_key)
    pipe.Add(timestreamflagging.noiseflagging.FlagOscillatingChannels,
             ts_key=ts_key, flag_key=flag_key, threshold = 5)

@core.pipesegment
def FlagSaturatedBolosMars(pipe, flag_key, ts_key, thresh=1.15):
    '''
    Remove timestreams that saturate on Mars.
    '''
    pipe.Add(timestreamflagging.FlagSaturatedBolos, ts_key=ts_key,
             flag_key=flag_key, thresh=thresh)

