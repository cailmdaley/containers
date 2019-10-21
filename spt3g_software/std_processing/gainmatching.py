from spt3g import core
from spt3g.calibration.template_groups import get_template_groups
from spt3g.todfilter.dftutils import get_dft_of_ts_map
import numpy as np

@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def match_gains(frame, ts_key='RawTimestreams_I', flag_key=None,
                gain_match_key='GainMatchingCoeff', bolo_props=None, 
                freq_range=[0.1, 1.0], manual_gain_factor=None):
    if frame.type == core.G3FrameType.Scan:
        # collect the flag list if it exists
        if flag_key is not None:
            flags = frame[flag_key]
        else:
            flags = []

        # find good bolometer pairs to difference; will attempt to
        # difference all valid polarization pairs in the
        # timestreams
        pixel_tgroups = get_template_groups(bolo_props, per_band = True,
                                            per_pixel = True, per_wafer = True,
                                            include_keys=True)

        # take dfts of all bolometers
        ffts, freqs = get_dft_of_ts_map(frame[ts_key])
        freqs = freqs / core.G3Units.Hz
        f_range = (freqs>freq_range[0]) & (freqs<freq_range[1])
            
        frame[gain_match_key] = core.G3MapDouble()
        for group, bolos_in_pixel in pixel_tgroups.items():
            good_bolos_in_pixel = [bolo for bolo in bolos_in_pixel \
                                   if bolo not in flags and \
                                   bolo in frame[ts_key]]

            # calculate gain-matching coefficients
            if len(good_bolos_in_pixel) == 2:                    
                if manual_gain_factor == None:
                    fftx = ffts[good_bolos_in_pixel[0]]
                    ffty = ffts[good_bolos_in_pixel[1]]

                    fftx_inrange = fftx[f_range]
                    ffty_inrange = ffty[f_range]
                    XX = np.sum(np.abs(fftx_inrange)**2)
                    YY = np.sum(np.abs(ffty_inrange)**2)
                    ReXY = np.sum(np.real(np.conj(fftx_inrange)*ffty_inrange))
                    fX = np.sqrt( (XX + YY + 2*ReXY) / (2*XX + 2 * ReXY * np.sqrt(XX / YY)) )
                    fY = np.sqrt( (XX + YY + 2*ReXY) / (2*YY + 2 * ReXY * np.sqrt(YY / XX)) )

                else:
                    fX = manual_gain_factor[0]
                    fY = manual_gain_factor[1]

                frame[gain_match_key][good_bolos_in_pixel[0]] = fX
                frame[gain_match_key][good_bolos_in_pixel[1]] = fY
        
        return frame

def apply_gain_matching(frame, ts_key='RawTimestreams_I',
                        gain_match_key='GainMatchingCoeff',
                        out_key='GainMatchedTimestreams'):
    '''
    A pipeline segment that simply multiplies timestreams by gain-matching
    coefficients.
    
    Parameters
    ----------
    frame : G3Frame
        A frame.
    ts_key : str
        The key with timestreams on which to apply gain matching.
    gain_match_key : str
        The key containing the gain-matching coefficients.
    out_key : str
        The key to store the output, gain-matched timestreams.

    Returns
    -------
    None
    '''
    if not frame.type == core.G3FrameType.Scan:
        return

    out = core.G3TimestreamMap()
    for bolo, in_ts in frame[ts_key].iteritems():
        if bolo in frame[gain_match_key]:
            out_ts = core.G3Timestream(in_ts)
            out_ts *= frame[gain_match_key][bolo]
            out[bolo] = out_ts
    frame[out_key] = out


@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def difference_pairs(frame, ts_key='RawTimestreams_I',
                     gain_match_key='GainMatchingCoeff',
                     pair_diff_key='PairDiffTimestreams',
                     bolo_props=None):
    if frame.type == core.G3FrameType.Scan and \
       gain_match_key in frame:
        pixel_tgroups = get_template_groups(bolo_props, per_band = True,
                                            per_pixel = True, per_wafer = True,
                                            include_keys=True)
        frame[pair_diff_key] = core.G3TimestreamMap()
        bolo_with_gains = list(frame[gain_match_key].keys())

        for group, bolos_in_pixel in pixel_tgroups.items():
            good_bolos_in_pixel = [bolo for bolo in bolos_in_pixel \
                                   if bolo in bolo_with_gains]

            if len(good_bolos_in_pixel) == 2:
                if good_bolos_in_pixel[0] in frame[gain_match_key].keys() and \
                   good_bolos_in_pixel[1] in frame[gain_match_key].keys() and \
                   good_bolos_in_pixel[0] in frame[ts_key] and \
                   good_bolos_in_pixel[1] in frame[ts_key]:
                    coeff1 = frame[gain_match_key][good_bolos_in_pixel[0]]
                    coeff2 = frame[gain_match_key][good_bolos_in_pixel[1]]
                    ts1 = frame[ts_key][good_bolos_in_pixel[0]]
                    ts2 = frame[ts_key][good_bolos_in_pixel[1]]
                    frame[pair_diff_key][group] = coeff1*ts1 - coeff2*ts2
                
        return frame

@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def sum_pairs(frame, ts_key='RawTimestreams_I',
                     gain_match_key='GainMatchingCoeff',
                     pair_diff_key='PairDiffTimestreams',
                     bolo_props=None):
    if frame.type == core.G3FrameType.Scan and \
       gain_match_key in frame:
        pixel_tgroups = get_template_groups(bolo_props, per_band = True,
                                            per_pixel = True, per_wafer = True,
                                            include_keys=True)
        frame[pair_diff_key] = core.G3TimestreamMap()
        bolo_with_gains = list(frame[gain_match_key].keys())

        for group, bolos_in_pixel in pixel_tgroups.items():
            good_bolos_in_pixel = [bolo for bolo in bolos_in_pixel \
                                   if bolo in bolo_with_gains]

            if len(good_bolos_in_pixel) == 2:
                if good_bolos_in_pixel[0] in frame[gain_match_key].keys() and \
                   good_bolos_in_pixel[1] in frame[gain_match_key].keys() and \
                   good_bolos_in_pixel[0] in frame[ts_key] and \
                   good_bolos_in_pixel[1] in frame[ts_key]:
                    coeff1 = frame[gain_match_key][good_bolos_in_pixel[0]]
                    coeff2 = frame[gain_match_key][good_bolos_in_pixel[1]]
                    ts1 = frame[ts_key][good_bolos_in_pixel[0]]
                    ts2 = frame[ts_key][good_bolos_in_pixel[1]]
                    frame[pair_diff_key][group] = coeff1*ts1 + coeff2*ts2
                
        return frame
