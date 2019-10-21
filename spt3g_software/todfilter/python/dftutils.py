import numpy as np
import copy
from spt3g import core
from spt3g.todfilter import fft_filter_mem_friendly, MaskedPolyHpf, average_psd_bins
from spt3g.todfilter import rolling_mean_filter_ts, VarianceAdder, MadVarianceAdder, mean_smooth_ts
from spt3g.todfilter import get_psd_pybinding_grr

hz = core.G3Units.hz



@core.usefulfunc
def get_fft_padding_length( ts_len, other_factors = [1, 3, 5]):
    '''
    Returns the smallest length that is greater ts_len and can be expressed as:
    2**n * other_factor, where n is a natural number and other_factor is in other_factors

    FFT algorithms work faster on composites of small numbers.
    
    '''

    v = []
    ts_len = float(ts_len)
    for of in other_factors:
        if (ts_len < of):
            v.append(of)
        else:
            v.append(int(2**np.ceil(np.log(ts_len/of)/np.log(2)))*of)
    return min(v)


@core.indexmod
def AddFftPaddingLengthKey(frame, ts_map_key, padding_length_key):
    '''
    Adds a key for how long we should pad any filters we specify

    ts_map_key: timestream map so we know how long a timestream is

    padding_length_key: where to store the fft padding length
    '''

    if frame.type == core.G3FrameType.Scan:
        if padding_length_key in frame:
            return
        ts_m = frame[ts_map_key]
        if len(ts_m) == 0:
            return
        ts_len = len(ts_m[ts_m.keys()[0]])
        frame[padding_length_key] = get_fft_padding_length(ts_len)

@core.indexmod
def AddMaxDerivative(frame, ts_map_key, output_key, 
                     pre_smooth_scale = None):
    '''
    Adds the maximum value of the first derivative of the timestreams in the Scan.
    Before doing this it applies a mean smoothing to the data over pre_smooth_scale samples.

    The goal of this is to find latching detectors.


    pre_smooth_scale : Number of samples to smooth over, if None no smoothing is applied
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    ts_m = frame[ts_map_key]
    out_m = core.G3MapDouble()
    for k in ts_m.keys():
        ts = ts_m[k]
        if not pre_smooth_scale is None:
            ts = mean_smooth_ts(ts, pre_smooth_scale)
        out_m[k] = np.nan_to_num(np.max(np.abs(np.diff(ts))) * ts.sample_rate)
    frame[output_key] = out_m

def lowpass_func(f, cutoff_freq = 30):
    return np.exp( -1.0 * (f / cutoff_freq)**6.0 )

@core.indexmod
def LowPassFilterSpecifier(frame, 
                           cutoff_freq,
                           ts_map_key, #used to get the length
                           sample_rate_override_key = None,
                           output_filter_field = 'LowPassFilter', 
                           input_filter_field = None, 
                           already_specified_key = 'LowPassFilterSpecified', 
                           padding_key = 'FftPadding'):
    '''
    Specifies the low pass filter used by the fft filter.

    cutoff_freq, the cutoff frequency of the LPF.
    sample_rate_override_key: if not None, specifies a key in the frame to use as the sample rate
    ts_map_key: -> G3TimestreamMap: Used to find the length and sample rate of the filter.
    output_filter_field -> G3VectorDouble: where we stored the filter
    input_filter_field -> If specified multiplies the low pass filter by the input filter.
    already_specified_key:  In case we have low pass filter segments in a row and out of sloth
                            include the specifier multiple times we store this key the first time
                            we specify the filter and then if it is present later we 
                            don't specify the filter again
    padding_key: If specified, maps to an integer that is the amount of padding used on the fft filter.
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    if already_specified_key in frame and frame[already_specified_key]:
        return 
    frame[already_specified_key] = True
    ts_map = frame[ts_map_key]

    if padding_key in frame:
        filt_size = frame[padding_key]
    else:
        filt_size = len(ts_map[ts_map.keys()[0]])
    #if there is a previous filter grab it
    if not input_filter_field is None:
        filt = copy.copy(frame[input_filter_field])
    else:
        filt = np.zeros(filt_size) + 1

    if sample_rate_override_key != None:
        sample_rate = frame[sample_rate_override_key]
    else:
        sample_rate = ts_map[ts_map.keys()[0]].sample_rate

    #deprecated fix
    if np.isinf(sample_rate):
        sample_rate = frame["SampleRate"]

    freqs = np.fft.fftfreq(filt_size, 1.0/sample_rate)
    filt = core.G3VectorComplexDouble(filt * lowpass_func(freqs, cutoff_freq))
    frame[output_filter_field] = filt



@core.indexmod
def FftFilterCpp(frame, in_ts_key, filter_path, out_ts_key):
    '''
    Applies a fourier space filter to the data.
    
    in_ts_key  [-> G3TimestreamMap] : input timestream map key
    filter_path  [-> G3VectorComplexDouble] : The complex filter to apply to the timestream.  The frequency of each bin matches the
                                       numpy.fft.fftfreq definition of frequency
    out_ts_key [-> G3TimestreamMap] : output timestream map key
    
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    out_data = core.G3TimestreamMap()
    fft_filter_mem_friendly(frame[in_ts_key],
                            frame[filter_path],
                            out_data, False, None)
    frame[out_ts_key] = out_data


@core.indexmod
def FftFilterPython(frame, in_ts_key, filter_path, out_ts_key):
    '''
    Applies a fourier space filter to the data.
    
    in_ts_key  [-> G3TimestreamMap] : input timestream map key
    filter_path  [-> G3VectorDouble] : The complex filter to apply to the timestream.  The frequency of each bin matches the
                                       numpy.fft.fftfreq definition of frequency
    out_ts_key [-> G3TimestreamMap] : output timestream map key
    
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    out_map = fft_filter_ts_map_py(frame[in_ts_key], frame[filter_path])
    frame[out_ts_key] = out_map

@core.usefulfunc
def fft_filter_ts_map_py(ts_map, ts_filter, window_function = np.hamming):
    '''
    This is a function for applying a frequency based filter written in python

    This is mostly kept around for debugging purposes, you are encouraged to use:

    spt3g.todfilter.fft_filter_mem_friendly
    '''
    if len(ts_map) == 0:
        return {}, np.array([])
    import scipy.fftpack
    ts_len = len(ts_map[ts_map.keys()[0]])
    padding_width = len(ts_filter)
    out_map = core.G3TimestreamMap()
    for k in ts_map.keys():
        in_arr = ts_map[k] * window_function(ts_len)
        filtered_ts = scipy.fftpack.fft(in_arr, n=padding_width) * ts_filter
        out_map[k] = core.G3Timestream(np.real(scipy.fftpack.ifft(filtered_ts))[:ts_len]) / window_function(ts_len)
    return out_map

FftFilter = FftFilterCpp


@core.usefulfunc
def get_dft_of_ts_map(ts_map, 
                      window_function = np.hanning,
                      normalize_power = True,
                      pad = True, padding_other_factors = [1,3,5]):
    '''

    Applies a window function and performs the DFT on an input timestream.  
    If you are estimating a psd, use the get_psd_of_ts_map, it's faster.

    ts_map: the G3TimestreamMap input
    
    window_function: the window function to apply to the data.  (if you don't want a window_function pass 'lambda x: 1' to it.
    
    normalize_power:  Applies a normalization factor of 1.0/(np.mean(window_function(ts_len)**2.0)**0.5) to the output DFT to make it almost normalized

    pad:  FFTs tend to run faster if they are multiples of 2**n * small_number.  This will pad things to make it run faster

    padding_other_factors:  See get_fft_padding_length documentation

    returns: { k: dft(ts_map[k]), ... }, fft_freqs
    
    Uses numpy normalization of DFT.  Only handles forward transforms
    '''
    if len(ts_map) == 0:
        return {}, np.array([])

    import scipy.fftpack
    ts_len = len(ts_map[ts_map.keys()[0]])
    if pad:
        padding_width = get_fft_padding_length(ts_len, other_factors = padding_other_factors)
    else:
        padding_width = ts_len
    if normalize_power:
        norm_fac = 1.0/(np.mean(window_function(ts_len)**2.0)**0.5)
    else:
        norm_fac = 1.0
    out_d = {}

    sample_rate = ts_map[ts_map.keys()[0]].sample_rate
    for k in ts_map.keys():
        in_arr = ts_map[k] * window_function(ts_len)
        out_d[k] = scipy.fftpack.fft(in_arr, n=padding_width) * norm_fac
    return out_d, scipy.fftpack.fftfreq(padding_width, 1.0/sample_rate)


@core.usefulfunc
def get_psd_of_ts_map(ts_map, 
                      window_function = np.hanning,
                      pad = True, padding_other_factors = [1,3,5],
                      ):
    '''
    Returns the PSD for every detector in the timestream map.

    Output units are  Input Timestream Units ^2 / Hz


    Remember to incorporate the units system when interpretting the results.

    The psd units are (Input Units)^2 / Hz  you will want to multiply your result by 
    core.G3Units.Hz if you want to be able to plot it.  If this line makes no sense 
    please read the doc/units document.  

    Like, right now.
    
    '''

    if len(ts_map) == 0:
        return {}, np.array([])
        
    ts_map = copy.copy(ts_map)

    import scipy.fftpack
    ts_len = len(ts_map[ts_map.keys()[0]])
    if pad:
        padding_width = get_fft_padding_length(ts_len, other_factors = padding_other_factors)
    else:
        padding_width = ts_len
        
    out_d = {}    
    sample_rate = ts_map[ts_map.keys()[0]].sample_rate
    wf = core.G3VectorDouble(window_function(ts_len))

    out_size = padding_width//2 + 1

    #sets to divide by the bandwidth
    norm_fac = 1 / sample_rate 

    freqs = scipy.fftpack.fftfreq(padding_width, 1.0/sample_rate)[:out_size]
    freqs[-1] = np.abs(freqs[-1]) #so all the freqs are positive

    out_d = get_psd_pybinding_grr(ts_map, padding_width, wf, norm_fac)
    return out_d, freqs

@core.indexmod
def AddBinnedPsdInfo(frame, 
                     ts_map_key,
                     bins,
                     psd_label = 'BinAveragedPsds',
                     psd_lower_bound_label = 'BinAveragedPsdsLowerBound',
                     psd_upper_bound_label = 'BinAveragedPsdsUpperBound'
                 ):
    if frame.type != core.G3FrameType.Scan:
        return
    ts_map = frame[ts_map_key]
    psd_vals, freqs = get_psd_of_ts_map(ts_map)
    freq_sep = freqs[1]-freqs[0]

    lb, ub = zip(*bins)
    lower_bounds = core.G3VectorDouble(lb)
    upper_bounds = core.G3VectorDouble(ub)

    out_psd_av = core.G3MapVectorDouble()

    average_psd_bins(freq_sep, psd_vals, lower_bounds, upper_bounds, out_psd_av)
    
    frame[psd_label] = out_psd_av
    frame[psd_lower_bound_label] = lower_bounds
    frame[psd_upper_bound_label] = upper_bounds

    return frame

