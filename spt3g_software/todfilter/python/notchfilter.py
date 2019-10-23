import numpy as np
import scipy.stats as stats
import scipy.signal as signal
from operator import itemgetter
from spt3g import core
from spt3g.todfilter import fft_freqs
from spt3g.todfilter import notch_filter_lls
from spt3g.todfilter import make_empty_filter_mask
from spt3g.todfilter import poly_filter_ts_data_with_abscissa
from spt3g.todfilter import fft_filter_mem_friendly_w_multi_filters
from spt3g.todfilter.dftutils import get_psd_of_ts_map
from spt3g.calibration.template_groups import get_template_groups



'''

Potentially useful pipeline modules and segments contained in this script:

  1. GetMultiScanSpectra:
  
       This pipeline module concatenates the timestreams from
       multiple scan frames and returns the PSDs of the long timestreams.
       It can be useful when a user wants to examine PSDs with a
       higher frequency resolution than that can be provided by
       just one scan frame.
  
  
  2. NotchFilterFFTVersion
  
       This pipeline segment applies FFT-based notch filters to timestreams.
       It has the ability to concatenate multiple scan frames and then
       perform the notching on the long timestreams, which can be useful
       if the bandwidth a user wants to remove is much narrower than
       the sample frequency interval provided by just one scan.
       However, the concatenation process needs the data from
       scan frames corresponding to turnarounds as well,
       and this makes the pipeline segment not mock observation-friendly
       because mock observations ignore turnarounds
       (i.e. it is hard to simulate the effects of the filtering).
       In addition, concatenation increases memory usage.
       If doing the filtering on a scan-by-scan basis suffices,
       please use the next pipeline module.
  
  
  3. NotchFilterLLSVersion
  
       This pipeline module applies LLS-based notch filters to timestreams,
       which is basically the same as the masked high pass filter we use.
       It does not have the ability to concatenate multiple scan frames,
       and it is meant to be used on a scan-by-scan basis.
       Some advantages it has over the previous version are
       (1) it is mock observation-friendly
       (2) it allows masking point sources during filtering.

'''





class GetMultiScanSpectra(object):

    '''
    This module is meant to get PSDs from multiple scans.
    
    It first bundles multiple scans together.    
    With the bundled data,
    it performs a linear least squares filter
    that fits two sets of polynomials.
    One polynomial has an x axis of time,
    and the other has an x axis of Cosecant(elevation).
    
    After performing this filter it estimates the PSDs.
    The results are stored in two members of this class:
    self.out_psds and self.freqs.
    
    
    Init Arguments:
    ---------------
    
    timestream_map_key:
      where the timestream map lives
    
    elevation_key:
      where the telescope elevation data lives
    
    abscissa:
      defaults to 'cscel', if you change this to 'el',
      it will use elevation instead of csc(elevation) for the fit
    
    poly_order:
      the order of the polynomial using csc(el) as abscissa
    
    index_poly_order:
      the order the polynomial using time as abscissa
    
    n_scan_frames_to_cache:
      the number of scan frames to cache before performing a psd,
      when it gets to the end of processing,
      it will just take a psd with what is available.
    
    '''
    
    def __init__(self, timestream_map_key,
                 elevation_key,
                 abscissa = 'cscel',
                 poly_order = 7,
                 index_poly_order = 5,
                 n_scan_frames_to_cache = 40):
        
        valid_abscissa_opts = ['el', 'cscel']
        assert(abscissa in valid_abscissa_opts)
        self.abscissa = abscissa
        self.n_scan_frames_to_cache = n_scan_frames_to_cache
        self.cached_ts_maps = []
        self.cached_el = []
        
        self.timestream_map_key = timestream_map_key
        self.elevation_key = elevation_key
        
        self.poly_order = poly_order
        self.index_poly_order = index_poly_order
        
        self.out_psd = None
        self.out_freqs = None
    
    def __call__(self, frame):
        
        if not self.out_psd is None:
            return
        
        if frame.type != core.G3FrameType.Scan and \
           frame.type != core.G3FrameType.EndProcessing:
            return
        
        # we cache the frames
        if frame.type == core.G3FrameType.Scan:
            self.cached_ts_maps.append( frame[self.timestream_map_key] )
            self.cached_el.append( frame[self.elevation_key] )
        
        if (len(self.cached_ts_maps) >= self.n_scan_frames_to_cache or
            frame.type == core.G3FrameType.EndProcessing):
            if len(self.cached_el) == 0:
                return
            
            # if we've cached enough we concatenate the frames
            ts_map = core.G3TimestreamMap.concatenate(self.cached_ts_maps)
            self.cached_ts_maps = []
            elevation = core.G3Timestream.concatenate(self.cached_el)
            self.cached_el = []
            
            # we then filter to remove the sky signal
            # because it causes a huge 1/f^alpha knee
            if (self.abscissa == 'cscel'):
                elevation[:] = 1.0 / np.sin( np.abs(elevation) )
            ts_map_filt = poly_filter_ts_data_with_abscissa(
                              elevation, ts_map,
                              make_empty_filter_mask(ts_map),
                              self.poly_order,
                              self.index_poly_order);
            ts_map = None
            
            # actually get the psd
            psd, freqs = get_psd_of_ts_map(ts_map_filt)
            
            self.out_psd = psd
            self.out_freqs = freqs





@core.indexmod
class NotchFilterLLSVersion(object):

    '''
    This module applies notch filters to timestreams,
    and the filtering is based on linear least squares fitting.
    
    Sines and cosines at certain frequencies are fitted to and
    subtracted from timestreams. This method makes it possible to
    mask point sources while notching certain frequency components,
    which is just like the masked high-pass filter.
    
    By default, it is assumed that one of the input g3 files to
    the pipeline that contains this module stores calibration frames
    that specify which lines should be notched from which detectors.
    For an example of the format of the calibration frames, please see
    /home/weiquan/projects/notch_filter/
    line_libs/2018/line_library_for_2018_data.g3 on Amundsen/Scott.
    
    However, the user can also specify the lines to be notched
    in a dictionary and pass it to this module. For more information,
    please see the docstrings for the parameters below.
    
    An advantage of this version over the FFT version below is that
    the former does not require using the data from turnarounds,
    which the latter does in order to deal with
    potential oscillations introduced near the two ends of a timestream
    by the Fourier transform-based filtering.
    This is an advantage because mock observations ignore turnarounds,
    which makes it hard to simulate the effects of the FFT version.
    
    When using this module, please make sure that
    the input timestreams have already had a poly-filter applied
    (at least mean removed)!
    
    
    Parameters:
    -----------
    
    input_tsm_key: str
        The key in a scan frame that
        points to the timestreams to be filtered.
    
    output_tsm_key: str
        The key that will be used to store
        the filtered timestreams in the scan frame.
    
    use_masks: bool, optional
        Whether point source masks will be applied to
        the timestreams during the filtering. Default is True.
    
    masks_key: str, optional
        The key in a scan frame that stores the masks
        to be used during the filtering. Default is 'FilterMask'.
    
    specify_lines_as_args: bool, optional
        Whether the lines to be notched will be specified
        by passing dictionaries as arguments to this pipeline module
        instead of by including g3 files that contain information on lines
        in the input files to be supplied to the analysis script.
        Default is False. Please also see the example at the end of
        the docstring.
    
    line_dict_for_inc_az_scan: dict, optional
        If specify_lines_as_args is set to True,
        a dictionary specifying lines to be notched from
        timestreams of increasing-az scans needs to be supplied here.
        An example format for this type of dictionary is
        {"w201": [(1.39*core.G3Units.Hz, 1.44*core.G3Units.Hz), (...)]},
        where the keys can be individual bolometers' names or
        names of groups of bolometers, and the value for each key is a list
        of tuples, each of which has two elements that define
        the width of certain line.
        As for the keys,
        if notch_lines_by_group is True, then the possible options are
        those returned by the function
        calibration.template_groups.get_template_groups(),
        i.e. 90.0_w201, w181, and so on.
        If notch_lines_by_group is False, then the keys should be
        bolometer names.
        Default is None.
    
    line_dict_for_dec_az_scan: dict, optional
        If specify_lines_as_args is set to True,
        a dictionary specifying lines to be notched from
        timestreams of decreasing-as scans is supplied here.
        This can be the same dictionary as above
        if the lines do not depend on the scan direction.
        Default is None.
    
    notch_lines_by_group: bool, optional
        Whether certain set of lines will be removed from
        a group of timestreams collectivel, or
        different lines will be removed from different timestrems.
        If specify_lines_as_args is set to True,
        the value of this boolean needs to be consistent with the
        keys in the line dictionaries supplied.
        Otherwise (i.e. if g3 files storing line libraries are to be used),
        the value of this boolean will be overwritten by
        that stored in the line libraries.
        Default is True.
    
    by_band, by_wafer, by_board, by_squid: bools, optional
        Booleans specifying how bolometers should be grouped
        when lines are notched. These booleans will be supplied to
        calibration.template_groups.get_template_groups
        to obtain the appropriate groups. 
        If specify_lines_as_args is set to True,
        these booleans need to be consistent with the
        keys in the dictionaries supplied.
        Otherwise, the values of these booleans will be overwritten by
        those stored in the line libraries contained in the g3 files.
        Defaults are the only the by_wafer is set to True.
    
    show_log_messages: bool, optional
        Whether to show some log messages that
        report the progress of the filtering. Default is False.
    
    
    Examples:
    ---------
    
        line_dict = {"w181": [(1.39*core.G3Units.Hz, 1.45*core.G3Units.Hz)]}
        pipeline.Add(todfilter.notchfilter.NotchFilterLLSVersion,
                     input_tsm_key="PolyFilteredTimestreams",
                     output_tsm_key="NotchFilteredTimestreams",
                     specify_lines_as_args=True,
                     line_dict_for_inc_az_scan=line_dict,
                     line_dict_for_dec_az_scan=line_dict,
                     notch_lines_by_group=True,
                     by_wafer=True)
    
    '''
    
    def __init__(self,
                 input_tsm_key,
                 output_tsm_key,
                 use_masks=True,
                 masks_key="FilterMask",
                 specify_lines_as_args=False,
                 line_dict_for_inc_az_scan=None,
                 line_dict_for_dec_az_scan=None,
                 notch_lines_by_group=True,
                 by_band=False,
                 by_wafer=True,
                 by_board=False,
                 by_squid=False,
                 show_log_messages=False):
    
        self.input_tsm_key  = input_tsm_key
        self.output_tsm_key = output_tsm_key
        self.use_masks = use_masks
        self.masks_key = masks_key
        
        self.notch_lines_by_group = notch_lines_by_group
        self.boloproma = None
        self.wirma     = None
        self.by_band   = by_band
        self.by_wafer  = by_wafer
        self.by_board  = by_board
        self.by_squid  = by_squid
        
        self.specify_lines_as_args = specify_lines_as_args
        if self.specify_lines_as_args:
            self.line_dict_for_inc_az_scan = line_dict_for_inc_az_scan
            self.line_dict_for_dec_az_scan = line_dict_for_dec_az_scan
        else:
            self.line_libs_for_inc_az_scan = []
            self.line_libs_for_dec_az_scan = []
            self.line_lib_for_inc_az_scan  = None
            self.line_lib_for_dec_az_scan  = None
        
        if show_log_messages:
            core.set_log_level(core.G3LogLevel.LOG_INFO,
                               unit="notchfiltertimestreams")
    
    
    def figure_out_scan_direction(self, az_timestream):
    
        az_changes   = np.diff(az_timestream)
        most_changes = stats.mode(az_changes)[0][0]
        
        if most_changes > 0.0:
            scan_direction = "IncreasingAz"
        elif most_changes < 0.0:
            scan_direction = "DecreasingAz"
        else:
            scan_direction = "Strange"
        
        return scan_direction
    
    
    def line_lib_valid_for_this_obs(self, lib, obs_frame):
    
        if (lib["ObservationType"] == obs_frame["SourceName"])    and \
           (lib["ValidFrom"]       <= obs_frame["ObservationID"]) and \
           (lib["ValidThrough"]    >= obs_frame["ObservationID"]):
            is_valid = True
        else:
            is_valid = False
        
        return is_valid
    
    
    def get_freqs_to_notch_from_lib(self, group, lib,
                                    ts_sampling_rate, ts_n_samples):
    
        line_locs_key        = None
        line_left_edges_key  = None
        line_right_edges_key = None
        
        for key in lib.keys():
            if "LineLocations" in key:
                line_locs_key = key
            elif "LineLeftEdges" in key:
                line_left_edges_key = key
            elif "LineRightEdges" in key:
                line_right_edges_key = key
        
        if group not in lib[line_locs_key].keys():
            return None
        else:
            possible_freq_values = \
                np.asarray(fft_freqs(ts_sampling_rate, ts_n_samples))
            if ts_n_samples % 2 == 0:
                possible_freq_values[ts_n_samples//2] *= -1
            
            all_frequencies_to_notch = []
            for line_loc in lib[line_locs_key][group]:
                line_left_edge   = lib[line_left_edges_key][group][line_loc]
                line_right_edge  = lib[line_right_edges_key][group][line_loc]
                indices_to_notch = \
                    np.where((possible_freq_values >= line_left_edge ) &
                             (possible_freq_values <= line_right_edge))[0]
                freqs_to_notch = possible_freq_values[indices_to_notch]
                for freq_to_notch in freqs_to_notch:
                    all_frequencies_to_notch.append(freq_to_notch)
            
            return core.G3VectorDouble(all_frequencies_to_notch)
    
    
    def get_freqs_to_notch_from_dict(self, group, dct,
                                     ts_sampling_rate, ts_n_samples):
    
        if group not in dct.keys():
            return None
        else:
            possible_freq_values = \
                np.asarray(fft_freqs(ts_sampling_rate, ts_n_samples))
            if ts_n_samples % 2 == 0:
                possible_freq_values[ts_n_samples//2] *= -1
            
            all_frequencies_to_notch = []
            for tpl in dct[group]:
                line_left_edge  = tpl[0]
                line_right_edge = tpl[1]
                indices_to_notch = \
                    np.where((possible_freq_values >= line_left_edge ) &
                             (possible_freq_values <= line_right_edge))[0]
                freqs_to_notch = possible_freq_values[indices_to_notch]
                for freq_to_notch in freqs_to_notch:
                    all_frequencies_to_notch.append(freq_to_notch)
            
            return core.G3VectorDouble(all_frequencies_to_notch)
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Calibration:
            
            if "BolometerProperties" in frame:
                self.boloproma = frame["BolometerProperties"]
            elif "InformationOnLinesContained" in frame.keys():
                if   "LineLocationsIncreasingAzScan" in frame.keys():
                    self.line_libs_for_inc_az_scan.append(frame)
                elif "LineLocationsDecreasingAzScan" in frame.keys():
                    self.line_libs_for_dec_az_scan.append(frame)
            return
        
        
        if frame.type == core.G3FrameType.Wiring:
            
            self.wirma = frame["WiringMap"]
            return
        
        
        if frame.type == core.G3FrameType.Observation:
            
            if self.specify_lines_as_args:
                return
            
            # - Find the appropriate line libraries for this observation
            
            for candidate in self.line_libs_for_inc_az_scan:
                if self.line_lib_valid_for_this_obs(candidate, frame):
                    self.line_lib_for_inc_az_scan = candidate
                    break
            for candidate in self.line_libs_for_dec_az_scan:
                if self.line_lib_valid_for_this_obs(candidate, frame):
                    self.line_lib_for_dec_az_scan = candidate
                    break
            
            del self.line_libs_for_inc_az_scan
            del self.line_libs_for_dec_az_scan
            
            
            # - Make sure the libraries exist
            
            assert self.line_lib_for_inc_az_scan is not None, \
                   "Line library for increasing-az scans not found!"
            assert self.line_lib_for_dec_az_scan is not None, \
                   "Line library for decreasing-az scans not found!"
            
            core.log_info("\n", "A line library for increasing-az scans and",
                          "\n", "one for decreasing-az scans were both cached.",
                          "\n", unit="notchfiltertimestreams")
            return
        
        
        if frame.type == core.G3FrameType.Scan:
            
            core.log_info("\n", "A scan frame has arrived.",
                          "\n", "Necessary data will be gathered and",
                          "\n", "supplied to the filters.",
                          "\n", unit="notchfiltertimestreams")
            
            # - Retrive the correct library 
            
            az_timestream  = frame["RawBoresightAz"]
            scan_direction = self.figure_out_scan_direction(az_timestream)
            
            if scan_direction == "IncreasingAz":
                if self.specify_lines_as_args:
                    line_dict_for_this_scan = self.line_dict_for_inc_az_scan
                else:
                    line_lib_for_this_scan = self.line_lib_for_inc_az_scan
            elif scan_direction == "DecreasingAz":
                if self.specify_lines_as_args:
                    line_dict_for_this_scan = self.line_dict_for_dec_az_scan
                else:
                    line_lib_for_this_scan = self.line_lib_for_dec_az_scan
            else:
                core.log_info("\n", "Not clear what direction the scan is...",
                              "\n", "Notching will not occur for this frame.",
                              "\n", unit="notchfiltertimestreams")
                frame[self.output_tsm_name] = frame[self.input_tsm_name]
                return
            
            
            # - Decide on modes of notching
            
            input_tsm  = frame[self.input_tsm_key]
            
            if self.use_masks:
                core.log_info("\n", "Point source masks will be used",
                              "\n", "during the filtering.",
                              "\n", unit="notchfiltertimestreams")
                masks_for_notching = frame[self.masks_key]
            else:
                core.log_info("\n", "Point source masks will not be used",
                              "\n", "during the filtering.",
                              "\n", unit="notchfiltertimestreams")
                masks_for_notching = make_empty_filter_mask(input_tsm)
            
            if not self.specify_lines_as_args:
                self.notch_lines_by_group = line_lib_for_this_scan["by_group"]
            
            if self.notch_lines_by_group:
                if self.specify_lines_as_args:
                    groups_and_bolonames = \
                        get_template_groups(
                            self.boloproma, wiring_map=self.wirma,
                            per_band=self.by_band, per_wafer=self.by_wafer,
                            per_board=self.by_board, per_squid=self.by_squid,
                            include_keys=True, exclude_empty=True)
                else:
                    groups_and_bolonames = \
                        get_template_groups(
                            self.boloproma, wiring_map=self.wirma,
                            per_band=line_lib_for_this_scan["by_band"],
                            per_wafer=line_lib_for_this_scan["by_wafer"],
                            per_board=line_lib_for_this_scan["by_board"],
                            per_squid=line_lib_for_this_scan["by_squid"],
                            include_keys=True, exclude_empty=True)
            else:
                groups_and_bolonames = \
                    {key: key for key in input_tsm.keys()}
                ## The variable name in this case is not really applicable,
                ## but it is made to have the same format as above so that
                ## the later processing will be easier.
            
            
            # - Apply notching
            
            output_tsm = core.G3TimestreamMap()
            
            # -- Find out potential groups of bolos for notching
            
            if self.specify_lines_as_args:
                available_groups = line_dict_for_this_scan.keys()
            else:
                for key in line_lib_for_this_scan.keys():
                    if "LineLocations" in key:
                        available_groups = line_lib_for_this_scan[key].keys()
                        break
            
            for group in available_groups:
                
                # -- Ignore inapplicable groups
                
                if group not in groups_and_bolonames.keys():
                    core.log_info("\n", group, "seems like",
                                  "an invalid name!",
                                  "\n", "Line notching",
                                  "cannot be done for this group.",
                                  "\n", unit="notchfiltertimestreams")
                    continue
                
                available_bolonames = \
                    set(groups_and_bolonames[group]) & set(input_tsm.keys())
                
                if len(available_bolonames) == 0:
                    core.log_info("\n", group, "does not",
                                  "have any timestreams,",
                                  "\n", "Line notching",
                                  "cannot be done for this group.",
                                  "\n", unit="notchfiltertimestreams")
                    continue
                
                # -- Process applicable groups of timestreams
                
                input_tsm_this_group = core.G3TimestreamMap()
                for boloname in available_bolonames:
                    input_tsm_this_group[boloname] = input_tsm[boloname]
                
                if self.specify_lines_as_args:
                    frequencies_to_notch = \
                        self.get_freqs_to_notch_from_dict(
                            group, line_dict_for_this_scan,
                            input_tsm.sample_rate, input_tsm.n_samples)
                else:
                    frequencies_to_notch = \
                        self.get_freqs_to_notch_from_lib(
                            group, line_lib_for_this_scan,
                            input_tsm.sample_rate, input_tsm.n_samples)
                
                if len(frequencies_to_notch) == 0:
                    continue
                
                f_strs = [str(np.round(f/core.G3Units.Hz, 4))
                          for f in sorted(list(frequencies_to_notch))]
                n_f    = len(f_strs)
                f_strs = " ".join(f_strs)
                core.log_info("\n", "Line notching is happening for",
                              "\n", "this group of detectors:", group,
                              "\n", "  * The frequencies [Hz] to notch:",
                              "\n", f_strs,
                              "\n", "  * In total", n_f, "values.",
                              "\n", unit="notchfiltertimestreams")
                
                output_tsm_this_group = \
                    notch_filter_lls(
                        input_tsm_this_group,
                        masks_for_notching,
                        frequencies_to_notch)
                
                for boloname, ts in output_tsm_this_group.items():
                    output_tsm[boloname] = ts
            
            
            # - Copy timestreams that were not filtered
            
            missing_bolonames = set(input_tsm.keys()) - set(output_tsm.keys())
            for boloname in missing_bolonames:
                output_tsm[boloname] = input_tsm[boloname]
            
            
            # - Record the filtered timestream map
            
            frame[self.output_tsm_key] = output_tsm
            
            return





@core.pipesegment
def NotchFilterFFTVersion(
        pipeline,
        input_timestream_map_key,
        output_timestream_map_key,
        desired_frequency_resolution_in_mHz=20.0,
        use_premade_line_libraries=True,
        specify_lines_as_arg=False,
        line_dictionary=None,
        lines_from_median_psds_only=True,
        by_band=False,
        by_wafer=True,
        by_board=False,
        by_squid=False,
        elevation_key="OnlineBoresightEl",
        poly_order_in_elevation=0,
        poly_order_in_time_per_frame=3,
        ratio_of_threshold_to_model_psd=2.00,
        ratio_of_baseline_to_model_psd=1.20,
        power_law_fit_start_freq=0.06,
        power_law_fit_end_freq=0.60,
        white_noise_estimate_start_freq=4.5,
        white_noise_estimate_end_freq=7.5,
        low_frequency_cut=0.99,
        maximum_number_of_lines=50,
        drop_calframes_of_lines=True,
        show_log_messages=False,
        make_plots_of_lines_found=False,
        number_of_figures_to_generate=50,
        directory_to_save_figures="."):
    
    '''
    This is a pipeline segment that attemps to remove spectral lines
    from bolometer timestreams by using an FFT-based notch filter.
    
    When using this pipeline segment, please first note the following:
    
    This pipeline segment concatenates timestream maps
    from multiple scan frames, and, given the way the concatenation works,
    the place in which this pipeline segment can be put
    in a mapmaking pipeline is restricted:
    
    The concatenation assumes that
    the timestream maps are contiguous in time,
    so, scan frames corresponding to turnarounds must not be removed
    prior to this pipeline segment.
    
    The concatenation assumes that all the timestream maps 
    have the same set of keys, so the removal of flagged bolometers
    must not happen prior to this module in case
    different timestream maps had different bolometers removed.
    
    Now, please read on for more information about the code:
    
    The filter can be either a "dynamic" one or a "static" one.
    
    If the argument use_premade_line_libraries is set to False,
    the filter will be in the dynamic mode,
    which means that an automatic line-finding process will occur
    for each group of scan frames, and any found lines will be notched.
    
    If the argument use_premade_line_libraries is set to True,
    a premade library of lines needs to be supplied,
    and the fixed set of lines recorded in the library will be removed from
    the timestreams throughout an entire observation
    without the line-finding process being invoked.
    
    By default, it is assumed that the premade line library is contained
    in one of the input g3 files to the pipeline that contains this module.
    For an example of the format of the calibration frames, please see
    /home/weiquan/projects/notch_filter/
    line_libs/2018/line_library_for_2018_data.g3 on Amundsen/Scott.
    
    However, the user can also specify the lines to be notched
    in a dictionary and pass it to this module. For more information,
    please see the docstrings for the parameters below.
    
    This pipeline segment mainly consists of two modules:
    
    The first module is FindLinesFromPSDs,
    which attemps to find lines from the PSDs of timestreams
    from each group of scan frames and then records
    the information on the lines in a calibration frame.
    This calibration frame, along with the scan frames,
    will then be passed to the next module, RemoveLinesFromTimestreams.
    If the filter is in the static mode,
    the first module will simply gather the line library supplied and
    pass it and the scan frames to the second module
    without actually doing the line-finding.
    
    There are three main ingredients in the first module:
    
    The user can specify how much resolution in the frequency domain is desired.
    Based on the desired resolution, this module will concatenate timestreams
    from multiple scan frames so that the resultant timestreams become
    long enough to provide the needed resolution when FFT is performed.
    
    After concatenating the timestreams, this module will filter them
    by fitting/subtracting two polynomials to/from them.
    One polynomial is in the form of
    some timestream units (e.g. power) as a function of time,
    and the other polynomial is in the form of
    some timestream units (e.g. power) as a function of
    cosecant of the telescope elevation.
    The purpose of using the two polynomials is to
    reduce long-timescale drifts seen in timestreams and to
    reduce the effects due to the changes in the telescope elevation
    (and thus changes in atmospheric loading) that happened during those scans.
    
    Then, the module finds lines from the PSD of each timestream
    by constructing a model PSD and by regarding any PSD values
    that deviate from the model substantially as lines.
    The model PSD of each timestream is a combination of two terms.
    The first term is a power-law term,
    the amplitude and exponent of which are obtained by
    a linear least squares fitting to the PSD values in a low frequency range,
    and the second term is a constant,
    which is equal to the median PSD value in a high frequency range.
    
    As for the second module, RemoveLinesFromTimestreams,
    it finds the calibration frames generated or inserted by the first module,
    concatenates timestreams from each group of scans again,
    removes lines from the long timestreams through FFT,
    and finally splits the timestreams back into the original chunks.
    
    Explanation on each argument is below.
    
    
    Parameters:
    -----------
    
    input_timestream_map_key: str
        The key name of the timestream map that
        contains timestreams to be processed by this notch filter.
    
    output_timestream_map_key: str
        The key name of the timestream map that
        will contain the filtered timestreams.
    
    desired_frequency_resolution_in_mHz: float, optional
        As described above, this notch filter module uses FFT,
        and the user can specify how much resolution
        in the frequency space is needed
        by supplying certain number to this argument.
        The number should be in the units of milli Hertz.
        Then, the module FindLinesFromPSDs
        will concatenate timestreams from multiple scan frames
        so that the timestreams become long enough
        to result in the desired resolution in the frequency space
        when FFT is performed.
        Even in the case where one constant-speed scan
        can provide enough frequency resolution,
        that scan will still be concatenated with
        the previous turnaround and the next one so that
        potential edge effects that occur after the filtering
        can be contained within the turnarounds.
        Default is 20.0, which probably means one constant-speed scan suffices.
    
    use_premade_line_libraries: bool, optional
        If True, the static version of the filter will be used.
        If False, the dynamics version. Default is True.
    
    specify_lines_as_arg: bool, optional
        If use_premade_line_libraries is True, and this is also True,
        then the line libraries will not be supplied through
        input g3 files but through the argument below.
        Default is False. Please also see the example at the end of
        the docstring.
    
    line_dictionary: dict, optional
        A dictionary that contains the information on
        which lines to be notched. The format of this dictionary should be
        e.g. {"w201": [(1.39*core.G3Units.Hz, 1.44*core.G3Units.Hz), (...)]},
        where the keys can be individual bolometers' names or
        names of groups of bolometers, and the value for each key is a list
        of tuples, each of which has two elements that define
        the width of certain line.
        As for the keys, if lines_from_median_psds_only is True,
        then the possible options are those returned by the function
        calibration.template_groups.get_template_groups(),
        i.e. 90.0_w201, w181, and so on.
        If lines_from_median_psds_only is False,
        the keys should be bolometer names.
        Default is None.
    
    lines_from_median_psds_only: bool, optional
        If False, in the dynamics case, lines will be be found
        from every detector timestream's PSD, and, in the static case,
        the library of lines specifies lines from individual detectors.
        Then, different notch filters will be applied to
        different detectors' timestreams.
        If True, in the dynamics case, lines will be found
        from by-group median PSDs only, and, in the static case,
        the library of lines specifies lines from groups of detectors.
        Then, one filter will be applied to the timestreams of
        all the detectors in the same group.
        If use_premade_line_libraries is set to True,
        and specify_lines_as_arg is set to False,
        this boolean will be overwritten by the corresponding value stored
        in the line libraries in the input g3 files.
        If use_premade_line_libraries is set to True,
        and specify_lines_as_arg is set to True,
        the bool needs to be consistent with the line dictionary supplied
        (i.e. the keys must not be individual bolometers' names).
        Default is True.
    
    by_band, by_wafer, by_board, by_squid: bools, optional
        If lines_from_median_psds_only is set to True, 
        These booleans specify how timestreams should be grouped
        during the line-finding and line-removing processes
        (these booleans will be supplied to the function
        calibration.template_groups.get_template_groups()
        to get appropriate groups).
        If use_premade_line_libraries is set to True,
        but specify_lines_as_arg is set to False,
        these booleans will be overwritten by the values stored
        in the line libraries in the input g3 files.
        If use_premade_line_libraries is set to True,
        and specify_lines_as_arg is also set to True,
        these booleans need to be consistent with the line dictionary supplied
        (i.e. if only by_wafer is True, then the keys should be
        names like 'w201', 'w181', etc.).
        Defaults are such that only by_wafer is True.
    
    elevation_key: str, optional
        The key name in a scan frame that stores the timestream of
        the telescope's elevation. The elevation timestreams are used to
        poly-filter concatenated bolometer timestreams
        with the purpose of reducing any effects the changes in
        the telescope's elevation during turnarounds may have on timestreams.
        Default is 'OnlineBoresightEl'.
    
    poly_order_in_elevation: int, optional
        The order of the elevation domain polynomial that will be used
        for fitting/subtraction to/from the concatenated bolometer timestreams.
        Unlike poly_order_time_domain_per_frame,
        the number here is independent of the number of scan frames
        from which timestreams are concatenated.
        Default is 0 because the default desired_frequency_resolution_in_mHz
        corresponds to not concatenating more than one constant-speed scans
        (no change in the elevation).
    
    poly_order_in_time_per_frame: int, optional
        The order of the time domain polynomial that is used
        for fitting/subtraction to/from the concatenated bolometer timestreams.
        This number is a per-scan-frame number,
        so, if 7 is used here, for example,
        and if bolometer timestream maps from 5 scan frames are concatenated,
        then a time domain polynomial of order 35 will be filtered out.
        If this value is increased from the default value,
        then the default value for the argument power_law_fit_start_freq below
        should probably also be increased so that the power-law fit occurs
        in an appropriate frequency range. Default is 3.
    
    ratio_of_threshold_to_model_psd: float, optional
        The FindLinesFromPSDs module finds lines
        from the PSD of each timestream
        by constructing a model PSD and by regarding any PSD values
        that deviate from the model substantially as lines.
        The criterion of being substantial is specified by this argument.
        If the value is 2.00, for example, that means that,
        in certain frequency bin, if the actual PSD value is larger
        than the model PSD value by a factor 2.00,
        then the former is regarded as part of a line.
        2.00 seems to be not too bad if lines are to be found
        from by-wafer median PSDs only. If lines are to be found
        from every detector timestream's PSD,
        then a higher value like changed to 15 is recommended
        because individual bolometers' PSDs are noisier.
        For more information,
        please see the function find_lines_from_PSDs below.
        Default is 2.00.
    
    ratio_of_baseline_to_model_psd: float, optional
        When making an estimate of how wide each line is,
        FindLinesFromPSDs starts from the peak of a line and
        check the PSD values in the neighboring frequency bins
        to see at which bins the PSD values drop below some baseline values.
        The criterion for detemining the baseline values is
        specified by this argument.
        If it is 1.20, for example, the two frequency bins
        that have PSD values smaller than 1.20 times the model
        and that are closest to the bin at which the peak of a line occurs
        are used to define the left and right edges of a line.
        1.20 seems to be not too bad if lines are to be found
        from by-wafer median PSDs only. If lines are to be found from
        every detector timestream's PSD,
        then a higher value like 5 is recommended.
        Defauls is 1.20.
    
    power_law_fit_start_freq: float, optional
        Part of the line-finding process involves constructing a model PSD.
        The model PSD of each timestream is a combination of two terms.
        The first term is a power-law term,
        the amplitude and exponent of which are obtained by
        least-squares fitting to the PSD values in a low frequency range.
        This argument specifies the lower bound of this low frequency range
        in the units of Hz. Default is 0.06.
        If a lower value is desired,
        then the the values of the arguement poly_order_time_domain_per_frame
        above may need to be decreased so that the filtering does not
        remove too much of the '1/f' part of a PSD.
    
    power_law_fit_end_freq: float, optional
        Similar to the argument above,
        this one specifies the upper bound of the low frequency range
        in the units of Hz. Default is 0.60.
    
    white_noise_estimate_start_freq: float, optional
        The second term in the model PSD is a constant,
        which is equal to the median PSD in a high frequency range.
        This argument specifies the lower bound of this high frequency range
        in the units of Hz. Default value is 3.5.
    
    white_noise_estimate_end_freq: float, optional
        Similar to the argument right above,
        this one specifies the upper bound of the high frequency range
        in the units of Hz. The default value is 7.5.
    
    low_frequency_cut: float, optional
        Any lines found at frequencies below this value will be ignored.
        The number should be expressed in units of Hz.
        This argument exists because the code can sometimes regard peaks
        shown in a noisy low frequency part of a PSD as lines,
        but it's not clear whether those peaks are really lines.
        Default is 0.99.
    
    maximum_number_of_lines: int, optional
        The code stops finding lines once this many lines have been found.
        This is to prevent a situation where for some reason
        the constructed model is very bad,
        and a lot of outlier PSD values are found. Default is 50.
    
    drop_calframes_of_lines: bool, optional
        If set to False, the calibration frames inserted by
        FindLinesFromPSDs will not be dropped
        at the end of this pipeline segment.
        Keeping these calibration frames may be useful for the purpose of,
        for example, checking some statistics on the lines found.
        Default is True.
    
    show_log_messages: bool, optional
        If set to True, various log messages indicating
        the progress of the filtering will be shown. Default is False.
    
    make_plots_of_lines_found: bool, optional
        If set to True, figures that show the lines found will be generated,
        which may serve as a useful check to see
        how well the line-finding process worked. Default is False.
    
    number_of_figures_to_generate: int, optional
        If figures of found lines are to be made, this argument specifies
        the approximate number of figures that will be made.
        For example, if 10 figures are to be made,
        but lines were found from 10000 PSDs, then figures will be made for
        only 1 in 100 PSDs and their lines. Default is 50.
    
    directory_to_save_figures: str, optional
        The location where the figures will be saves.
        The file names will be automatically chosen by the code.
        Default is '.'.
    
    
    Examples:
    ---------
    
    line_dict = {"w181": [(1.39*core.G3Units.Hz, 1.45*core.G3Units.Hz)]}
    pipeline.Add(todfilter.notchfilter.NotchFilterFFTVersion,
    input_timestream_map_key="RawTimestreams_I",
    output_timestream_map_key="NotchedRawTimestreams_I",
    use_premade_line_libraries=False,
    specify_lines_arg=True,
    line_dictionary=line_dict,
    lines_from_median_psds_only=True,
    by_wafer=True)
    
    '''
    
    if show_log_messages:
        core.set_log_level(core.G3LogLevel.LOG_INFO,
                           unit="notchfiltertimestreams")
    
    pipeline.Add(FindLinesFromPSDs,
                 input_timestream_map_key=input_timestream_map_key,
                 desired_freq_resol_in_mHz=desired_frequency_resolution_in_mHz,
                 use_premade_line_libraries=use_premade_line_libraries,
                 specify_lines_as_arg=specify_lines_as_arg,
                 line_dictionary=line_dictionary,
                 lines_from_median_psds_only=lines_from_median_psds_only,
                 by_band=by_band, by_wafer=by_wafer,
                 by_board=by_board, by_squid=by_squid,
                 elevation_key=elevation_key,
                 poly_order_in_elevation=poly_order_in_elevation,
                 poly_order_in_time_per_frame=poly_order_in_time_per_frame,
                 ratio_of_threshold_to_model=ratio_of_threshold_to_model_psd,
                 ratio_of_baseline_to_model=ratio_of_baseline_to_model_psd,
                 power_law_fit_start=power_law_fit_start_freq,
                 power_law_fit_end=power_law_fit_end_freq,
                 white_noise_estimate_start=white_noise_estimate_start_freq,
                 white_noise_estimate_end=white_noise_estimate_end_freq,
                 low_frequency_cut=low_frequency_cut,
                 maximum_number_of_lines=maximum_number_of_lines,
                 make_plots_of_lines_found=make_plots_of_lines_found,
                 number_of_figures_to_generate=number_of_figures_to_generate,
                 directory_to_save_figures=directory_to_save_figures)
    
    pipeline.Add(RemoveLinesFromTimestreams,
                 input_timestream_map_key=input_timestream_map_key,
                 output_timestream_map_key=output_timestream_map_key)
    
    if drop_calframes_of_lines:
        pipeline.Add(lambda frame: "InformationOnLinesContained" not in frame)





''' The individual pieces called by NotchFilterFFTVersion are defined below. '''



class FindLinesFromPSDs(object):

    '''
    This module attemps to find spectral lines from timestreams
    (assuming premade line libraries are not to be supplied).
    
    It first keeps caching timestream maps from consecutive scan frames
    until both of the following conditions are met:
    
    The timestreams from these scan frames cover a long enough
    time interval so that the desired frequency resolution will be achieved
    when FFT is performed on the concatenated timestreams.
    
    The first and last scan frames corrrespond to turnarounds.
    This is related to the fact that line-notching process used
    in the RemoveLinesFromTimestreams module
    consists of multipling a timestream by a window function,
    performing FFT, setting some Fourier coefficients to zero,
    performing inverse FFT, and dividing the timestream by the window function,
    and this process seems to sometimes introduce very noticeable oscillations
    near both ends of a timestream. By making sure that
    the first and last scan frames of a group of concatenated scan frames
    correspond to turnarounds, we can confine those oscillations
    within the turnarounds, which is better than leaving them
    near an end of a timestreams from a constant-speed scan
    because the data from turnarounds will eventually not be used
    in the map-making anyways.
    
    Then, the concatenated timestream maps will be filtered by
    linear least squares-based poly-filters.
    For more information on the filtering, please see the docstrings
    in the piepeline segment NotchFilterFFTVersion above.
    
    After that, lines will be found from the concatenated timestreams,
    and the innformation will be recorded in a calibration frame.
    For more informatin on the line-finding procedure and this calbration frame,
    please see the docstings in the functions find_lines_from_PSDs below.
    
    This calibration frame, along with the scan frames from which
    timestreams were cached, will be passed to the next processing step.
    For example, if this module caches timestream maps from 3 scan frames
    and finds lines from the concatenated timestream maps,
    the pipeline module that comes right after this one will
    see 3 scan frames followed by a calibration frame.
    
    If use_premade_line_libraries is set to True, then
    the filtering and line-finding processes will not occur, but
    the concatenation process will still be done so that the line libraries
    will be inserted at the correct locations among scan frames.
    
    
    Paramters:
    ----------
    
    For explanations on the arguments, please see the docstings
    in the piepeline segment NotchFilterFFTVersion above.
    
    '''
    
    def __init__(self,
                 input_timestream_map_key,
                 desired_freq_resol_in_mHz=20.0,
                 use_premade_line_libraries=False,
                 specify_lines_as_arg=False,
                 line_dictionary=None,
                 lines_from_median_psds_only=True,
                 by_band=True,
                 by_wafer=True,
                 by_board=False,
                 by_squid=False,
                 elevation_key="OnlineBoresightEl",
                 poly_order_in_elevation=0,
                 poly_order_in_time_per_frame=3,
                 ratio_of_threshold_to_model=2.00,
                 ratio_of_baseline_to_model=1.20,
                 power_law_fit_start=0.06,
                 power_law_fit_end=0.60,
                 white_noise_estimate_start=4.5,
                 white_noise_estimate_end=7.5,
                 low_frequency_cut=0.99,
                 maximum_number_of_lines=50,
                 make_plots_of_lines_found=False,
                 number_of_figures_to_generate=50,
                 directory_to_save_figures="."):
        
        # These variables are useful for several purposes
        # such as grouping timestreams, checking observation type,
        # and so on.
        
        self.bolometer_properties = None
        self.wiring_map           = None
        self.observation_frame    = None
        
        # These variables are related to keeping track of
        # what is inside various caches
        
        self.counter                      = 0
        self.cached_timestream_maps       = []
        self.cached_elevation_timestreams = []
        self.cached_scan_frame_numbers    = []
        self.cached_scan_frame_types      = []
        self.n_data_pts_each_timestream   = []
        self.sampling_rates               = []
        self.total_n_data_pts             = 0
        self.need_to_clear_cache          = False
        
        # These variables are related to deciding
        # whether more frames need to be cached or not
        
        self.current_freq_resol_in_mHz = np.inf
        self.desired_freq_resol_in_mHz = desired_freq_resol_in_mHz
        
        # These variables are related to the process of
        # concatenating and filtering timestreams
        
        self.timestream_map_key          = input_timestream_map_key
        self.elevation_data_key          = elevation_key
        self.poly_order_t_domain_per_fr  = poly_order_in_time_per_frame
        self.poly_order_elevation_domain = poly_order_in_elevation
        
        # These variables specify the parameters
        # used in the line-finding process
        
        self.from_group_median_psd_only  = lines_from_median_psds_only
        self.by_band                     = by_band
        self.by_wafer                    = by_wafer
        self.by_board                    = by_board
        self.by_squid                    = by_squid
        self.ratio_of_threshold_to_model = ratio_of_threshold_to_model
        self.ratio_of_baseline_to_model  = ratio_of_baseline_to_model
        self.power_law_fit_start_freq    = power_law_fit_start
        self.power_law_fit_end_freq      = power_law_fit_end
        self.white_noise_estimate_start  = white_noise_estimate_start
        self.white_noise_estimate_end    = white_noise_estimate_end
        self.low_frequency_cut           = low_frequency_cut
        self.maximum_number_of_lines     = maximum_number_of_lines
        
        # These variables are for making figures as a debugging tool
        
        self.make_plots_of_lines_found   = make_plots_of_lines_found
        self.n_figures_to_generate       = number_of_figures_to_generate
        self.directory_to_save_figures   = directory_to_save_figures
        
        # These variables are related to the case where
        # the module will just insert line libraries at
        # appropriate locations and not do the line-finding
        
        self.just_insert_libs = use_premade_line_libraries
        if self.just_insert_libs:
            self.specify_lines_as_arg = specify_lines_as_arg
            if self.specify_lines_as_arg:
                self.line_dictionary      = line_dictionary
                self.line_frame_from_dict = None
            else:
                self.candidate_libs = []
                self.correct_lib    = None            
            self.just_return_concatenated = True
            # If line libraries already exist,
            # there is no need to filter timestreams
            # inside the function that concatenates and filters timestreams.
        else:
            self.just_return_concatenated = False
    
    
    def convert_line_dict_to_calframe(self, line_dictionary):
    
        calframe = core.G3Frame(core.G3FrameType.Calibration)
        
        calframe["InformationOnLinesContained"] = True
        calframe["ObservationType"] = self.observation_frame["SourceName"]
        calframe["ValidFrom"]       = self.observation_frame["ObservationID"]
        calframe["ValidThrough"]    = self.observation_frame["ObservationID"]
        
        calframe["LineLocations"]  = core.G3MapVectorString()
        calframe["LineLeftEdges"]  = core.G3MapMapDouble()
        calframe["LineRightEdges"] = core.G3MapMapDouble()
        
        calframe["by_group"] = self.from_group_median_psd_only
        calframe["by_band"]  = self.by_band
        calframe["by_wafer"] = self.by_wafer
        calframe["by_board"] = self.by_board
        calframe["by_squid"] = self.by_squid
        
        for key, bandwidths in line_dictionary.items():
            calframe["LineLocations"][key]  = core.G3VectorString()
            calframe["LineLeftEdges"][key]  = core.G3MapDouble()
            calframe["LineRightEdges"][key] = core.G3MapDouble()
            for tpl in bandwidths:
                left   = tpl[0]
                right  = tpl[1]
                center = str(np.mean([left, right]))
                calframe["LineLocations"][key].append(center)
                calframe["LineLeftEdges"][key][center]  = left
                calframe["LineRightEdges"][key][center] = right
        
        return calframe
    
    
    def process_timestreams(self):
    
        core.log_info(
            "\n", "In the FindLinesFromPSDs module,",
            "\n", "enough timestream maps have been cached",
            "\n", "because the desired frequency resolution is",
            np.round(self.desired_freq_resol_in_mHz, 2), "mHz",
            "\n", "and because the first and last maps",
            "correspond to turnarounds.",
            "\n", "Processing data from scan frames",
            str(self.cached_scan_frame_numbers[0]), "to",
            str(self.cached_scan_frame_numbers[-1])+"!",
            "\n", unit="notchfiltertimestreams")
        
        concatenated_and_filtered_timestream_maps = \
            get_concatenated_and_filtered_timestreams(
                self.cached_timestream_maps,
                self.cached_elevation_timestreams,
                self.timestream_map_key,
                self.elevation_data_key,
                just_return_concatenated=\
                    self.just_return_concatenated,
                poly_order_time_domain_per_frame=\
                    self.poly_order_t_domain_per_fr,
                poly_order_elevation_domain=\
                    self.poly_order_elevation_domain)
        
        if concatenated_and_filtered_timestream_maps is None:
            core.log_info(
                "\n", "In the FindLinesFromPSDs module,",
                "\n", "There was an error in concatenating",
                "these timestream maps,",
                "\n", "and line notching will not occur for them.",
                "\n", unit="notchfiltertimestreams")
            calframe_for_lines = core.G3Frame(core.G3FrameType.Calibration)
            calframe_for_lines["SkipLineNotching"] = True
        
        elif self.just_insert_libs:
            core.log_info(
                "\n", "In the FindLinesFromPSDs module,",
                "\n", "a line library is inserted",
                "\n", "after the previous scan frame.",
                "\n", unit="notchfiltertimestreams")
            if self.specify_lines_as_arg:
                calframe_for_lines = self.line_frame_from_dict
            else:
                calframe_for_lines = self.correct_lib
        
        else:
            calframe_for_lines = \
                find_lines_from_PSDs(
                    concatenated_and_filtered_timestream_maps,
                    self.n_data_pts_each_timestream,
                    self.bolometer_properties,
                    self.wiring_map,
                    from_by_group_median_PSD_only=\
                        self.from_group_median_psd_only,
                    by_band=self.by_band, by_wafer=self.by_wafer,
                    by_board=self.by_board, by_squid=self.by_squid,
                    factor_for_thresholds=self.ratio_of_threshold_to_model,
                    factor_for_baseline=self.ratio_of_baseline_to_model,
                    power_law_fit_starting_frequency=\
                        self.power_law_fit_start_freq,
                    power_law_fit_ending_frequency=\
                        self.power_law_fit_end_freq,
                    white_noise_estimate_starting_frequency=\
                        self.white_noise_estimate_start,
                    white_noise_estimate_ending_frequency=\
                        self.white_noise_estimate_end,
                    low_frequency_cut=self.low_frequency_cut,
                    maximum_number_of_lines=self.maximum_number_of_lines,
                    make_plots_of_lines_found=self.make_plots_of_lines_found,
                    n_figures_to_generate=self.n_figures_to_generate,
                    directory_to_save_figures=self.directory_to_save_figures,
                    observation_frame=self.observation_frame,
                    scan_frame_numbers=self.cached_scan_frame_numbers)
        
            calframe_for_lines["InformationOnLinesContained"] = True
            obs_id = self.observation_frame["ObservationID"]
            calframe_for_lines["ValidFrom"]    = obs_id
            calframe_for_lines["ValidThrough"] = obs_id
            calframe_for_lines["TimestreamMapName"] = self.timestream_map_key
        
        return calframe_for_lines
    
    
    def __call__(self, frame):
        
        if frame.type == core.G3FrameType.Calibration:
            if "NominalBolometerProperties" in frame:
                self.bolometer_properties = frame["NominalBolometerProperties"]
                return
            if "BolometerProperties" in frame:
                self.bolometer_properties = frame["BolometerProperties"]
                return
            if "InformationOnLinesContained" in frame:
                if self.just_insert_libs:
                    if not self.specify_lines_as_arg:
                        self.candidate_libs.append(frame)
                return
        
        if frame.type == core.G3FrameType.Observation:
            self.observation_frame = frame
            if self.just_insert_libs:
                if self.specify_lines_as_arg:
                    self.line_frame_from_dict = \
                        self.convert_line_dict_to_calframe(self.line_dictionary)
                    del self.line_dictionary
                else:
                    for lib in self.candidate_libs:
                        if (lib["ObservationType"] == frame["SourceName"]) and \
                           (lib["ValidFrom"]    <= frame["ObservationID"]) and \
                           (lib["ValidThrough"] >= frame["ObservationID"]):
                            self.correct_lib = lib
                            del self.candidate_libs
                            break
                    assert self.correct_lib is not None, \
                           "An appropriate line library to use was not found!"
            return
        
        if frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame["WiringMap"]
            return
        
        if (frame.type == core.G3FrameType.Scan) or \
           (frame.type == core.G3FrameType.EndProcessing):
            self.counter += 1
            
            if self.need_to_clear_cache:
                core.log_info(
                    "\n", "In the FindLinesFromPSDs module",
                    "\n", "clearing the cache because it's full!",
                    "\n", "(The data from the last scan frame will be kept",
                    "\n", " because it corresponds to a turnaround and",
                    "\n", " will be used as a buffer for the beginning of",
                    "\n", " the next set of timestreams.)",
                    "\n", unit="notchfiltertimestreams")
                
                self.cached_timestream_maps = \
                    [self.cached_timestream_maps[-1]]
                self.cached_elevation_timestreams = \
                    [self.cached_elevation_timestreams[-1]]
                self.cached_scan_frame_numbers = \
                    [self.cached_scan_frame_numbers[-1]]
                self.cached_scan_frame_types = \
                    [self.cached_scan_frame_types[-1]]
                self.n_data_pts_each_timestream = \
                    [self.n_data_pts_each_timestream[-1]]
                self.sampling_rates = \
                    [self.sampling_rates[-1]]
                
                # Using a for loop to loop through these caches and update them
                # doesn't seem to work because the caches seem to go back to
                # the original state after we exit from the loop...
                
                
                self.total_n_data_pts = \
                    self.n_data_pts_each_timestream[0]
                self.current_freq_resol_in_mHz = \
                    1e3 / (self.total_n_data_pts / self.sampling_rates[0])
                self.need_to_clear_cache = False
            
            
            # - Now, 3 situations are possible here...
            
            ## - === SITUATION 1 ===
            ## - If there is not enough number of frames in the cache yet,
            ## - another frame is cached 
            ## - (assuming EndProcessing has not arrived).
            ## - Then, a check is made to see if more is needed.
            
            if frame.type == core.G3FrameType.Scan:
                
                self.cached_timestream_maps.\
                    append(frame[self.timestream_map_key])
                self.cached_elevation_timestreams.\
                    append(frame[self.elevation_data_key])
                self.cached_scan_frame_numbers.\
                    append(self.counter)
                self.n_data_pts_each_timestream.\
                    append(frame[self.timestream_map_key].n_samples)
                self.sampling_rates.\
                    append(frame[self.timestream_map_key].sample_rate)
                if "Turnaround" in frame:
                    self.cached_scan_frame_types.append("Turnaround")
                else:
                    self.cached_scan_frame_types.append("Scan")
                
                self.total_n_data_pts = \
                    np.sum(self.n_data_pts_each_timestream)
                avg_sampling_rate = \
                    np.mean(self.sampling_rates) / core.G3Units.Hz
                self.current_freq_resol_in_mHz = \
                    1000.0 / (self.total_n_data_pts / avg_sampling_rate)
                
                core.log_info(
                    "\n", "In the FindLinesFromPSDs module,",
                    "\n", "timestream map from scan frame",
                    str(self.counter).rjust(2),
                    "\n", "is being added to the cache...",
                    "\n", "Number of timestream maps cached so far:",
                    len(self.cached_timestream_maps),
                    "\n", "These are the corresponding scan frames:",
                    self.cached_scan_frame_numbers,
                    "\n", "List of length of each cached timestream:",
                    self.n_data_pts_each_timestream,
                    "\n", "list of type of each scan frame:",
                    self.cached_scan_frame_types,
                    "\n", "Sum of the length of each cached timestream:",
                    self.total_n_data_pts,
                    "\n", "If the timestreams were concatenated now,",
                    "\n", "the frequency resolution would be",
                    np.round(self.current_freq_resol_in_mHz, 2),
                    "milli Hertz.",
                    "\n", unit="notchfiltertimestreams")
                
                
                ### - If the cache still does not have enough data,
                ### - more data will be added.
                ### - Otherwise, a line-finding process starts.
                
                resolution_not_good  = \
                    self.current_freq_resol_in_mHz > \
                    self.desired_freq_resol_in_mHz
                no_buffer_at_one_end = \
                    (self.cached_scan_frame_types[0]  != "Turnaround") or \
                    (self.cached_scan_frame_types[-1] != "Turnaround")
                
                if resolution_not_good or no_buffer_at_one_end:
                    return [frame]
                else:
                    calframe_for_lines = self.process_timestreams()
                    self.need_to_clear_cache = True
                    return [frame, calframe_for_lines]
            
            
            ## - === SITUATION 2 ===
            ## - If the Endprocessing has arrived,
            ## - and there is only one scan frame in the cache
            ## - (and the observation has more than one scan frame),
            ## - that is just the leftover from the last line-finding process,
            ## - so there is nothing to do.
            
            elif (frame.type == core.G3FrameType.EndProcessing) and \
                 (len(self.cached_timestream_maps) <= 1)        and \
                 (self.counter > 2):
                
                core.log_info(
                    "\n", "In the FindLinesFromPSDs module,",
                    "\n", "EndProcessing frame has arrived.",
                    "\n", "There is one timestream map in the cache,",
                    "\n", "but it probably has slready been processed,",
                    "\n", "so, there is nothing to do!",
                    "\n", unit="notchfiltertimestreams")
                return [frame]
            
            
            ## - === SITUATION 3 ===
            ## - If the Endprocessing has arrived,
            ## - but the cache contains data that have not been processed yet,
            ## - then the line-finding process will occur
            ## - regardless of whether the cache has enough number of frames
            ## - to meet the two conditions.
            
            elif (frame.type == core.G3FrameType.EndProcessing) and \
                 (len(self.cached_timestream_maps) > 1):
                
                core.log_info(
                    "\n", "In the FindLinesFromPSDs module,",
                    "\n", "EndProcessing frame has arrived,",
                    "\n", "and there are no more scan frames,",
                    "\n", "so the data in the cache are processed!",
                    "\n", "Processing data from scan frames",
                    str(self.cached_scan_frame_numbers[0]), "to",
                    str(self.cached_scan_frame_numbers[-1])+"!",
                    "\n", unit="notchfiltertimestreams")
                calframe_for_lines = self.process_timestreams()
                return [calframe_for_lines, frame]





class RemoveLinesFromTimestreams(object):

    '''
    This module attemps to remove the spectral lines
    found by the previous module, FindLinesFromPSDs,
    and thus it needs to be used together with it.
    
    The module FindLinesFromPSDs concatenates timestream maps
    from multiple scan frames, finds lines from the long timestreams,
    and records the information on the lines in a calibration frame.
    Then, it passes this calibration frame followed by
    the corresponding scan frames to the next processing step in the pipeline.
    
    This module, RemoveLinesFromTimestreams, first caches
    all the scan frames it receives until the aforementioned
    calibration frame arrives. Then, it concatenates the timestream maps
    from the cached scan frames and uses the infomration on lines
    to apply notch filters.
    
    After that, it splits each filtered timestream into the original chunks
    and adds these chucks to the corresponding scan frames.
    
    
    Parameters:
    -----------
    
    input_timestream_map_key: str
    output_timestream_map_key: str
    
    '''
    
    def __init__(self, 
                 input_timestream_map_key,
                 output_timestream_map_key):
        
        self.original_timestream_map_key = input_timestream_map_key
        self.filtered_timestream_map_key = output_timestream_map_key
        self.cached_scan_frames          = []
        self.cached_scan_frame_numbers   = []
        self.n_data_pts_each_timestream  = []
        self.scan_frame_counter          = 0
        self.need_to_clear_cache         = False
        self.boloproma = None
        self.wirma     = None
    
    
    def __call__(self, frame):
    
        if self.need_to_clear_cache:
            core.log_info(
                "\n", "In the RemoveLinesFromTimestreams module,",
                "\n", "the cache needs to be cleaned because",
                "\n", "the processing was done for the current set of frames",
                "\n", "and a new set of scan frames is about to arrive.",
                "\n", unit="notchfiltertimestreams")
            
            self.cached_scan_frames = \
                [self.cached_scan_frames[-1]]
            self.cached_scan_frame_numbers = \
                [self.cached_scan_frame_numbers[-1]]
            self.n_data_pts_each_timestream = \
                [self.n_data_pts_each_timestream[-1]]
            
            self.need_to_clear_cache = False
        
        
        if frame.type == core.G3FrameType.Scan:
            self.scan_frame_counter += 1
            self.cached_scan_frame_numbers.append(self.scan_frame_counter)
            self.cached_scan_frames.append(frame)
            self.n_data_pts_each_timestream.\
                append(frame[self.original_timestream_map_key].n_samples)
            
            core.log_info(
                "\n", "In the RemoveLinesFromTimestreams module,",
                "\n", "scan frame", self.scan_frame_counter,
                "is added to the cache.", 
                "\n", "Number of scan frames in cache:", 
                len(self.cached_scan_frames),
                "\n", "These are the corresponding scan frames:",
                self.cached_scan_frame_numbers,
                "\n", "Length of the timestream from each scan frame:",
                self.n_data_pts_each_timestream,
                "\n", unit="notchfiltertimestreams")
            
            return []
        
        
        if frame.type == core.G3FrameType.Wiring:
            self.wirma = frame["WiringMap"]
            return frame
        
        
        if frame.type == core.G3FrameType.Calibration:
            if "NominalBolometerProperties" in frame:
                self.boloproma = frame["NominalBolometerProperties"]
                return frame
            
            if "BolometerProperties" in frame:
                self.boloproma = frame["BolometerProperties"]
                return frame
            
            if "SkipLineNotching" in frame:
                core.log_info(
                    "\n", "In the RemoveLinesFromTimestreams module,",
                    "\n", "a frame having information on lines is received,",
                    "\n", "but it looks that line-finding was not done,",
                    "\n", "so there is nothing to notch.",
                    "\n", "Each scan frame will get a new key",
                    self.filtered_timestream_map_key+",",
                    "\n", "but the value stored is the same as",
                    "that stored in", self.original_timestream_map_key+".",
                    "\n", unit="notchfiltertimestreams")
                
                for i in range(len(self.cached_scan_frames)):
                    if self.filtered_timestream_map_key      \
                    not in self.cached_scan_frames[i].keys():
                        self.cached_scan_frames[i]           \
                        [self.filtered_timestream_map_key] = \
                            self.cached_scan_frames[i]       \
                            [self.original_timestream_map_key]
                
                self.need_to_clear_cache = True
                
                return self.cached_scan_frames[:-1] + \
                       [self.calframe_for_info_on_lines]
            
            
            if "InformationOnLinesContained" in frame.keys():
                if self.scan_frame_counter == 0:
                    # Any calframe on lines that arrived
                    # before any scan frames arrive
                    # should not trigger line-notching process
                    # because there are no timestreams in the cache yet.
                    return
                
                core.log_info(
                    "\n", "In the RemoveLinesFromTimestreams module,",
                    "\n", "a frame having information on lines is received.", 
                    "\n", "The line-removal process is about to start.",
                    "\n",  unit="notchfiltertimestreams")
                
                line_frame = frame
                
                timestream_maps = \
                    [each_frame[self.original_timestream_map_key] \
                     for each_frame in self.cached_scan_frames]
                
                elevation_timestreams = \
                    [each_frame["RawBoresightEl"] \
                     for each_frame in self.cached_scan_frames]
                
                concatenated_timestream_maps, mean_of_each_timestream = \
                    get_concatenated_and_filtered_timestreams(
                        timestream_maps,
                        elevation_timestreams,
                        self.original_timestream_map_key,
                        "RawBoresightEl (doesn't matter)",
                        poly_order_time_domain_per_frame=0,
                        poly_order_elevation_domain=0,
                        return_mean_map_too=True)
                    
                    # The elevation timestreams were used in the
                    # FindLinesFromPSDs module with the purpose of
                    # reducing effects of el steps on timestreams
                    # so that the periodic el steps don't create
                    # lines in the PSDs that may be flagged.
                    # However, that doesn't seem necessary here.
                    # If the mean is not removed from a timestream
                    # prior to the filtering, though, the filtered version
                    # can still look bad in some cases
                    # (notching a line from a timestream
                    #  that has a non-zero mean sometimes results in
                    #  a timestream that has a new line at
                    #  a slightly different frequency),
                    # so poly order of zero is used here.
                
                notch_filtered_timestream_map = \
                    remove_lines_from_timestreams(
                        line_frame,
                        concatenated_timestream_maps,
                        self.n_data_pts_each_timestream,
                        line_frame["by_group"],
                        self.boloproma, self.wirma,
                        line_frame["by_band"],
                        line_frame["by_wafer"],
                        line_frame["by_board"],
                        line_frame["by_squid"])
                
                for bolo, ts in notch_filtered_timestream_map.items():
                    notch_filtered_timestream_map[bolo] = \
                        ts + mean_of_each_timestream[bolo]
                                
                deconcatenated_timestream_maps = \
                    split_timestreams(
                        self.n_data_pts_each_timestream,
                        notch_filtered_timestream_map)
                
                core.log_info(
                    "\n", "In the RemoveLinesFromTimestreams module,",
                    "\n", "filtered timestream maps are added to",
                    "\n", "each scan frame.",
                    "\n", unit="notchfiltertimestreams")
                
                for i in range(len(self.cached_scan_frames)):
                    if self.filtered_timestream_map_key      \
                    not in self.cached_scan_frames[i]:
                        self.cached_scan_frames[i]           \
                        [self.filtered_timestream_map_key] = \
                            deconcatenated_timestream_maps[i]
                    
                    start_time_original_timestream_map = \
                        self.cached_scan_frames[i]       \
                        [self.original_timestream_map_key].start
                    end_time_original_timestream_map = \
                        self.cached_scan_frames[i]     \
                        [self.original_timestream_map_key].stop
                    start_time_new_timestream_map = \
                        deconcatenated_timestream_maps[i].start
                    end_time_new_timestream_map =   \
                        deconcatenated_timestream_maps[i].stop
                    
                    if start_time_original_timestream_map != \
                       start_time_new_timestream_map:
                        core.log_info(
                            "\n", "- The two timestream maps in",
                            "scan frame", self.cached_scan_frame_numbers[i],
                            "\n", "  now have different start times!",
                            "\n", "   (original:",
                            str(start_time_original_timestream_map)+",",
                            "\n", "    new     :",
                            str(start_time_new_timestream_map)+")",
                            "\n", "  The new value will be changed",
                            "to the old one.",
                            "\n", unit="notchfiltertimestreams")
                        deconcatenated_timestream_maps[i].start = \
                            self.cached_scan_frames[i]            \
                            [self.original_timestream_map_key].start
                    
                    if end_time_original_timestream_map != \
                       end_time_new_timestream_map:
                        core.log_info(
                            "\n", "- The two timestream maps in",
                            "scan frame", self.cached_scan_frame_numbers[i],
                            "\n", "  now have different end times!",
                            "\n", "   (original:",
                            str(start_time_original_timestream_map)+",",
                            "\n", "    new     :",
                            str(start_time_new_timestream_map)+")",
                            "\n", "  The new value will be changed",
                            "to the old one.",
                            "\n", unit="notchfiltertimestreams")
                        deconcatenated_timestream_maps[i].stop = \
                            self.cached_scan_frames[i]           \
                            [self.original_timestream_map_key].stop
                    
                    # There are some cases where the start time of
                    # a pre-filtered timestream differs from that of
                    # its post-filtered version by severa nanoseconds.
                
                self.need_to_clear_cache = True
                
                return  self.cached_scan_frames[:-1] + [line_frame]
        
        
        if frame.type == core.G3FrameType.EndProcessing:
            core.log_info("\n", "In the RemoveLinesFromTimestreams module,",
                          "\n", "actually, the EndProcessing fr. has arrived,", 
                          "\n", "so the job of this notch filter is done!",
                          "\n", " (The one remaining frame in the cache",
                          "\n", "  will be returned.)",
                          "\n", unit="notchfiltertimestreams")
            
            return [self.cached_scan_frames[-1], frame]





def get_concatenated_and_filtered_timestreams(
        timestream_maps,
        elevation_timestreams,
        timestream_map_key,
        elevation_data_key,
        just_return_concatenated=False,
        poly_order_time_domain_per_frame=3,
        poly_order_elevation_domain=0,
        return_mean_map_too=False):
    
    '''
    This function is mainly intended to be used as part of
    the pipeline segment NotchFilterFFTVersion.
    The function does the following three things:
    
    It concatenates bolometer timestream maps stored in several scan frames
    that are contiguous in time. In addition, the timestreams of
    the telescope elevation stored in those frames are also concatenated.
    
    After that, it fits/subtracts two polymonials to/from
    the concatenated bolometer timestreams.
    One polynomial is in time domain,
    and the other is in (cosecant of) elevation domain
    (e.g. pW(t) for the first one, and pW(csc(el)) for the second).
    
    The purpose of also using the elevation domain polynomial for filtering is
    to attempt to reduce effects on a bolometer timestream due to changes
    in the telescope elevation (and thus changes in atmosphere loading).
    The effects may be significant if bolometer timestreams
    from multitle scan frames
    that include multiple elevation steps are concatenated.
    
    
    Parameters:
    -----------
    
    timestream_maps: list
        A list of G3TimestreamMaps that contain bolometer timestreams
        to be concatenated and filtered.
    
    elevation_timestreams: list
        A list of G3Timestreams that store the telescope elevation
        during the same time period.
    
    timestream_map_key: str
        The key name of the bolometer timestream maps to be concatenated.
        This argument is needed jusf for the purpose of
        producing some log messages.
    
    elevation_data_key: str
        The key name of thes telescope elevation data.
        Similar to above, this argument is needed
        just for the purpose of producing some log messages.
    
    just_return_concatenated: bool, optional
        If True, the function returns the simply concatenated
        bolometer timestream maps and does not do any filtering.
        If False, the function returns the concatenated and filtered
        timestream maps. Default is False.
    
    poly_order_time_domain_per_frame: int, optional
        The order of the time domain polynomial that will be used
        for the fitting/subtraction to/from
        the concatenated bolometer timestreams.
        This number is a per-frame number, so, if 7 is used here, 
        and if bolometer timestream maps from 5 scan frames are concatenated,
        then a time domain polynomial of order 35 will be 
        fitted/subtracted to/from each long timestream. Default is 3.
    
    poly_order_elevation_domain: int, optional
        The order of the elevation domain polynomial that will be used
        for fitting/subtraction to/from the concatenated bolometer timestreams.
        Unlike the previous argument, this number is not a per frame number.
        Default is 0.
    
    return_mean_map_too: bool, optional
        If true, a map double that records the mean of each timestream
        before it is filtered will be returned. This map will be useful
        when the means are added back to notch-filtered timestreams.
        Default is False.
    
    
    Returns:
    --------
    
    concatenated_timestream_maps: G3TimestreamMap
        Simply concatenated bolometer timestream maps
        (if just_return_concatenated is True).
    
    filtered_timestream_map: G3TimestreamMap
        Concatenated and also polynomial-filtered bolometer timestream maps
        (if just_return_concatenated is False).
    
    '''
    
    core.log_info(
        "\n", "- Concatenating", timestream_map_key, "from",
        len(timestream_maps), "frames...",
        "\n", unit="notchfiltertimestreams")
    
    core.log_info(
        "\n", "   (Number of bolometers in each timestream map:",
        "\n", "   ", [len(each_tsm) for each_tsm in timestream_maps],")",
        "\n", unit="notchfiltertimestreams")
    
    try:
        concatenated_timestream_maps = \
            core.G3TimestreamMap.concatenate(timestream_maps)
    
    except RuntimeError:
        core.log_info(
            "\n", "   (Trying to concatenate timestreams caused a",
            "RuntimeError.",
            "\n", "    It's possible that the timestreams are not", 
            "contiguous in time.)",
            "\n", unit="notchfiltertimestreams")
        
        # There were some observations where timestreams from
        # some consecutive scan frames were not actually cotinguous in time
        # (the timestamp of the end of one timestream is separated from
        # the timestamp of the start of the next timestream by more than
        # one sample). The function that concatenates timestreams throws
        # an error if this discontinuity is encountered. In that case,
        # line notching will not be done.
        
        return None
    
    
    mean_of_each_timestream = core.G3MapDouble()
    for bolo, ts in concatenated_timestream_maps.items():
        mean_of_each_timestream[bolo] = np.nanmean(ts)
    
    if just_return_concatenated:
        if return_mean_map_too:
            return concatenated_timestream_maps, mean_of_each_timestream
        else:
            return concatenated_timestream_maps
    
    
    concatenated_elevation_data     = \
        core.G3Timestream.concatenate(elevation_timestreams)
    concatenated_elevation_data[:]  = \
        1.0 / np.sin(concatenated_elevation_data/core.G3Units.rad)
    
    core.log_info(
        "\n", "- Filtering the concatenated data based on",
        elevation_data_key+" and so on...",
        "\n", "   (Orders of polynomials fitted and subtracted:",
        "\n", "     time domain = ",
        str(len(timestream_maps)*poly_order_time_domain_per_frame)+", ",
        "elevation domain = "+str(poly_order_elevation_domain)+")",
        "\n", unit="notchfiltertimestreams")
    
    filtered_timestream_map = \
        poly_filter_ts_data_with_abscissa(
            concatenated_elevation_data,
            concatenated_timestream_maps,
            make_empty_filter_mask(concatenated_timestream_maps),
            poly_order_elevation_domain,
            len(timestream_maps)*poly_order_time_domain_per_frame)
    
    if return_mean_map_too:
        return filtered_timestream_map, mean_of_each_timestream
    else:
        return filtered_timestream_map





def find_lines_from_PSDs(
        timestream_map,
        n_data_pts_each_timestream,
        bolometer_properties_map,
        wiring_map,
        from_by_group_median_PSD_only=True,
        by_band=False,
        by_wafer=True,
        by_board=False,
        by_squid=False,
        factor_for_thresholds=2.00,
        factor_for_baseline=1.20,
        power_law_fit_starting_frequency=0.06,
        power_law_fit_ending_frequency=0.60,
        white_noise_estimate_starting_frequency=4.5,
        white_noise_estimate_ending_frequency=7.5,
        low_frequency_cut=0.99,
        maximum_number_of_lines=50,
        save_PSDs=False,
        make_plots_of_lines_found=False,
        n_figures_to_generate=50,
        directory_to_save_figures=".",
        observation_frame=None,
        scan_frame_numbers=None):
    
    '''
    This function attemps to identify spectral lines from timestreams.
    The line-finding method used here is as follows:
    
    PSD of a timestream is calculated, and a model is constructed from the PSD.
    The model is a combination of two terms, one of which is a constant,
    and the other of which is a power law term.
    
    The constant is equal to the median PSD values in the frequency range
    specified by the arguments white_noise_estimate_starting_frequency and
    white_noise_estimate_ending_frequency.
    
    The amplitude and exponent of the power law term are obtained
    from a linear least squares fit to the logarithm of frequency values
    and PSD values in the frequency range specified by the arguments
    power_law_fit_starting_frequency and power_law_fit_ending_frequency.
    
    After the model is constructed, all the actual PSD values
    that are larger than the model values by a factor of the number
    specified by the argument factor_for_thresholds are flagged as outliers.
    
    For each outlier that corresponds to a local maximum,
    the PSD values at the neighboring frequency bins are used
    to estimate how wide this peak is.
    The first frequency bin on the left and the one on the right
    at which the PSD values drop below the number specified
    by the argument factor_for_baseline times
    the model values define the width of the peak.
    
    Then, spectral lines found from all the bolometers are
    recorded in a calibration frame that has the following keys:
    
    "FrequencyBins" points to a G3VectorDouble that looks like
    [frequency bin 1, frequency bin 2, ...]
    
    "LineLocations" points to a G3MapVectorString that looks like
    {bolometer 1: [line location 1, line location 2, ...],
     bolometer 2: [line location 1, line location 2, ...],
     bolometer 3: ...}
    
    "LineLeftEdges" points to a G3MapMapDouble that looks like
    {bolometer 1: {line location 1: loc. of left edge of line 1,
                   line location 2: loc. of left edge of line 2,
                   line location 3: ...},
     bolometer 2: ...}
    
    "LineRightEdges" points to a G3MapMapDouble that looks like
    {bolometer 1: {line location 1: loc. of right edge of line 1,
                   line location 2: loc. of right edge of line 2,
                   line location 3: ...},
     bolometer 2: ...}
    
    "LineRelativeHeights" points to a G3MapMapDouble that looks like
    {bolometer 1: {line location 1: height of line 1,
                   line location 2: height of line 2,
                   line location 3: ...},
     bolometer 2: ...}
    
    
    Parameters:
    -----------
    
    timestream_map: G3TimestreamMap
        The timestream map from which lines are to be found.
    
    n_data_pts_each_timestream: list
        A list of integers that represent the length of each of the
        timestreams that were concatenated.
    
    from_by_group_median_PSD_only: bool, optional
        If True, lines will be found from only the median PSDs
        of each group of timestreams instead of
        from the PSD of every detector timestream. Default is True.
    
    by_band, by_wafer, by_board, by_squid: bools, optional
        These arguments specify how timestreams should be grouped
        if the argument above is True. Defaults are such that
        only by_wafer is True.
    
    bolometer_properties_map: BolometerPropertiesMap
    wiring_map: WiringMap
        In case lines are to be found from median PSDs of
        groups of timestreams, these are needed
        to separate detector into groups.
    
    factor_for_thresholds: float, optional
        As mentioned above, this is a number specifying
        how much deviation a PSD value needs to have from the model value
        in order to be regarded as part of a line. Default is 2.00.
    
    factor_for_baseline: float, optional
        As mentioned above, this is a number used to
        estimate how wide each line is. Default is 1.20.
    
    power_law_fit_starting_frequency,
    power_law_fit_ending_frequency  : floats, optional
        These arguments define the frequency range used for
        the power-law fit described above. The numbers are in the units of Hz.
        Defaults are 0.06 and 0.60.
    
    white_noise_estimate_starting_frequency,
    white_noise_estimate_ending_frequency  : floats, optional
        These arguments define the frequency range used to get
        the constant term in the model described above.
        The numbers are in the units of Hz. Defaults are 4.5 and 7.5.
    
    low_frequency_cut: float, optional
        Any outliers found at frequencies below this value will be ignored.
        This value should be expressed in the units of Hertz. Default is 0.99.
    
    maximum_number_of_lines: int, optional
        The code stops finding lines once
        this many lines have been found from a PSD. Default is 50.
    
    save_PSDs: bool, optional
        If True, in addition to the calibration frame that contains information
        on lines, another frame that records the median PSD from each timestream
        will be created and returned. Default is False.
    
    make_plots_of_lines_found: bool, optional
        If True, figures will be made that show the lines found
        from the timestreams of select bolometers, which can be useful to see
        the performance of the line-finder. Default is False.
        All the arguments below are related to making figures.
    
    n_figures_to_generate: int, optional
        Approximate number of figures to generate. Figures will be made
        for only one in N bolometers, where N is roughly equal to 
        total number of bolometers / n_figures_to generate. Default is 50.
    
    directory_to_save_figures: str, optional
        The location where the figures will be saved.
        The file name of each figure will be something like
        "ra0hdec-44.75_scan_001_lines_found_from_2019.xyz.png"
        Default is '.'.
    
    observation_frame: G3Frame, optional
    scan_frame_numbers: list, optional
        The title of one of the plots in a given figure for a given bolometer
        needs information on which scans of which observation
        the timestream is from. The scan_frame_numbers needs to be a list
        of monotonically increasing integers. Defaults are None.
    
    
    Returns:
    --------
    
    calframe_for_lines: G3Frame
        A calibration frame that contains the information on lines found.
    
    frame_for_psds: G3Frame
        If save_PSDs is True, then another frame that contains PSDs
        will be returned as well.
    
    '''
    
    core.log_info(
        "\n", "- Finding lines from power spectral densities...",
        "\n", unit="notchfiltertimestreams")
    
    
    # - Take FFT first
    
    if len(n_data_pts_each_timestream) == 1:
        wf_callable = np.hanning
    else:
        custom_window_function = \
            make_tukey_window_w_turnarounds_tapered(n_data_pts_each_timestream)
        def wf_callable(ts_len):
            # The get_psd_of_ts_map function
            # does not accept an array as the window function,
            # so here is a dummy function that returns the
            # custom window function made above,
            # which is the window function used for the line-removal process.
            return custom_window_function
    
    core.log_info(
        "\n", "   (Taking FFT...)", "\n", unit="notchfiltertimestreams")
    
    if from_by_group_median_PSD_only:
        psd_map, frequency_values = \
            get_psd_of_ts_map(timestream_map,
                              window_function=wf_callable,
                              pad=False)        
        core.log_info(
            "\n", "   (Splitting PSD map by group...)",
            "\n", unit="notchfiltertimestreams")
        psd_map = split_map_by_group(
                      psd_map, bolometer_properties_map, wiring_map,
                      by_band=by_band, by_wafer=by_wafer,
                      by_squid=by_squid, by_board=by_board)
        core.log_info(
            "\n", "   (Calculating median PSD for each group...)",
            "\n", unit="notchfiltertimestreams")
        psd_map, n_detectors_map = get_median_for_each_split_map(psd_map)
    
    else:
        psd_map, frequency_values = \
            get_psd_of_ts_map(timestream_map,
                              window_function=wf_callable,
                              pad=False)
        n_dectors_map = {}
    
    if save_PSDs:
        frame_for_psds = core.G3Frame()
        for group_id, psd in psd_map.items():
            frame_for_median_psds["PSDofTimestreamsfrom"+group_id] = \
                core.G3VectorDouble(psd)
        for group_id, n_bolo in n_detectors_map.items():
            frame_for_median_psds["NumberofBolometersfrom"+group_id] = \
                core.G3Int(n_detectors)
        frame_for_psds["FrequencyValuesForPSDs"] = \
            core.G3VectorDouble(frequency_values)
    
    frequency_values = frequency_values / core.G3Units.Hz
    delta_f          = frequency_values[1] - frequency_values[0]
    
    
    # - Then find lines
    
    ## - Prepare a few things
    
    calframe_for_lines = core.G3Frame(core.G3FrameType.Calibration)
    calframe_for_lines["LineLocations"]  = core.G3MapVectorString()
    calframe_for_lines["LineLeftEdges"]  = core.G3MapMapDouble()
    calframe_for_lines["LineRightEdges"] = core.G3MapMapDouble()
    calframe_for_lines["LineRelativeHeights"] = core.G3MapMapDouble()
    calframe_for_lines["by_group"] = from_by_group_median_PSD_only
    calframe_for_lines["by_band"]  = by_band
    calframe_for_lines["by_wafer"] = by_wafer
    calframe_for_lines["by_board"] = by_board
    calframe_for_lines["by_squid"] = by_squid
    
    core.log_info(
        "\n", "   (Constructing a model for each PSD and",
        "\n", "    Finding lines and estimating their width",
        "based on these values: ",
        "\n", "    threshold = "+str(factor_for_thresholds)+" * model,",
        "baseline = "+str(factor_for_baseline)+" * model...)",
        "\n", unit="notchfiltertimestreams")
    
    f_s_pl = power_law_fit_starting_frequency
    f_e_pl = power_law_fit_ending_frequency
    i_s_pl = int(np.argmin(np.absolute(frequency_values - f_s_pl)))
    i_e_pl = int(np.argmin(np.absolute(frequency_values - f_e_pl)))
    
    f_s_wn = white_noise_estimate_starting_frequency
    f_e_wn = white_noise_estimate_ending_frequency
    i_s_wn = int(np.argmin(np.absolute(frequency_values - f_s_wn)))
    i_e_wn = int(np.argmin(np.absolute(frequency_values - f_e_wn)))
    
    
    n_bolometers_encountered = 0
    n_figures_generated      = 0
    
    for bolo, real_psd in psd_map.items():
        # if from_by_group_median_PSD_only,
        # then the variable each_bolometer actually refers to a group
        
        n_bolometers_encountered += 1
        
        calframe_for_lines["LineLocations"][bolo]  = core.G3VectorString()
        calframe_for_lines["LineLeftEdges"][bolo]  = core.G3MapDouble()
        calframe_for_lines["LineRightEdges"][bolo] = core.G3MapDouble()
        calframe_for_lines["LineRelativeHeights"][bolo] = core.G3MapDouble()
        
        ## - Construct a model for each PSD
        
        model_psd = np.zeros(len(frequency_values))
        
        wnl = white_noise_level  = np.median(real_psd[i_s_wn:i_e_wn+1])
        
        frequency_values_for_fit = frequency_values[i_s_pl:i_e_pl+1]
        psd_values_for_fit       = real_psd[i_s_pl:i_e_pl+1]
        power_law_fit_results    = stats.linregress(
                                       np.log10(frequency_values_for_fit),
                                       np.log10(psd_values_for_fit))
        exp = power_law_fit_results[0]
        amp = np.power(10, power_law_fit_results[1])
        
        power_law_psd = amp * np.power(frequency_values[1:], exp)
        power_law_psd = np.insert(power_law_psd, 0, power_law_psd[0])
        # To avoid the power-law part become inf at f = 0.
        
        
        if from_by_group_median_PSD_only:
        
            kernel_size  = int(power_law_fit_ending_frequency*1.8 / delta_f)
            if kernel_size % 2 == 0:
                kernel_size += 1
            medfiltered_psd = signal.medfilt(real_psd, kernel_size=kernel_size)
            
            # The power-law-plus-constant-term model does not seem to 
            # fit a median PSD from a wafer very well because
            # the PSD has a decreasing trend as the frequency increases,
            # so, the smoothed PSD obtained from a median filter will be used
            # to replace the constant-term part.
            #
            # To stitch the power-law part and the result from the
            # median filter together, the kernel size should not be
            # more than twice as large as the frequency range from 0 Hz to
            # the number specified by power_law_fit_ending_frequency.
            # Otherwise, the flat part that exists in the beginning of
            # the result from the median filter
            # makes it hard for the power law part to join smoothly.
            
            i_switch = int(i_e_pl * 1.25)
            model_psd[1:i_switch] = power_law_psd[1:i_switch] - \
                                    power_law_psd[i_switch]   + \
                                    medfiltered_psd[i_switch]
            model_psd[i_switch:]  = medfiltered_psd[i_switch:]
            
            # The power-law PSD and the smoothed PSD are stitched together
            # at a frequency near power_law_fit_ending_frequency.
        
        else:
            # The power-law-plus-constant-term model seems to work better
            # for the PSD from one detector than the median PSD,
            # so the model will still be used.
            
            one_over_f_knee_freq = np.power(0.5*wnl/amp, 1/exp)
            
            ik = knee_freq_index = \
                np.argmin(np.absolute(frequency_values - one_over_f_knee_freq))
            model_psd[1:ik] = power_law_psd[1:ik] + 0.5 * wnl
            model_psd[ik:]  = wnl
            
            # The power-law part and the constant part
            # are stitched together at the frequency where
            # the power-law term is equal to half the constant term.
            # The 0.5 used in the definition of "1/f knee" is arbitrary.
            # Having a method that finds the optimal value would be nice.
        
        
        ## - Flag large PSD values based on the model
        
        outliers   = []
        thresholds = model_psd * factor_for_thresholds
        
        for i in range(len(frequency_values)-1):
            if (real_psd[i] > thresholds[i]) and \
               (frequency_values[i] > float(low_frequency_cut)):               
                if (real_psd[i] > real_psd[i+1]) and \
                   (real_psd[i] > real_psd[i-1]):
                    outliers.append((
                        frequency_values[i],
                        real_psd[i]/model_psd[i]))
        
        outliers = sorted(outliers, key=itemgetter(1), reverse=True)
        
        frequency_ranges_of_lines_found_so_far = []
        baseline_values = model_psd * factor_for_baseline
        
        for counter, each_outlier in enumerate(outliers, 1):
            if counter == int(maximum_number_of_lines):
                break
            
            location    = each_outlier[0]
            peak_height = each_outlier[1]
            i_location  = int(np.where(frequency_values == location)[0][0])
            
            can_skip_this_one = False
            for each_range in frequency_ranges_of_lines_found_so_far:
                left_edge  = each_range[0]
                right_edge = each_range[1]
                if (location > left_edge) and (location < right_edge):
                    can_skip_this_one = True
                    break
            
            if can_skip_this_one:
                continue
            
            # If an outlier lies within the width of an existing line,
            # this outlier is not regarded as a new line.
            
            # For each line,
            # PSD values in the neighboring frequency bins
            # are compared with the corresponding baseline values
            # for the purpose of estimating how wide the line is.
            
            j = 0
            while real_psd[i_location+j] > baseline_values[i_location+j]:
                if (i_location + j) == (len(real_psd) - 1):
                    break
                j += 1
            right_edge_of_this_line = frequency_values[i_location+j]
            
            j = 0
            while real_psd[i_location-j] > baseline_values[i_location-j]:
                if (i_location - j) == 0:
                    break
                j += 1
            left_edge_of_this_line  = frequency_values[i_location-j]
            
            frequency_ranges_of_lines_found_so_far.\
                append((left_edge_of_this_line, right_edge_of_this_line))
            
            hz = core.G3Units.Hz
            calframe_for_lines["LineLocations"][bolo].append(str(location))
            calframe_for_lines["LineLeftEdges"][bolo][str(location)] = \
                left_edge_of_this_line * hz
            calframe_for_lines["LineRightEdges"][bolo][str(location)] = \
                right_edge_of_this_line * hz
            calframe_for_lines["LineRelativeHeights"][bolo][str(location)] = \
                peak_height
        
        
        ## - Make plots of lines found if requested
        
        if make_plots_of_lines_found:
        
            make_plots_for_this_one = False
            
            one_in_how_many_bolos_to_check = int(np.round(
                len(psd_map.keys()) / n_figures_to_generate))
            
            if (len(psd_map.keys()) < n_figures_to_generate) or \
               (n_bolometers_encountered % one_in_how_many_bolos_to_check == 0):
                if bolo in timestream_map.keys():
                    if False not in np.isfinite(timestream_map[bolo]):
                        make_plots_for_this_one = True
                else:
                    make_plots_for_this_one = True
            
            if make_plots_for_this_one:
                n_figures_generated += 1
                if from_by_group_median_PSD_only:
                    core.log_info(
                        "\n", "   (Making plots of lines for "+bolo,
                        "[figure "+str(n_figures_generated)+"] ...)",
                        "\n", unit="notchfiltertimestreams")
                else:
                    core.log_info(
                        "\n", "   (Making plots of lines for "+\
                        bolometer_properties_map[bolo].physical_name,
                        "[figure "+str(n_figures_generated)+"] ...)",
                        "\n", unit="notchfiltertimestreams")
                
                if from_by_group_median_PSD_only:
                    timestream = timestream_map.values()[0]
                    no_timestream_plot = True
                else:
                    timestream = timestream_map[bolo]
                    no_timestream_plot = False
                
                plot_lines_found(
                    bolo,
                    timestream,
                    frequency_values, delta_f,
                    psd_map[bolo], model_psd,
                    factor_for_baseline, factor_for_thresholds,
                    i_s_pl, i_e_pl, i_s_wn, i_e_wn,
                    low_frequency_cut,
                    calframe_for_lines, bolometer_properties_map,
                    observation_frame, scan_frame_numbers,
                    directory_to_save_figures,
                    no_timestream_plot=no_timestream_plot)
    
    if save_PSDs:
        return frame_for_psds, calframe_for_lines
    else:
        return calframe_for_lines





def split_map_by_group(
        map_to_split,
        bolo_props_map,
        wir_map,
        by_band=False,
        by_wafer=True,
        by_board=False,
        by_squid=False):
    
    '''
    This function splits a map (G3MapVectorDouble, G3TimestreamMap, etc.)
    by template groups, which can be some combination of
    bands, wafers, IceBoards, and SQUIDs.
    
    Parameters:
    -----------
    
    map_to_split: certain G3Map object 
        In the context of the notch filter, this will be a PSD map,
        i.e. G3MapVectorDouble.
    
    bolo_props_map: BolometerPropertiesMap
    wir_map: WiringMap
        These are needed to group bolometers.
    
    by_band, by_wafer, by_board, by_squid: bools, optional
        Booleans specifying how detectors should be grouped.
        Defaults are such that only by_wafer is True.
    
    
    Returns:
    --------
    
    dictionary_of_maps: dict
        A python dictionary whose keys correspond to groups
        and whose values correspond maps
        associated with the appropriate groups.
    
    '''
    
    groups = get_template_groups(
                 bolo_props_map, wir_map,
                 per_band=by_band, per_wafer=by_wafer,
                 per_board=by_board, per_squid=by_squid,
                 include_keys=True, exclude_empty=True)
    
    dictionary_of_maps = \
        {each_group: type(map_to_split)() for each_group in groups.keys()} 
    
    tg_inv = {}
    for group, bolos in groups.items():
        for bolo in bolos:
            tg_inv[bolo] = group
    
    available_bolos = set(map_to_split.keys()) & set(tg_inv.keys())
    for bolo in available_bolos:
        dictionary_of_maps[tg_inv[bolo]][bolo] = map_to_split[bolo]
    
    bad_groups = []
    for group, map_obj in dictionary_of_maps.items():
        if len(map_obj) == 0:
            bad_groups.append(each_group)
    for group in bad_groups:
        dictionary_of_maps.pop(each_group, None)
    dictionary_of_maps.pop("All", None)
    dictionary_of_maps.pop("0_", None)
    
    return dictionary_of_maps





def get_median_for_each_split_map(split_maps):
    
    '''
    This function takes a dictionary of maps such as PSD maps
    that are grouped by certain attribute of bolometers (wafer, SQUID, etc.)
    and calculates the median for each group.
    
    
    Parameters:
    -----------
    
    split_maps: dict
        A dictionary whose keys are groups (e.g. wafers) and
        whose values are maps of timestreams, PSDs, etc.
        associated with the corresponding groups.
    
    
    Returns:
    --------
    
    by_group_median_map: G3Map object
        A map that contains the median PSD, etc. from each group.
    
    '''
    
    by_group_median_map = \
        {group: np.nanmedian(map_obj.values(), axis=0)
         for group, map_obj in split_maps.items()}
    
    n_detectors_map = \
        {group: len(map_obj.keys())
         for group, map_obj in split_maps.items()}
    
    return by_group_median_map, n_detectors_map





def plot_lines_found(
        boloname,
        timestream,
        frequency_values, delta_f,
        actual_psd, model_psd,
        factor_for_baseline, factor_for_thresholds,
        i_s_pl, i_e_pl, i_s_wn, i_e_wn,
        low_frequency_cut,
        information_on_lines, bolometer_properties,
        obs_fr, scan_frames_info,
        directory_to_save_figures, no_timestream_plot=False):
    
    '''
    This function makes plots of lines found from a PSD,
    and it is intended to be called by the find_lines_from_PSDs function above.
    
    It makes a plot showing the corresponding timestream (if applicable)
    and then another plot showing the PSD with the found lines indicated.
    Finally, it creates one smaller plot showing each of the lines found.
    
    
    Parameters:
    -----------
    
    Since all the arguments needed are defined
    in the find_lines_from_PSDs function, and since this plotting code is
    probably not useful unless called by that particular function,
    explaining each argument does not seem very meaningful so will be omitted.
    
    '''
    
    import matplotlib.pyplot as pyplot
    import matplotlib
    matplotlib.rcParams["mathtext.default"] = "regular"
    
    # - Prepare the title of the figure and the name of the png file
    
    source     = obs_fr["SourceName"]
    obsid      = str(obs_fr["ObservationID"])
    start_date = str(obs_fr["ObservationStart"]).split(":")[0].split("-")[0:2]
    start_date = start_date[1] + " " + start_date[0]
    fr_first   = "{:03d}".format(scan_frames_info[0])
    fr_last    = "{:03d}".format(scan_frames_info[-1])
    
    
    # - Organize data to be plotted
    
    sampling_rate = timestream.sample_rate / core.G3Units.Hz
    delta_t       = 1.0 / sampling_rate
    time_values   = np.arange(len(timestream)) * delta_t
    start_time    = timestream.start
        
    if timestream.units == core.G3TimestreamUnits.Power:
        ts_units  = core.G3Units.pW / 1000.0
        tsu_str   = "Timestream [fW]"
        psd_units = core.G3Units.aW * core.G3Units.aW / core.G3Units.Hz
        psdu_str  = r"$PSD [{aW}^2 / Hz]$"
    elif timestream.units == core.G3TimestreamUnits.Tcmb:
        ts_units  = core.G3Units.mK
        tsu_str   = "Timestream [mK]"
        psd_units = core.G3Units.mK * core.G3Units.mK / core.G3Units.Hz
        psdu_str  = r"$PSD [${mK}^2 / Hz]$"
    elif timestream.units == core.G3TimestreamUnits.Counts:
        ts_units  = 1.0
        tsu_str   = "Timestream [Counts]"
        psd_units = 1.0 / core.G3Units.Hz
        psdu_str  = r"$PSD [{Counts}^2 / Hz]$"
    else:
        ts_units  = 1.0
        tsu_str   = "Timestream [Arb.]"
        psd_units = 1.0
        psdu_str  = "PSD [Arb.]"
    
    timestream /= ts_units
    actual_psd /= psd_units
    model_psd  /= psd_units
    
    
    # - Decide the size of the figures
    
    number_of_lines = len(information_on_lines["LineLocations"][boloname])
    n_lines_orig    = number_of_lines
    if number_of_lines > 24:
        number_of_lines = 24
    n_additional_rows = int(np.ceil(number_of_lines / 4.0))
    n_total_rows      = 2 + n_additional_rows
    figure_height     = 3 * n_total_rows
    figure_for_this_bolo = pyplot.figure(figsize=(14, figure_height))
    
    
    # - Plot the timestream
    
    plot_of_timestream = figure_for_this_bolo.add_subplot(n_total_rows, 2, 1)
    
    if not no_timestream_plot:
        plot_of_timestream.plot(time_values, timestream,
                                color="black", linewidth=0.1, alpha=0.9)
        plot_of_timestream.set_ylabel(tsu_str, fontsize=11)
        plot_of_timestream.set_xlabel("Time [s]", fontsize=11)
        plot_of_timestream.set_title(
            "Poly filtered timestream from "+\
             bolometer_properties[boloname].physical_name.upper()+\
             "/"+boloname+"\n"+\
             "("+source+", "+obsid+", "+str(start_time).split(".")[0]+", "+\
             "scan frames "+fr_first+" to "+fr_last+")", fontsize=11)
    else:
        plot_of_timestream.set_xticklabels([])
        plot_of_timestream.set_yticklabels([])
        plot_of_timestream.set_title("Timestream N/A", fontsize=11)
    
    
    # - Make a histogram of the timestream values
    
    k2, p_normality = stats.mstats.normaltest(timestream)
    if p_normality == 0.0:
        log10_p = "-inf"
    else:
        log10_p = str(np.round(np.log10(p_normality), 2))
    pctl_98 = np.nanpercentile(timestream, 98)
    pctl_02 = np.nanpercentile(timestream,  2)
    xlimpos = np.max([np.abs(pctl_98), np.abs(pctl_02)])
    
    histogram_of_ts = figure_for_this_bolo.add_subplot(n_total_rows, 2, 2)
    if not no_timestream_plot:
        n, bins, patches = histogram_of_ts.hist(
            timestream, bins=100, range=[-xlimpos, xlimpos],
            histtype="step", color="black", alpha=0.9, linewidth=1.0)
        histogram_of_ts.set_xlim(left=-xlimpos, right=xlimpos)
        histogram_of_ts.set_xlabel(tsu_str, fontsize=11)
        histogram_of_ts.set_ylim(top=1.1*np.max(n))
        histogram_of_ts.set_yticklabels([])
        histogram_of_ts.set_title(
            "Histogram of timestream values"+"\n"+\
            "(Results of scipy's normaltest: "+\
            "log10(p-value) = "+log10_p+")", fontsize=11)
    else:
        histogram_of_ts.set_xticklabels([])
        histogram_of_ts.set_yticklabels([])
        histogram_of_ts.set_title("Histogram N/A", fontsize=11)
    
    
    # - Plot the lines found from the PSD
    
    plot_of_psd = figure_for_this_bolo.add_subplot(n_total_rows, 1, 2)
    plot_of_psd.plot(frequency_values, actual_psd,
                     color="black", linewidth=0.3, alpha=0.9)
    for counter, line \
    in  enumerate(information_on_lines["LineLocations"][boloname]):
        if counter == 0:
            plot_of_psd.axvline(
                float(line), label="Lines found",
                color="orange", linewidth=1.0, alpha=0.5)
        else:
            plot_of_psd.axvline(
                float(line), color="orange", linewidth=2.0, alpha=0.1)
    
    plot_of_psd.plot(frequency_values, model_psd, label="Model PSD",
                     color="green", linewidth=1.0, alpha=0.4)
    plot_of_psd.plot(frequency_values, factor_for_thresholds*model_psd,
                     label="Thresholds for peak finding "+\
                           "("+str(factor_for_thresholds)+" * model)",
                     color="red", linewidth=1.0, alpha=0.4)
    plot_of_psd.plot(frequency_values, factor_for_baseline*model_psd,
                     label="Baseline values for peak width estimation "+\
                           "("+str(factor_for_baseline)+" * model)", 
                     color="yellow", linewidth=1.0, alpha=1.0)
    plot_of_psd.plot(frequency_values[i_s_pl:i_e_pl+1],
                     actual_psd[i_s_pl:i_e_pl+1],
                     label="Data used for power law fit and "+\
                           "noise floor estimate",
                     color="blue", linewidth=0.3, alpha=1.0)
    plot_of_psd.plot(frequency_values[i_s_wn:i_e_wn+1],
                     actual_psd[i_s_wn:i_e_wn+1],
                     color="blue", linewidth=0.3, alpha=1.0)
    plot_of_psd.axvline(low_frequency_cut,
                        label="Frequency below which no peak-finding occurred",
                        color="grey", linewidth=0.3, alpha=1.0)
    plot_of_psd.set_yscale("log")
    plot_of_psd.set_xscale("log")
    plot_of_psd.set_ylim(top=1.20*np.max(actual_psd),
                         bottom=0.10*np.median(actual_psd))
    plot_of_psd.set_xlim(left=2e-2, right=80)
    plot_of_psd.tick_params(axis="both", labelsize=10)
    plot_of_psd.grid(axis="y", which="major", linewidth=0.1)
    plot_of_psd.grid(axis="x", which="both", linewidth=0.1)
    plot_of_psd.legend(loc="lower left", fontsize=10, framealpha=0.0)
    plot_of_psd.set_ylabel(psdu_str, fontsize=11)
    plot_of_psd.set_xlabel("Frequency [Hz]", fontsize=11)
    if no_timestream_plot:
        title = "\n"+"Plot of "+str(n_lines_orig)+" "+\
                "lines found from "+\
                "median of PSDs of all timestreams associated with "+\
                boloname.upper()+\
                " (brightest, up to 24, lines "+\
                "shown in smaller plots below)"
    else:
        title = "\n"+"Plot of "+str(n_lines_orig)+" "+\
                "lines found from PSD of the timestream "+\
                "(brightest (up to 24) lines"+\
                "shown in smaller plots below)"

    plot_of_psd.set_title(title, fontsize=11)
    
    
    # - Plot each line
    
    i = 0
    locs_and_heights = \
        [(float(location), float(height)) \
         for location, height in                          \
         information_on_lines["LineRelativeHeights"][boloname].items()]
    locs_and_heights = \
        sorted(locs_and_heights, key=itemgetter(1), reverse=True)
    
    for line in [pair[0] for pair in locs_and_heights]:
        i += 1
        if i == 25: break
        le = information_on_lines["LineLeftEdges"][boloname][str(line)]
        re = information_on_lines["LineRightEdges"][boloname][str(line)]
        lh = information_on_lines["LineRelativeHeights"][boloname][str(line)]
        le /= core.G3Units.Hz
        re /= core.G3Units.Hz
        li = int(np.argmin(np.absolute(frequency_values - le)))
        ri = int(np.argmin(np.absolute(frequency_values - re)))
        li = li - 1 * (ri-li)
        ri = ri + 1 * (ri-li)
        if li < 0: li = 0
        if ri > (len(frequency_values) - 1): ri = len(frequency_values) - 1
        
        frequency_values_to_plot = frequency_values[li:ri]
        psd_values_to_plot       = actual_psd[li:ri]
        
        plot_for_this_line = \
            figure_for_this_bolo.add_subplot(n_total_rows, 4, 8+i)
        plot_for_this_line.plot(
            (frequency_values_to_plot-line)/delta_f,
            psd_values_to_plot,
            color="black", alpha=0.7, linewidth=0.4,
            marker=".", markersize=5.0)
        plot_for_this_line.axvline(
            (line-line)/delta_f,
            color="red", alpha=0.5, linewidth=0.4)
        plot_for_this_line.axvline(
            (le-line)/delta_f,
            label="Estimated"+"\n"+"width",
            color="blue", alpha=0.5, linewidth=0.4)
        plot_for_this_line.axvline(
            (re-line)/delta_f,
            color="blue", alpha=0.5, linewidth=0.4)
        plot_for_this_line.tick_params(axis="both", labelsize=10)
        if lh > 50: plot_for_this_line.set_yscale("log")
        
        plot_for_this_line.legend(
            loc="upper right", fontsize=9, framealpha=0.0)
        plot_for_this_line.set_xlabel(
            "(f - f_line)/("+str(np.round(delta_f*1000, 4))+" mHz)",
            fontsize=11)
        plot_for_this_line.set_ylabel(
            psdu_str, fontsize=11)
        plot_for_this_line.set_title(
            "f_line = "+str(np.round(line, 4))+" Hz ",
            fontsize=11)
    
    
    # - Finish up
    
    figure_for_this_bolo.tight_layout(h_pad=1.0, w_pad=1.0)
    figure_for_this_bolo.savefig(
        directory_to_save_figures+"/"+\
        source+"-"+obsid+"-"+boloname+"-"
        "plots_of_spectral_lines_found_in_the_timestream_from_"+\
        "scan_frames_"+fr_first+"_to_"+fr_last+".png",
        bbox_inches="tight")
    pyplot.close(figure_for_this_bolo)





def remove_lines_from_timestreams(
        calframe_for_information_on_lines,
        original_timestream_map,
        n_data_pts_each_timestream,
        by_group_filters_only,
        bolometer_properties_map,
        wiring_map,
        by_band,
        by_wafer,
        by_board,
        by_squid):
    
    '''
    This function attempts to remove lines from timestreams.
    
    It first uses the information on lines recorded in a calibration frame
    generated by the function find_lines_from_PSDs above
    to construct a filter function in the frequency domain for each timestream.
    
    The value of the filter function is equal to 0.0
    if the corresponding frequency bin lies within the width of a line
    and is equal to 1.0 otherwise.
    
    Then, the filters are applied to timestreams,
    and the resulttant timestreams are returned.
    
    
    Parameters:
    -----------
    
    calframe_for_information_on_lines: G3Frame
        A calibration frame generated by the function find_lines_from_PSDs.
        As a result, this function assumes a specific type of data structure
        in which information on lines is recorded.
    
    original_timestream_map: G3TimestreamMap
        The timestream map that contains timestreams needing notch filtering.
    
    n_data_pts_each_timestream: list
        A list that contains information on how many consecutive scan frames
        were joined to obtain the timestreams in the original_timestream_map.
        For example, if data from three scan frames
        (Turnaround + Scan + Turnaround) were used,
        this list may look like [832, 7584, 832].
    
    by_group_filter_only: bool
        If True, the calibration frame that contains information on lines
        actually contains lines found from median PSDs of groups of timestreams,
        i.e., not from every detector timestream's PSD.
        As a result, instead of applying ~10,000 different filters to
        ~10,000 different timestreams, the same filter is applied to
        all the timestreams associated with one particular group.
    
    bolometer_properties_map: BolometerPropertiesMap
    wiring_map: WiringMap
        These are needed to separate timestrems into groups.
    
    by_band, by_wafer, by_board, by_squid: bools
        These booleans specify how timestream should be grouped.
    
    
    Returns:
    --------
    
    filtered_timestream_map: G3TimestreamMap
        The notch-filtered version of the original_timestream_map.
    
    '''
    
    core.log_info(
        "\n", "- Removing spectral lines from timestreams...",
        "\n", unit="notchfiltertimestreams")
    
    
    # - Construct a window function
    
    core.log_info(
        "\n", "   (Constructing a window function based on",
        "\n", "    how many scan frames' worth timestreams", 
        "\n", "    were concatenated to obtain the input timestream map.)",
        "\n", unit="notchfiltertimestreams")
    
    window_function = \
        make_tukey_window_w_turnarounds_tapered(n_data_pts_each_timestream)
    
    
    # - Construct filters
    
    core.log_info(
        "\n", "   (Constructing a filter function for each timestream...)",
        "\n", unit="notchfiltertimestreams")
    
    frequency_values = \
        fft_freqs(original_timestream_map.sample_rate,
                  original_timestream_map.n_samples) \
                  [:original_timestream_map.n_samples//2+1]
    frequency_values[-1] = np.absolute(frequency_values[-1])
    frequency_values     = np.asarray(frequency_values)
    # Keep the positive frequencies only
    
    filter_function_map = core.G3MapVectorComplexDouble()
    
    for key in calframe_for_information_on_lines.keys():
        if "LeftEdge" in key:
            line_left_edge_key = key
        elif "RightEdge" in key:
            line_right_edge_key = key
        elif "Location" in key:
            line_location_key = key
    
    for bolo, lines \
    in  calframe_for_information_on_lines[line_location_key].items():
        # if by_group_filter_only, then
        # each_bolometer actually refers to each group.
        
        filter_function = np.array([1.0] * len(frequency_values))
        
        for each_line in lines:
            f_left_edge  = \
                calframe_for_information_on_lines \
                [line_left_edge_key][bolo][each_line]
            f_right_edge = \
                calframe_for_information_on_lines \
                [line_right_edge_key][bolo][each_line]
            indices_to_notch = \
                np.where((frequency_values >= f_left_edge ) &
                         (frequency_values <= f_right_edge))[0]
            for each_index in indices_to_notch:
                filter_function[each_index] = 0.0
        
        f_strs = [str(np.round(f/core.G3Units.Hz, 4))
                  for f in sorted(list(
                      np.asarray(frequency_values)\
                      [np.where(filter_function==0.0)]))]
        n_f    = len(f_strs)
        f_strs = " ".join(f_strs)
        try:
            bid = bolo + "/" + bolometer_properties_map[bolo].physical_name
        except:
            bid = bolo
        if len(f_strs) != 0:
            core.log_info("\n", "   (Lines to be notched for", bid+":",
                          "\n", "    * frequencies [Hz] to notch:",
                          "\n", "   ", f_strs,
                          "\n", "    * in total", n_f, "values.)",
                          "\n", unit="notchfiltertimestreams")
        
        if int(original_timestream_map.n_samples) % 2 == 0:
            additional_filter_values = filter_function[1:-1]
        if int(original_timestream_map.n_samples) % 2 == 1:
            additional_filter_values = filter_function[1:]
        additional_filter_values = additional_filter_values[::-1]
        complete_filter_function = \
            np.concatenate((filter_function, additional_filter_values))
        # Expand the filter function to the negative frequency range as well.
        # Depending on whether the number of the data points of the timestream
        # is even or odd, the procedure differs slightly.
        
        complete_filter_function = \
            core.G3VectorComplexDouble(complete_filter_function)
        
        filter_function_map[bolo] = complete_filter_function
    
    if by_group_filters_only:
        by_group_timestream_maps = \
            split_map_by_group(original_timestream_map,
                               bolometer_properties_map, wiring_map,
                               by_band=by_band, by_wafer=by_wafer,
                               by_board=by_board, by_squid=by_squid)
        
        new_filter_function_map = core.G3MapVectorComplexDouble()
        
        for each_group, timestream_map in by_group_timestream_maps.items():
            if each_group not in filter_function_map.keys():
                 continue
            for each_bolo in timestream_map.keys():
                new_filter_function_map[each_bolo] = \
                    filter_function_map[each_group]
        
        filter_function_map = new_filter_function_map
        
        # Create a new filter function map with keys being bolometer names
        # out of the original filter function map whose keys are group names
        # so that the function fft_filter_mem_friendly_w_multi_filters
        # can be used.
    
    
    # - Apply the filters
    
    core.log_info(
        "\n", "   (Applying the filter functions...)",
        "\n", unit="notchfiltertimestreams")
    
    filtered_timestream_map = core.G3TimestreamMap()
    
    missing_bolos_in_filter_map = \
        set(original_timestream_map.keys()) - \
        set(filter_function_map.keys())
    
    filter_len = len(filter_function_map.values()[0])
    
    for bolo in missing_bolos_in_filter_map:
        filter_function_map[bolo] = np.array([1.0] * filter_len)
    
    fft_filter_mem_friendly_w_multi_filters(
        original_timestream_map,
        window_function,
        filter_function_map,
        filtered_timestream_map)
    
    
    return filtered_timestream_map





def make_tukey_window_w_turnarounds_tapered(n_data_pts_each_timestream):

    '''
    This function constructs a Tukey window for a long timestream
    that resulted from concatenating timestreams stored in several
    contiguous scan frames.
    
    If the number of timestreams that were concatenated is odd,
    then the first and last timestream are assumed to correspond to turnarounds,
    and the tapered parts of the window will be made to
    align with these two timestreams.
    
    
    Parameters:
    -----------
    
    n_data_pts_each_timestream: list
        A list of integers that represent the length of each of
        the timestreams that were concatenated.
    
    
    Returns:
    --------
    
    window_function: G3VectorDouble
        An array that represents the custom Tukey window.
    
    '''
    
    if (len(n_data_pts_each_timestream) % 2 == 0) or \
       (len(n_data_pts_each_timestream) == 1):
        total_n_data_pts = np.sum(n_data_pts_each_timestream)
        window_function  = signal.tukey(total_n_data_pts)
    else:
        window_function = np.array([])
        buffer_length   = np.min([n_data_pts_each_timestream[0],
                                  n_data_pts_each_timestream[-1]])
        total_length    = np.sum(n_data_pts_each_timestream)
        taper_fraction  = (2 * buffer_length) / total_length
        window_function = signal.tukey(total_length, alpha=taper_fraction)
        
        # If the input timestream map contains timestreams from
        # an odd number of scan frames, which it should
        # if the concatenation process in the FindLinesFromPSDs module
        # works properly, then the first and last scan frames
        # should correspond to turnarounds.
        # The window function used in this case is a Tukey window
        # with the two tapered cosine parts corresponding to
        # these two turnarounds.
    
    if window_function[0] == 0.0:
        window_function[0] = window_function[1]
    if window_function[-1] == 0.0:
        window_function[-1] = window_function[-2]
    
    # After certain frequency components are removed and
    # the inverse FFT taken, the timestream will be divided by
    # the window function, and avoiding division by zero
    # seems like a not-crazy idea.
    
    window_function = core.G3VectorDouble(window_function)
    
    return window_function





def split_timestreams(
        list_of_length_of_each_chunk,
        long_timestream_map):
    
    '''
    This function deconcatenates concatenated timestream maps
    and returns a list of timestream maps 
    whose lengths match those of the original individual chunks.
    
    
    Parameters:
    -----------
    
    list_of_length_of_each_chunk: list
        A list containing the length of each of the chunks into which
        a long timestream needs to be split.
    
    long_timestream_map: G3TimestreamMap
        The timestream map that contains timestreams to be split.
    
    
    Returns:
    --------
    
    deconcatenated_timestream_maps: list
        A list of G3TimestreamMaps that resulted from 
        splitting the long_timestream_map into multiple chunks.
    
    '''
    
    core.log_info(
        "\n", "- Splitting timestreams into the original chunks...",
        "\n", "   (The chunks are:", list_of_length_of_each_chunk, ")",
        "\n", unit="notchfiltertimestreams")
    
    list_of_indices_used_for_slicing = []
    deconcatenated_timestream_maps   = []
    
    for i in range(len(list_of_length_of_each_chunk)):
        start_index_for_this_chunk = \
            int(np.sum(list_of_length_of_each_chunk[0:i]))
        stop_index_for_this_chunk  = \
            int(np.sum(list_of_length_of_each_chunk[0:i+1]))
        list_of_indices_used_for_slicing.\
            append([start_index_for_this_chunk, stop_index_for_this_chunk])
        timestream_map_for_this_chunk = core.G3TimestreamMap()
        deconcatenated_timestream_maps.append(timestream_map_for_this_chunk)
    
    for each_bolo_name, its_concatenated_ts in long_timestream_map.items():
        for i, two_indices in enumerate(list_of_indices_used_for_slicing):
            start_index = two_indices[0]
            stop_index  = two_indices[1]
            deconcatenated_timestream_maps[i][each_bolo_name] = \
                its_concatenated_ts[start_index:stop_index]
    
    return deconcatenated_timestream_maps


