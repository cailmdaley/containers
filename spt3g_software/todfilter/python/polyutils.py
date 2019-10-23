import numpy as np
import copy

from spt3g import core
from spt3g.todfilter import poly_and_mhpf_filter_ts_data, make_empty_filter_mask, FilterMask
from spt3g.todfilter import get_mad_std, simple_common_mode_filter_cpp, masked_common_mode_filter_cpp
from spt3g.todfilter.util import mask_ts_map

from spt3g.calibration.template_groups import get_template_groups, get_template_groups_inverse
import copy

def poly_filter_g3_timestream_map( ts_map, poly_order):
    '''
    Applies a polynomial filter of order poly_order to the timestream map ts_map.
    '''

    return poly_and_mhpf_filter_ts_data(ts_map, make_empty_filter_mask(ts_map), 0., -1., poly_order);

def make_dynamic_mask_for_sources(ts_map, filter_order = 1, threshold = 5):
    '''
    creates a dynamic point source mask from the ts_map.

    It performs a poly filter of filter_order, estimates the noise, and then masks anything 'threshold' standard
    deviations in amplitude.

    Don't use this when estimating power spectra.

    '''
    pfilt = poly_filter_g3_timestream_map(ts_map, filter_order)
    dynomask = FilterMask()
    for k in pfilt.keys():
        ts = pfilt[k]
        ts /= get_mad_std(ts)
        mask = core.G3VectorInt( np.asarray( np.abs(ts) > threshold, dtype = 'int') )
        dynomask.pixel_mask[k] = mask
        has_masked = int(np.any(mask))
        dynomask.has_masked_pixels[k] = has_masked
    return dynomask

def DynamicSourceFiltering(frame, ts_map_key, out_mask_key, 
                           pre_filter_poly_order = 1, threshold = 5.0):
    '''
    Generates a dynamic point source mask.

    Poly filters the data with a poly of order: pre_filter_poly_order
    Calculates the standard deviation using a MAD estimator and adds any point
    "threshold" sigma away from the filtered timestream to the filter mask.

    '''

    if frame.type != core.G3FrameType.Scan:
        return
    mask = make_dynamic_mask_for_sources(frame[ts_map_key], pre_filter_poly_order, threshold)
    frame[out_mask_key] = mask
                                         
def ReplaceTsMapWithDynamicMask(frame,
                                in_dynamic_mask_key,
                                in_timestream_key,
                                out_ts_map_key):
    '''
    Fills in the point source mask into a timestream object.  
    This is entirely for debugging issues with the point source mask or dynamic point source masking
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    tsm = core.G3TimestreamMap()
    in_dynamic_mask = frame[in_dynamic_mask_key]
    in_ts_map = frame[in_timestream_key]

    for k in in_ts_map.keys():
        tsm[k] = core.G3Timestream(in_dynamic_mask.pixel_mask[k])
        tsm[k].start = in_ts_map[k].start
        tsm[k].stop = in_ts_map[k].stop

    frame[out_ts_map_key] = tsm

    
def simple_common_mode_filter(ts_map, ts_map_keys):
    '''
    applies a common mode filter over the keys ts_map_keys in ts_map.
    This applies the common mode filter in place, so ts_map is modified.
    '''
    if len(ts_map_keys) == 0:
        return
    ts_sum = copy.copy(ts_map[ts_map.keys()[0]]) * 0
    tsk = []
    for k in ts_map_keys:
        if k in ts_map:
            tsk.append(k)
    for k in tsk:
        ts_sum += ts_map[k]
    ts_sum /= float(len(tsk))
    for k in tsk:
        ts_map[k] = ts_map[k] - ts_sum

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def CommonModeFilterSlower(frame, in_ts_map_key, out_ts_map_key, 
                     per_band = True, per_wafer=False, per_squid = False, 
                     bolo_props=None, wiring_map=None):
    '''
    Applies a common mode (mean) filter.  It applies it over the template group specified.
    '''
    template_groups = get_template_groups(bolo_props, wiring_map,
                                          per_band, per_wafer, per_squid)
    ts_map = copy.copy(frame[in_ts_map_key])
    for tg in template_groups:
        bolos = [bolo for bolo in tg if bolo in ts_map.keys()]
        simple_common_mode_filter(ts_map, bolos)
    frame[out_ts_map_key] = ts_map


class CommonModeFilter(object):
    def __init__(self, in_ts_map_key, out_ts_map_key,
                 per_band = True, per_wafer=False, per_squid = False,
                 bolo_props_key = 'BolometerProperties',
                 wiring_map_key = 'WiringMap',
                 filter_mask_key = None):
        self.in_ts_map_key = in_ts_map_key
        self.out_ts_map_key = out_ts_map_key
        self.per_band = per_band
        self.per_wafer = per_wafer
        self.per_squid = per_squid
        self.bolo_props_key = bolo_props_key
        self.wiring_map_key = wiring_map_key
        self.filter_mask_key = filter_mask_key
        self.bolo_props = None
        self.wiring_map = None
        self.template_groups = None
    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            if frame.type == core.G3FrameType.Calibration:
                if self.bolo_props_key in frame:
                    self.bolo_props = frame[self.bolo_props_key]
            elif frame.type == core.G3FrameType.Wiring:
                self.wiring_map = frame[self.wiring_map_key]
            if ( not self.bolo_props is None and
                 not self.wiring_map is None):
                if any((self.per_wafer,self.per_band,self.per_squid)):
                    self.template_groups = get_template_groups(
                        self.bolo_props, self.wiring_map, self.per_band, 
                        self.per_wafer, self.per_squid)
                else:
                    self.template_groups = [[b for b in self.bolo_props.keys()]]
            return
        if len(frame[self.in_ts_map_key].keys()) == 0:
            return False
        ts_map = copy.copy(frame[self.in_ts_map_key])
        if self.filter_mask_key is not None:
            ps_map = frame[self.filter_mask_key].pixel_mask
            ts_map_masked = mask_ts_map(ts_map, ps_map)
            for tg in self.template_groups:
                masked_common_mode_filter_cpp(ts_map, ts_map_masked, tg)
        else:
            for tg in self.template_groups:
                simple_common_mode_filter_cpp(ts_map, tg)
        frame[self.out_ts_map_key] = ts_map

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def GenerateGroupTimestreams(frame, in_ts_map_key, out_ts_map_key, 
                             per_band = True, per_wafer=False, per_squid = False, 
                             bolo_props=None, wiring_map=None):
    '''
    Creates sum timestreams of all the detectors in the template groups specified.
    '''
    template_groups = get_template_groups(bolo_props, wiring_map,
                                          per_band, per_wafer, per_squid, include_keys = True)
    ts_map = frame[in_ts_map_key]
    out_ts_map = core.G3TimestreamMap()
    
    for tg_id, ts_lst in template_groups.items():
        sval = sum([ ts_map[bid] if bid in ts_map else ts_map[ts_map.keys()[0]]*0 for bid in ts_lst])
        out_ts_map[tg_id] = sval
    frame[out_ts_map_key] = out_ts_map



@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def GenerateGroupTimestreamsWeighted(frame, in_ts_map_key, var_key, 
                                     out_ts_map_key, 
                                     per_band = True, per_wafer=False, per_squid = False, 
                                     bolo_props=None, wiring_map=None):
    '''
    Creates weighted sum timestreams of all the detectors in the template groups specified.

    The weights are taken to be the inverse variance of the variances specified at 
    the var_key which maps to a G3MapDouble
    '''

    template_groups = get_template_groups(bolo_props, wiring_map,
                                          per_band, per_wafer, per_squid, include_keys = True)
    ts_map = frame[in_ts_map_key]
    variances = frame[var_key]

    out_ts_map = core.G3TimestreamMap()
    for tg_id, ts_lst in template_groups.items():
        filt_lst = []
        sval=sum([ ts_map[bid]/(variances[bid]) if bid in ts_map else 0 * ts_map[ts_map.keys()[0]] for bid in ts_lst])
        out_ts_map[tg_id] = sval
    frame[out_ts_map_key] = out_ts_map
    


@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def GenerateGroupTimestreamsAbs(frame, in_ts_map_key, out_ts_map_key, 
                                per_band = True, per_wafer=False, per_squid = False, 
                                bolo_props=None, wiring_map=None):
    '''
    Creates sum timestreams of all the detectors in the template groups specified.
    '''
    template_groups = get_template_groups(bolo_props, wiring_map,
                                          per_band, per_wafer, per_squid, include_keys = True)
    ts_map = frame[in_ts_map_key]
    out_ts_map = core.G3TimestreamMap()
    for tg_id, ts_lst in template_groups.items():
        sval = ts_map[ts_map.keys()[0]] * 0
        sval += sum([ np.abs(ts_map[bid]) if bid in ts_map else ts_map[ts_map.keys()[0]]*0 for bid in ts_lst])
        out_ts_map[tg_id] = sval
            
    frame[out_ts_map_key] = out_ts_map



@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def SplitTimestreamsByTemplateGroup(frame, in_ts_map_key, out_ts_map_key_root_name, 
                                    per_band = True, per_wafer=False, per_squid = False, 
                                    bolo_props=None, wiring_map=None):
    '''
    Untested as of yet
    '''
    template_groups = get_template_groups(bolo_props, wiring_map,
                                          per_band, per_wafer, per_squid, include_keys = True)
    ts_map = frame[in_ts_map_key]

    out_ts_map_dic = {}

    for tg_id, ts_lst in template_groups.items():
        k = out_ts_map_key_root_name + tg_id
        if not k in out_ts_map_dic:
            out_ts_map_dic[k] = core.G3TimestreamMap()
        for ts in ts_lst:
            if ts in ts_map:
                out_ts_map_dic[k][ts] = ts_map[k]
    for k in out_ts_map_dic:
        frame[out_ts_map_key] = out_ts_map
