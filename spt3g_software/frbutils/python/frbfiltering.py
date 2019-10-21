from __future__ import print_function
from spt3g import core, mapmaker, todfilter, coordinateutils
from spt3g.frbutils import G3VectorFrbEventInfo, G3VectorFrbDetInfo, get_poisson_threshold
from spt3g.frbutils import filter_frb_events
from copy import copy
from spt3g.sptpol import usefulmodules
from spt3g.calibration import template_groups


#true means it should be filteredj
import numpy as np
#generic utils


def count_elements_in_lst(lst):
    od = {}
    for k in lst:
        if not k in od:
            od[k] = 0
        od[k] += 1
    return od

def get_wonky_dets(ev_lst, prob_cutoff = 0.99, count_cutoff = None, count_scaling = 1):
    unique_dets = set()

    det_lst = []
    for ev in ev_lst:
        det_id = tuple(sorted((ev.det_info[0].bid,
                               ev.det_info[1].bid)))
        det_lst.append(det_id)
        unique_dets.add(det_id)

    n_active_dets = len(unique_dets)
    counted_dets = count_elements_in_lst(det_lst)
    total_count = sum(counted_dets.values())

    print('first', count_cutoff)
    if count_cutoff is None:
        count_cutoff = count_scaling * get_poisson_threshold(total_count, n_active_dets, prob_cutoff)
        
    print('total_count', total_count, 'n_active_dets', n_active_dets)
    print('count_cutoff', count_cutoff)
    bad_dets = []
    for k in counted_dets.keys():
        if counted_dets[k] > count_cutoff:
            bad_dets.append(k)
    return bad_dets, counted_dets


def filter_events(frb_events, is_bad_func):
    return filter_frb_events(frb_events, core.G3VectorInt(map(is_bad_func, frb_events)))

def apply_dev_filter_to_ev(ev, is_bad_func):
    ev.filter_det_evs(core.G3VectorInt(map(is_bad_func, ev.det_info)))

def FilterEvents(frame, input_frb_event_key, output_frb_event_key,
                 filter_func):
    if frame.type != core.G3FrameType.Scan:
        return
    frb_events = frame[input_frb_event_key]

    frame[output_frb_event_key] = filter_events(frb_events, filter_func)

def FilterDetectors(frame, input_frb_event_key, output_frb_event_key,
                    filter_func):
    if frame.type != core.G3FrameType.Scan:
        return
    frb_events = copy(frame[input_frb_event_key])
    for i, ev in enumerate(frb_events):
        apply_dev_filter_to_ev(ev, filter_func)
        frb_events[i] = ev
    frame[output_frb_event_key] = filter_events(frb_events, lambda ev: len(ev.det_info) == 0)

class FrbFilterAroundPointSources(object):
    def __init__(self, point_source_mask_file):
        self.point_source_mask = mapmaker.mapmakerutils.load_spt3g_map(point_source_mask_file)['T']
    def __call__(self, ev):
        for di in ev.det_info:
            ind = self.point_source_mask.angle_to_pixel(di.ra, di.dec)
            if (ind == -1 or ind >= np.size(self.point_source_mask)):
                continue
            val = self.point_source_mask[ind]
            if val:
                return True
        return False

#specific filters
class FrbFilterGlitchyDetectors(object):
    def __init__(self, input_frb_event_key, output_frb_event_key, 
                 n_tses_key,
                 unlikely_cutoff = 0.99  ):
        self.input_frb_event_key = input_frb_event_key
        self.output_frb_event_key = output_frb_event_key
        self.n_tses_key = n_tses_key
        self.unlikely_cutoff = unlikely_cutoff

        self.n_tses = []
        self.ev_count = {}
        self.frames = []
        self.has_data = False

        self.previous_obs_id = None
        self.new_obs = False

    def __call__(self, frame):
        if not self.has_data:
            self.n_tses = []
            self.frames = []
            self.ev_count = {}

        if frame.type == core.G3FrameType.Observation:
            obs_id = frame['ObservationID']
            if ( (not self.previous_obs_id is None) and
                 self.previous_obs_id != obs_id):
                self.new_obs = True
            self.previous_obs_id = obs_id

        if self.has_data and (frame.type == core.G3FrameType.EndProcessing or
                              (frame.type == core.G3FrameType.Observation and
                              self.new_obs)):

            self.new_obs = False
            self.has_data = False
            total_evs = sum(self.ev_count.values())
            num_dets = max(self.n_tses)
            cutoff = get_poisson_threshold(total_evs, num_dets*0.75, self.unlikely_cutoff)
            #cutoff = get_poisson_threshold(total_evs, num_dets, self.unlikely_cutoff)
            
            print("glitchy cutoff", cutoff, total_evs, num_dets)
            bad_dets = set()
            for k, v in self.ev_count.items():
                if v > cutoff:
                    bad_dets.add(k)
            for fr in self.frames:
                FilterEvents(fr, self.input_frb_event_key, self.output_frb_event_key,
                             lambda ev: any([di.bid in bad_dets for di in ev.det_info]))
            self.frames.append(frame)
            frames = self.frames
            self.frames = [] 
            self.n_tses = []
            self.ev_count = {}
            return frames

        if frame.type == core.G3FrameType.Scan:
            self.has_data = True 
            self.n_tses.append(frame[self.n_tses_key])
            self.frames.append(frame)
            evs = frame[self.input_frb_event_key]
            for ev in evs:
                for det_info in ev.det_info:
                    bid = det_info.bid
                    if not bid in self.ev_count:
                        self.ev_count[bid] = 0
                    self.ev_count[bid] += 1
            return []



class FrbFilterGlitchyDetectorsGoodDets(object):
    def __init__(self, input_frb_event_key, output_frb_event_key, 
                 tses_key,
                 unlikely_cutoff = 0.99  ):
        self.input_frb_event_key = input_frb_event_key
        self.output_frb_event_key = output_frb_event_key
        self.tses_key = tses_key
        self.unlikely_cutoff = unlikely_cutoff

        self.n_tses = []
        self.ev_count = {}
        self.frames = []
        self.has_data = False

        self.previous_obs_id = None
        self.new_obs = False

    def __call__(self, frame):
        if not self.has_data:
            self.n_tses = []
            self.frames = []
            self.ev_count = {}

        if frame.type == core.G3FrameType.Observation:
            obs_id = frame['ObservationID']
            if ( (not self.previous_obs_id is None) and
                 self.previous_obs_id != obs_id):
                self.new_obs = True
            self.previous_obs_id = obs_id

        if self.has_data and (frame.type == core.G3FrameType.EndProcessing or
                              (frame.type == core.G3FrameType.Observation and
                              self.new_obs)):

            self.new_obs = False
            self.has_data = False
            total_evs = sum(self.ev_count.values())
            num_dets = np.mean(self.n_tses)


            if self.unlikely_cutoff == 1:
                cutoff = total_evs
            elif total_evs > 0:
                cutoff = get_poisson_threshold(total_evs, num_dets, self.unlikely_cutoff)
            else:
                cutoff = 0
            #cutoff = get_poisson_threshold(total_evs, num_dets, self.unlikely_cutoff)
            
            print("glitchy cutoff", cutoff, total_evs, num_dets)
            bad_dets = set()
            for k, v in self.ev_count.items():
                if v > cutoff:
                    bad_dets.add(k)
            for fr in self.frames:
                FilterEvents(fr, self.input_frb_event_key, self.output_frb_event_key,
                             lambda ev: any([di.bid in bad_dets for di in ev.det_info]))
            self.frames.append(frame)
            frames = self.frames
            self.frames = [] 
            self.n_tses = []
            self.ev_count = {}
            return frames

        if frame.type == core.G3FrameType.Scan:
            self.has_data = True 
            tes_lst  = frame[self.tses_key]

            self.n_tses.append(len(tes_lst))
            self.frames.append(frame)
            evs = frame[self.input_frb_event_key]
            for ev in evs:
                for det_info in ev.det_info:
                    bid = det_info.bid
                    if not bid in tes_lst:
                        continue
                    if not bid in self.ev_count:
                        self.ev_count[bid] = 0
                    self.ev_count[bid] += 1
            return []


def det_filt_neg_amp(det_info):
    return det_info.amplitude <= 0

def det_filt_pos_amp(det_info):
    return det_info.amplitude >= 0


def ev_filt_neg_amp(ev):
    return any([di.amplitude <= 0 for di in ev.det_info])

def ev_filt_pos_amp(ev):
    return any([di.amplitude >= 0 for di in ev.det_info])

def ev_filt_same_ind(ev):
    return ev.det_info[0].scan_index != ev.det_info[1].scan_index

def ev_filt_one_ind(ev):
    return abs(ev.det_info[0].scan_index - ev.det_info[1].scan_index) > 1


def ev_filt_bad_heavi(ev):
    return any(map(lambda di: di.heavi_amp > di.amplitude, ev.det_info))


#lower percent_of_amp means stricter filter
def get_ev_filter_bad_heaviside(percent_of_amp = 1.0): 
    return lambda ev: any(map(lambda di:abs(di.heavi_amp) > abs(percent_of_amp * di.amplitude), 
                              ev.det_info))


def get_ev_filter_bad_phase(amp_ratio = 0.2): 
    return lambda ev: any(map(lambda di: abs(di.q_amplitude/di.amplitude) > amp_ratio, 
                              ev.det_info))


def get_ev_filt_percent_above_5_sig(percent_above):
    def helper(ev):
        for di in ev.det_info:
            if di.n_samps_above_sig_5/float(ev.scan_len) > percent_above:
                return True
        return False
    return helper

def FilterAbove5Sig(frame, event_key_in, event_key_out, 
                    percent_above = 0.002):
    FilterEvents(frame, event_key_in, event_key_out, 
                 get_ev_filt_percent_above_5_sig(percent_above))

def FilterHeaviside(frame, event_key_in, event_key_out, percent_of_amp = 0.5):
    FilterEvents(frame, event_key_in, event_key_out, 
                 get_ev_filter_bad_heaviside(percent_of_amp))


#two detectors
def ev_filt_two_events(ev):
    return len(ev.det_info) != 2

def ev_filt_four_events(ev):
    return len(ev.det_info) != 4

#group filter
def get_ev_filter_squid_sig(max_cutoff, mean_cutoff):
    def h(ev):
        return all(map(lambda di: di.squid_sig > mean_cutoff or 
                       di.squid_sig_max > max_cutoff, ev.det_info))
    return h

def get_ev_filter_squid_live_num(min_num_live):
    return lambda ev: all(map(lambda di: di.n_live_on_squid < min_num_live, ev.det_info))

def get_ev_filter_wafer_sig(cutoff):
    return lambda ev: all(map(lambda di: di.wafer_sig > cutoff, ev.det_info))

def get_bad_pixel_filter( pix_lst):
    def helper(ev):
        pixel_pair = tuple(sorted( [ev.det_info[0].bid, ev.det_info[1].bid]))
        return pixel_pair in pix_lst
    return helper

def FrbFilterBadPixels(frame, pix_lst, event_key_in, event_key_out  ):
    FilterEvents(frame, event_key_in, event_key_out,
                 get_bad_pixel_filter(pix_lst))




@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def DetFiltSharedWafer(frame, input_frb_event_key, output_frb_event_key,
                       bolo_props = None):
    assert( not bolo_props is None )
    frb_events = copy(frame[input_frb_event_key])
    for i, ev in enumerate(frb_events):
        #find the highest sig event
        sig = -1e6
        bid = None
        for di in ev.det_info:
            if di.significance > sig:
                sig = di.significance
                bid = di.bid
        waf_id = bolo_props[bid].wafer_id
        filter_func = lambda di: bolo_props[di.bid].wafer_id != waf_id
        apply_dev_filter_to_ev(ev, filter_func)
        frb_events[i] = ev
    frame[output_frb_event_key] = filter_events(frb_events, lambda ev: len(ev.det_info) == 0)


def get_two_ev_sig_filter(min_sig, max_sig, max_small_sig = None):
    def filt_func_inner(ev):
        sigs = sorted([di.significance for di in ev.det_info])
        assert(len(sigs) == 2)
        return ( (sigs[0] < min_sig) or
                 (sigs[1] > max_sig) or
                 ( (not max_small_sig is None) and
                   (sigs[0] > max_small_sig) ) )
    return filt_func_inner

def FilterSig(frame, event_key_in, event_key_out,
              min_sig, max_sig, max_small_sig = None):
    if frame.type != core.G3FrameType.Scan:
        return
    if (not event_key_in in frame):
        return []
    FilterEvents(frame, event_key_in, event_key_out,
                 get_two_ev_sig_filter(min_sig = min_sig, max_sig = max_sig, 
                                       max_small_sig = max_small_sig))

def FilterToPixelPair(frame, event_key_in, event_key_out):
    FilterEvents(frame, event_key_in, event_key_out,
                 lambda ev: ev.det_info[0].pixel_id != ev.det_info[1].pixel_id)

def FilterBolosList(frame, event_key_in, event_key_out, bad_bolos = None):
    def hfunc(ev):
        for di in ev.det_info:
            if di.bid in bad_bolos:
                return True
        return False
    assert(not bad_bolos is None)
    FilterEvents(frame, event_key_in, event_key_out, hfunc)


def FilterSquidSig(frame, event_key_in, event_key_out,
                   sig):
    if frame.type != core.G3FrameType.Scan:
        return

    if not event_key_in in frame:
        return []
    FilterEvents( frame, 
                  input_frb_event_key = event_key_in,
                  output_frb_event_key = event_key_out,
                  filter_func = get_ev_filter_squid_sig( cutoff = sig)
             )

def get_two_ev_var_filter(max_var, min_var = 0):
    def helper(ev):
        for di in ev.det_info:
            if di.variance > max_var or di.variance < min_var:
                return True
        return False
    return helper
        
def FilterVar(frame, event_key_in, event_key_out, max_var, min_var = 0):
    FilterEvents( frame, 
                  input_frb_event_key = event_key_in,
                  output_frb_event_key = event_key_out,
                  filter_func = get_two_ev_var_filter(max_var, min_var)
                  )

def get_amp_filter_two_det(min_amp, max_amp, min_large_amp = 0):
    def helper(ev):
        min_det_amp = min( ( ev.det_info[0].amplitude, ev.det_info[1].amplitude))
        max_det_amp = max( ( ev.det_info[0].amplitude, ev.det_info[1].amplitude))
        return ( (min_det_amp < min_amp) or (max_det_amp > max_amp) or (max_det_amp < min_large_amp) )
    return helper

def FilterAmpTwoEv(frame, event_key_in, event_key_out, min_amp, max_amp, min_large_amp = 0 ):
    FilterEvents( frame,
                  input_frb_event_key = event_key_in,
                  output_frb_event_key = event_key_out,
                  filter_func = get_amp_filter_two_det(min_amp, max_amp, min_large_amp)
             )

def FilterToSharedSquid(frame, event_key_in, event_key_out):
    FilterEvents(frame, event_key_in, event_key_out,
                 lambda ev: (ev.det_info[0].module_id != ev.det_info[1].module_id or
                             ev.det_info[0].board_id != ev.det_info[1].board_id)
                 )


def FilterToSharedWafer(frame, event_key_in, event_key_out):
    FilterEvents(frame, event_key_in, event_key_out,
                 lambda ev: (ev.det_info[0].bid[:2] != ev.det_info[1].bid[:2])
                 )


def FilterShitDetectors(frame, 
                        event_key_in, 
                        event_key_out, 
                        percent_above_sig_5, 
                        cut_chunk = 20,
                        model_ll_key = 'ModelLogLike', 
                        baseline_ll_key = 'BaselineLogLike'):
    if frame.type != core.G3FrameType.Scan:
        return
    model_ll = frame[model_ll_key]
    base_ll = frame[baseline_ll_key]
    print(frame)
    evs = frame[event_key_in]
    is_bad_vec = core.G3VectorInt(np.zeros(len(evs), dtype = 'int'))
    for i, ev in enumerate(evs):
        is_bad = 0
        start_ind = max(ev.scan_index - cut_chunk, 0)
        stop_ind =  min(ev.scan_index + cut_chunk, ev.scan_len)
        for di in ev.det_info:
            sig = np.abs(np.array(model_ll[di.bid] - base_ll[di.bid]))
            sig[start_ind:stop_ind] = 0
            percent_bad = float(np.size(np.where(sig > 5))) / float(ev.scan_len)
            if percent_bad > percent_above_sig_5:
                is_bad = 1
                break
        is_bad_vec[i] = is_bad
    frame[event_key_out] = filter_frb_events(evs, is_bad_vec)

@core.pipesegment
def SptpolFrbFiltering(pipe, 
                       input_frb_event_key = 'FrbEvents',
                       n_tses_key = 'NumLiveDetectors',
                       glitchy_prob_cutoff = 0.97,
                       point_source_file = None,
                       del_big_ev_lists = False,
                       filter_to_waf = True ):
                          

    prev_frb_key = input_frb_event_key
    keys_to_del = []

    current_frb_key = 'FrbEvs'
    if filter_to_waf:
        keys_to_del.append(prev_frb_key)
        current_frb_key += 'Waf'
        pipe.Add(DetFiltSharedWafer, 
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key) 
        prev_frb_key = current_frb_key
        keys_to_del.append(prev_frb_key)

    current_frb_key += 'Gl'
    pipe.Add(FrbFilterGlitchyDetectors, 
             input_frb_event_key = prev_frb_key, 
             output_frb_event_key = current_frb_key,  
             n_tses_key = n_tses_key,
             unlikely_cutoff = glitchy_prob_cutoff)
    prev_frb_key = current_frb_key
    keys_to_del.append(prev_frb_key)

    if not point_source_file is None:
        keys_to_del.append(prev_frb_key)
        current_frb_key += 'Pnt'
        pipe.Add(FilterEvents, 
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = FrbFilterAroundPointSources(point_source_file)
             )
        prev_frb_key = current_frb_key
    if del_big_ev_lists:
        pipe.Add(core.Delete, keys = keys_to_del,
                 type = core.G3FrameType.Scan)

@core.pipesegment
def SptpolFrbFilteringSecond(pipe, 
                             input_frb_event_key = 'FrbEvsGlPnt',
                             co_sig = None,
                             is_neg = False,
                             squid_max_cutoff = 4, 
                             squid_mean_cutoff = 4, 
                             wafer_ll_cutoff = 4,
                             heavi_percent_of_amp = 1.0,
                             filt_same = True,
                             filt_evs_of_sign = False,
                             point_source_file = None,
                             include_four_cut = False,
                             phase_amp = None,
                             del_big_ev_lists = False):

    prev_frb_key = input_frb_event_key
    current_frb_key = input_frb_event_key
    keys_to_del = []
    keys_to_del.append(prev_frb_key)

    if not co_sig is None:
        current_frb_key += 'CoSig'
        pipe.Add(FilterDetectors, 
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = lambda det_info: det_info.significance < co_sig
                 )
        prev_frb_key = current_frb_key
        keys_to_del.append(prev_frb_key)


    keys_to_del.append(prev_frb_key)
    current_frb_key += 'Waf'
    pipe.Add(DetFiltSharedWafer, 
             input_frb_event_key = prev_frb_key, 
             output_frb_event_key = current_frb_key) 
    prev_frb_key = current_frb_key
    keys_to_del.append(prev_frb_key)



    current_frb_key += 'Hev'
    pipe.Add(FilterEvents, 
             input_frb_event_key = prev_frb_key, 
             output_frb_event_key = current_frb_key,  
             filter_func = get_ev_filter_bad_heaviside(percent_of_amp = heavi_percent_of_amp)
    )
    prev_frb_key = current_frb_key
    keys_to_del.append(prev_frb_key)


    if not filt_evs_of_sign:
        if is_neg:
            current_frb_key += 'Pos'
            pipe.Add(FilterDetectors, 
                     input_frb_event_key = prev_frb_key, 
                     output_frb_event_key = current_frb_key,  
                     filter_func = det_filt_pos_amp )
            prev_frb_key = current_frb_key
            keys_to_del.append(prev_frb_key)
        else:
            current_frb_key += 'Neg'
            pipe.Add(FilterDetectors, 
                     input_frb_event_key = prev_frb_key, 
                     output_frb_event_key = current_frb_key,  
                     filter_func = det_filt_neg_amp
                     )
            prev_frb_key = current_frb_key
            keys_to_del.append(prev_frb_key)
    else:
        if is_neg:
            current_frb_key += 'PosEv'
            pipe.Add(FilterEvents, 
                     input_frb_event_key = prev_frb_key, 
                     output_frb_event_key = current_frb_key,  
                     filter_func = ev_filt_pos_amp
                     )
            prev_frb_key = current_frb_key
            keys_to_del.append(prev_frb_key)
        else:
            current_frb_key += 'NegEv'
            pipe.Add(FilterEvents, 
                     input_frb_event_key = prev_frb_key, 
                     output_frb_event_key = current_frb_key,  
                     filter_func = ev_filt_neg_amp
                     )
            prev_frb_key = current_frb_key
            keys_to_del.append(prev_frb_key)

    if not include_four_cut:
        current_frb_key += 'Two'
        pipe.Add(FilterEvents, input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = ev_filt_two_events
                 )
        prev_frb_key = current_frb_key
        keys_to_del.append(prev_frb_key)
    else:
        current_frb_key += 'Four'
        pipe.Add(FilterEvents, input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = ev_filt_four_events
                 )
        prev_frb_key = current_frb_key
        keys_to_del.append(prev_frb_key)        

    if not point_source_file is None:
        current_frb_key += 'Pnt'
        pipe.Add(FilterEvents, 
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = FrbFilterAroundPointSources(point_source_file)
             )
        prev_frb_key = current_frb_key
        keys_to_del.append(prev_frb_key)
        

    current_frb_key += 'Sqs'
    pipe.Add(FilterEvents, input_frb_event_key = prev_frb_key, 
             output_frb_event_key = current_frb_key,  
             filter_func = get_ev_filter_squid_sig( max_cutoff = squid_max_cutoff, 
                                                    mean_cutoff = squid_mean_cutoff  )
         )
    prev_frb_key = current_frb_key

    current_frb_key += 'Waf'
    pipe.Add(FilterEvents, input_frb_event_key = prev_frb_key, 
             output_frb_event_key = current_frb_key,  
             filter_func = get_ev_filter_wafer_sig( cutoff = wafer_ll_cutoff)
    )
    prev_frb_key = current_frb_key


    if filt_same:
        current_frb_key += 'IndS'
        pipe.Add(FilterEvents,
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = ev_filt_same_ind)
        prev_frb_key = current_frb_key
    else:
        current_frb_key += 'IndT'
        pipe.Add(FilterEvents,
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = ev_filt_one_ind)
        prev_frb_key = current_frb_key

    if not phase_amp is None:
        current_frb_key += 'Pha'
        pipe.Add(FilterEvents, 
                 input_frb_event_key = prev_frb_key, 
                 output_frb_event_key = current_frb_key,  
                 filter_func = get_ev_filter_bad_phase(amp_ratio = phase_amp)
                 )
        prev_frb_key = current_frb_key

    pipe.Add(core.Rename, keys = {prev_frb_key:'FrbEventsFilteringOut'}, 
             type = core.G3FrameType.Scan)

    if del_big_ev_lists:
        pipe.Add(core.Delete, keys = keys_to_del,
                 type = core.G3FrameType.Scan)

def make_sptpol_deep_point_source_mask(pntsrc_file, out_fn):
    flat_map = core.G3Map
    mapmaker.pointsourceutils.make_point_source_map(flat_map, pntsrc_file, 
                                                    mask_oob_pixels = False)




class ConstructValidBidsList:
    def __init__(self, 
                 input_ts_lst = 'DetectorIdList',
                 output_key = 'GoodChannels',

                 variance_key = 'TsVariance',
                 ts_len_key = 'TimestreamLength',
                 n_samps_above_5_sig_key = 'NSampsAbove5Sig',

                 var_min = 0,
                 var_max = 1e12,
                 
                 percent_above_5_sig = 0.002,
                 
                 bad_channels = None,
                 bolo_props_key = 'BolometerProperties', 
                 wiring_map_key = 'WiringMap',
                 cal_response_key = 'CalibratorResponseSN',
                 el_nod_sn_key = 'ElnodSNSlopes',
                 
                 ):
        self.input_ts_lst = input_ts_lst
        self.output_key = output_key
        self.variance_key = variance_key
        self.ts_len_key = ts_len_key
        self.n_samps_above_5_sig_key = n_samps_above_5_sig_key
        self.var_min = var_min
        self.var_max = var_max
        self.percent_above_5_sig = percent_above_5_sig
        self.bad_channels = bad_channels
        self.bolo_props_key = bolo_props_key
        self.wiring_map_key = wiring_map_key
        self.cal_response_key = cal_response_key
        self.el_nod_sn_key = el_nod_sn_key

        self.bias_props = None
        self.wiring_map = None

        self.calsn = None
        self.elnod_sn = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame['WiringMap']

        if frame.type == core.G3FrameType.Calibration:
            self.bias_props = frame[self.bolo_props_key]
            self.calsn = frame[self.cal_response_key]
            self.elnod_sn = frame[self.el_nod_sn_key]
        if ((frame.type == core.G3FrameType.Wiring or 
             frame.type == core.G3FrameType.Calibration) and
            (not self.bias_props is None) and (not self.wiring_map is None) ):

            self.pixel_to_bolo_map = template_groups.get_template_groups(
                self.bias_props, per_band = False, per_pixel = True, include_keys = True)
            self.bolo_to_pixel_map = template_groups.get_template_groups_inverse(
                self.bias_props, per_band = False, per_pixel = True)

        if frame.type != core.G3FrameType.Scan:
            return
        
        bids = frame[self.input_ts_lst]

        variance = frame[self.variance_key]
        ts_len = frame[self.ts_len_key]
        n_samps_above_5_sig = frame[self.n_samps_above_5_sig_key]
        
        if self.bad_channels is None:
            bad_bids = set()
        else:
            bad_bids = set(self.bad_channels)

        for bid in bids:
            physical_id = self.bias_props[bid].physical_name
            if ( variance[bid] < self.var_min 
                 or variance[bid] > self.var_max
                 or float(n_samps_above_5_sig[bid])/ts_len > self.percent_above_5_sig
                 or physical_id == ''
                 or (physical_id[-1] != 'Y' and physical_id[-1] != 'X') 
                 or (physical_id[0] == 'F')
                 or (physical_id[0] == 'f')
                 or self.calsn[bid] < 20
                 or self.elnod_sn[bid] < 20
                 ):
                bad_bids.add(bid)

        tmp_bad_bids = set(bad_bids)
        for bid in tmp_bad_bids:
            pix = self.bolo_to_pixel_map[bid]
            for pix_bid in self.pixel_to_bolo_map[pix]:
                bad_bids.add(pix_bid)

        bids = set(bids)
        out_bids = []
        for bid in bids:
            if not bid in bad_bids:
                out_bids.append(bid)
        frame[self.output_key] = core.G3VectorString(out_bids)
