from spt3g import core, todfilter, calibration, coordinateutils, gcp, dfmux
from spt3g import mapmaker, todfilter, std_processing

from spt3g.frbutils import PolyLikelihoodFiller, DeltaEventHunter, FrbEventInfo, FrbDetInfo
from spt3g.frbutils import get_ts_poly_delta_heavi_ll
from spt3g.frbutils.impulseinjections import InjectFrbSignal, FastTransientSignalBeam
from spt3g.frbutils.impulseinjections import LookForInjectedEvents
from spt3g.timestreamflagging import FlagBadG3MapValue, RemoveFlaggedTimestreams
from spt3g.calibration.template_groups import get_template_groups, get_template_groups_inverse

from spt3g.sptpol import usefulmodules
from spt3g.sptpol import sptpoltime

from spt3g.sptpol.directidfreader import DirectIdfReader, MultipleIdfReader

from spt3g.util import genericutils as GU

from spt3g.frbutils.frbmontecarlo import CreateFakeFrbTimestream, read_xtalk_file

from spt3g.timestreamflagging import FlagNaNs, FlagTimestreamsWithoutProperties, FlagBadHousekeeping, FlagToBand, FlagMissing

#used to find where I am getting numpy warnings
#from spt3g.util.genericutils import stdTracer
#import sys
#sys.stdout = stdTracer()
#sys.stderr = stdTracer()

import os.path, copy, random
try:
    import cPickle as pickle
except ImportError:
    import pickle
import numpy as np

import operator

PROFILE_FRB_CODE = False
'''
Vectors return by value (not reference) when accessed.  So modifying the contents of vectors of 
FrbEventInfo/FrbDetInfo in place is not possible.  You need to modify the individual elements
and then assign the modified structure to the vector.
'''

'''
Cosmic ray ldf will come from the one file
Need to modify the find the injected event and cache it in the frame

detector characterization?  
detector count histogram? can get after the fact.
counts as a func of time?

'''


def fc_add_vec_doub(a,b):
    return core.G3VectorDouble(np.asarray(a) + np.asarray(b))
def fc_add(a,b):
    return a + b
def fc_ident(a,b):
    return a
def fc_extend(a,b):
    a_cpy = copy.copy(a)
    a_cpy.extend(b)
    return a_cpy
def fc_add_int(a,b):
    return core.G3Int(core.G3Int(a).value + core.G3Int(b).value)
def fc_add_dic(a,b):
    b = copy.copy(b)
    for k in a.keys():
        if k in b:
            b[k] += a[k]
        else:
            b[k] = a[k]
    return b

class FrameCombiner(object):
    def __init__(self, type, key_op_mapping):
        self.frame_type = type
        self.output_frame = core.G3Frame(type)
        self.kop = key_op_mapping
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            return [self.output_frame, frame]
        if frame.type != self.frame_type:
            return
        for k, op in self.kop.items():
            if not k in frame:
                #print('not in')
                continue
            if not k in self.output_frame:
                try:
                    self.output_frame[k] = frame[k]
                except:
                    if GU.is_int(frame[k]):
                        self.output_frame[k] = core.G3Int(frame[k])
            else:
                oframe_val = self.output_frame[k]
                del self.output_frame[k]
                try:
                    self.output_frame[k] = op(oframe_val, frame[k])
                except:
                    print(type(oframe_val), oframe_val)
                    print(k)
                    raise RuntimeError()
        return []

class SptpolFrbHuntFrameCombiner(object):
    def __init__(self):
        self.fc = FrameCombiner(type = core.G3FrameType.Scan,
                                key_op_mapping = {

                'FrbEvsWaf': fc_extend,

                'FrbEvsHevNegGlTwo': fc_extend,
                'FrbEvsHevNegGlTwoSqn': fc_extend,
                'FrbEvsHevNegGlTwoSqnSqs': fc_extend,
                'FrbEvsHevNegGlTwoSqnSqsWaf': fc_extend,
                
                'FrbEvsHevNegGlTwo': fc_extend,
                'FrbEvsHevNegGlTwoPnt': fc_extend,
                'FrbEvsHevNegGlTwoPntSqn': fc_extend,
                'FrbEvsHevNegGlTwoPntSqnSqs': fc_extend,
                'FrbEvsHevNegGlTwoPntSqnSqsWaf': fc_extend,
                'FrbEventsFilteringOut': fc_extend,
                'FrbEventsFilteringOutWithBids': fc_extend,
                'FrbEvsHevNegGlTwoPntSqnSqsWafWithBids': fc_extend,
                
                'HistOfLoglikes':fc_add_vec_doub,
                'HistOfLoglikesBins':fc_ident,
                'HistOfTimestreams': fc_add_vec_doub,
                'HistOfTimestreamsBins':fc_ident,
                'EffectiveDetectorLiveSamps':fc_add_int,
                'Unity':fc_add,
                'DetFrbEventCount': fc_add_dic,
                
                'FullFilteredFrbs': fc_extend, #for later analysis
                                })
    def __call__(self, frame):
        return self.fc(frame)

class AddScanNumber(object):
    def __init__(self):
        self.sn = 0
    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return
        frame['ScanNumber'] = self.sn
        self.sn += 1


def add_frb_files(file_list, out_fn):
    assert(len(file_list) > 0)
    
    f = core.G3File(file_list[0])
    cal_frame = None
    wiring_frame = None
    fc = SptpolFrbHuntFrameCombiner()
    for frame in f:
        if frame.type == core.G3FrameType.Calibration:
            cal_frame = frame
        elif frame.type == core.G3FrameType.Wiring:
            wiring_frame = frame
        fc(frame)
    for fn in file_list[1:]:
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename = fn)
        #pipe.Add(core.Dump)
        pipe.Add(fc)
        pipe.Run()
    scan_frame = fc.fc.output_frame
    
    frames = [wiring_frame, cal_frame, scan_frame]
    writer = core.G3Writer(filename = out_fn)
    for fr in frames:
        writer(fr)



def AddPixelIdSptpol(frame, bkey = 'BolometerProperties'):
    if frame.type == core.G3FrameType.Calibration:
        bprops = copy.copy(frame[bkey])
        for k in bprops.keys():
            phys_name = bprops[k].physical_name
            bprops[k].pixel_id = phys_name[:phys_name.rfind('.')]
            bprops[k].wafer_id = phys_name[:phys_name.find('.')]
        del frame[bkey]
        frame[bkey] = bprops


def plot_tes_spatial(xs,ys, vals, circle_radius = .1, color_map = 'rainbow'):
    import matplotlib.pyplot as plt
    import copy
    vals = copy.copy(vals)
    vmax =  max(vals.values())
    xmax = max(xs.values())
    xmin = min(xs.values())
    ymax = max(ys.values())
    ymin = min(ys.values())

    for k in vals:
        vals[k] /= vmax

    cmap = plt.get_cmap(color_map)
    fig = plt.gcf()


    for k in vals.keys():
        if not k in xs or not k in ys:
            continue
        circle = plt.Circle((xs[k],ys[k]),circle_radius, color = cmap(vals[k]))
        fig.gca().add_artist(circle)
    fig.axes[0].set_xlim((xmin-circle_radius*2,xmax+circle_radius*2))
    fig.axes[0].set_ylim((ymin-circle_radius*2,ymax+circle_radius*2))    



def plot_seasonal_time_thing( dmap ):
    import sptpoltime
    import pylab as pl

    pl.clf()
    for k in dmap:
        mjd = sptpoltime.datetime_to_mjd(sptpoltime.SptDatetime(k))
        pl.plot(mjd, dmap[k], '*')


def AddBidCount(frame, ev_key = None, out_key = "DetFrbEventCount"):
    if frame.type != core.G3FrameType.Scan:
        return
    assert(not ev_key is None and not out_key is None)
    evs = frame[ev_key]
    od = core.G3MapInt()
    for ev in evs:
        for di in ev.det_info:
            od[di.bid] = 1 + od.get(di.bid, 0)
    frame[out_key] = od


#sptpoltime.datetime_to_mjd(sptpoltime.SptDatetime('20130422_003228'))
def AddEventLiveTimeStats(frame, id_ts_map_key,
                          liveness_ts_map_key,
                          ts_len_key, det_number_key, det_id_lst_key,
                          detector_livetime_stats_key
):
    if frame.type == core.G3FrameType.Scan:
        ts_map = frame[liveness_ts_map_key]
        ndets = len(ts_map)
        if (ndets > 0):
            tslen = len(ts_map[ts_map.keys()[0]])
        else:
            tslen = 0
        frame[ts_len_key] = int(tslen)
        frame[det_number_key] = int(ndets)
        frame[detector_livetime_stats_key] = int(tslen) * int(ndets)

        ts_lst = frame[id_ts_map_key].keys()
        frame[det_id_lst_key] = core.G3VectorString(ts_lst)

class AddTesPosition(object):
    def __init__(self, pkl_file):
        d = pickle.load(open(pkl_file))
        self.om_x = core.G3MapDouble()
        self.om_y = core.G3MapDouble()
        for l in d:
            self.om_x[l[0]] = float(l[1])
            self.om_y[l[0]] = float(l[2])
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            frame['TesPosX'] = self.om_x
            frame['TesPosY'] = self.om_y

def get_histogram_bins(n_bins, min_bin, max_bin):
    return [ i * (max_bin-min_bin) / float(n_bins) + min_bin for i in range(n_bins+1)]

def HistogramLoglikelihoods(frame, model_key, baseline_key, bins, output_key, bin_key):
    core.log_debug('entering HistogramLoglikelihoods')
    if frame.type == core.G3FrameType.Scan:
        assert(model_key in frame)
        assert(baseline_key in frame)
        assert(output_key not in frame)
        assert(bin_key not in frame)
        output_histogram =  np.zeros( len(bins) - 1)
        core.log_debug('getting vals')
        mll_map = frame[model_key]
        bll_map = frame[baseline_key]
        for k in mll_map.keys():
            output_histogram +=  np.histogram(mll_map[k] - bll_map[k], bins = bins)[0]
        core.log_debug('storing vals')
        frame[output_key] = core.G3VectorDouble(output_histogram)
        frame[bin_key] = core.G3VectorDouble(bins)

def HistogramTimestreams(frame, ts_key, bins, output_key, bin_key):
    if frame.type == core.G3FrameType.Scan:
        output_histogram =  np.zeros( len(bins) - 1)
        ts_map = frame[ts_key]
        for k in ts_map.keys():
            output_histogram +=  np.histogram(ts_map[k], bins = bins)[0]
        frame[output_key] = core.G3VectorDouble(output_histogram)
        frame[bin_key] = core.G3VectorDouble(bins)

def FilterBumScans(frame, ts_key):
    if frame.type == core.G3FrameType.Scan:
        tss = frame[ts_key]
        if len(tss[tss.keys()[0]]) < 20:
            return None
        else:
            return frame

def FilterBadScan(frame, bad_scan_key):
    if frame.type != core.G3FrameType.Scan:
        return
    if bad_scan_key in frame and frame[bad_scan_key]:
        return []
    else:
        return

'''
things we need:
board_id, module_id: bprops wiringmap
ra, dec: azel calc
n_live_on_squid: ts_map, inverse_val
squid_sig: bprops ts_map
wafer_sig: bprops ts_map
'''


class AddSquidLiveTimes(object):
    def __init__(self, bolo_props_key, wiring_map_key, ts_map_key, squid_liveness_key, 
                 bolo_squid_liveness_key ):
        self.bolo_props_key = bolo_props_key
        self.wiring_map_key = wiring_map_key
        self.ts_map_key = ts_map_key
        self.squid_liveness_key = squid_liveness_key

        self.bolo_props = None
        self.wiring_map = None
        
        self.squid_template_group = None
        self.squid_template_group_inv = None
        self.bolo_squid_liveness_key = bolo_squid_liveness_key
        
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if self.squid_template_group == None:
                self.squid_template_group = get_template_groups(
                    self.bolo_props, self.wiring_map,
                    per_band = False, per_wafer=True, per_squid = True,
                    include_keys = True)
            squid_live_count = core.G3MapInt()
            ts_map = frame[self.ts_map_key]

            bolo_squid_live_count = core.G3MapDouble()
            for sq, bids in self.squid_template_group.items():
                nlive = 0
                for bid in bids:
                    if bid in ts_map:
                        nlive += 1
                squid_live_count[sq] = nlive
            frame[self.squid_liveness_key] = squid_live_count
            for sq, bids in self.squid_template_group.items():
                nlive = 0
                for bid in bids:
                    if bid in ts_map:
                        bolo_squid_live_count[bid] = squid_live_count[sq]
            frame[self.bolo_squid_liveness_key] = bolo_squid_live_count

        elif frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame[self.wiring_map_key]
        elif frame.type == core.G3FrameType.Calibration:
            self.bolo_props = frame[self.bolo_props_key]




def AddSampleRate(frame, sample_rate_key, ts_map_key):
    if frame.type != core.G3FrameType.Scan:
        return
    tses = frame[ts_map_key]
    frame[sample_rate_key] = tses[tses.keys()[0]].sample_rate

def AddUnityKey(frame, unity_key = 'Unity'):
    frame[unity_key] = 1

def AddNSampsAboveSig(
        frame,
        thresh, output_key,
        model_ll_key = 'ModelLogLike', 
        baseline_ll_key = 'BaselineLogLike'):
    if frame.type != core.G3FrameType.Scan:
        return
    out_map = core.G3MapInt()
    model_ll = frame[model_ll_key]
    base_ll = frame[baseline_ll_key]
    for k in model_ll.keys():
        mod = model_ll[k]
        base = base_ll[k]
        sig = np.abs(mod-base)
        n_above = np.size(np.where(sig > thresh))
        out_map[k] = n_above
    frame[output_key] = out_map


@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def AddFrbSquidSig(frame, ll_model_key, ll_baseline_key,
                   frb_events_key, bolo_props = None,wiring_map = None):
    squid_template_group = get_template_groups(
        bolo_props, wiring_map,
        per_band = False, per_wafer=True, per_squid = True,
        include_keys = True)
    squid_template_group_inv = get_template_groups_inverse(
        bolo_props, wiring_map,
        per_band = False, per_wafer=True, per_squid = True)
    frb_events = copy.copy(frame[frb_events_key])
    ll_model = frame[ll_model_key]
    ll_baseline = frame[ll_baseline_key]
    for i, ev in enumerate(frb_events):
        ev_bids = map(lambda dinf: dinf.bid, ev.det_info)
        for j, dinf in enumerate(ev.det_info):
            bid = dinf.bid
            squid = squid_template_group_inv[bid]
            sigs = []
            #sq_bids = []
            for sq_bid in squid_template_group[squid]:
                if not sq_bid in ev_bids and sq_bid in ll_model:
                    #sq_bids.append(sq_bid)
                    sigs.append(abs(ll_model[sq_bid][dinf.scan_index] - 
                                    ll_baseline[sq_bid][dinf.scan_index]))
            if len(sigs) > 0:
                dinf.squid_sig = np.mean(sigs)
                dinf.squid_sig_max = np.max(sigs)
            else:
                dinf.squid_sig = 0
                dinf.squid_sig_max = 0
            ev.det_info[j] = dinf
        #print( ev_bids, sq_bids, sigs, dinf.squid_sig, dinf.squid_sig_max)
        frb_events[i] = ev
    del frame[frb_events_key]
    frame[frb_events_key] = frb_events

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap')
def AddFrbWaferSig(frame,
                        ts_map_key, 
                        var_key,
                        frb_events_key, 
                        wafer_group_timestreams_key,
                        fit_poly_order, fit_model_len,
                        bolo_props = None, wiring_map = None,
           ):
    wafer_template_group_inv = get_template_groups_inverse(
        bolo_props, wiring_map,
        per_band = False, per_wafer=True, per_squid = False)

    frb_events = copy.copy(frame[frb_events_key])

    ts_map = frame[ts_map_key]
    variances = frame[var_key]
    wafer_sum_timestreams = frame[wafer_group_timestreams_key]

    for i, ev in enumerate(frb_events):
        #calculates the detector pointing
        det_ra = core.G3MapDouble()        
        det_dec = core.G3MapDouble()
        dets = core.G3VectorString( [ di.bid for di in ev.det_info])
        #handles calculating the sum timestreams
        sum_tses = core.G3TimestreamMap()
        for j, dinf in enumerate(ev.det_info):
            bid = dinf.bid
            wafer = wafer_template_group_inv[bid]
            if not wafer in sum_tses:
                sum_tses[wafer] = wafer_sum_timestreams[wafer]
            sum_tses[wafer] -= ts_map[bid]/(variances[bid])
        #get the variances:
        sum_variances = core.G3MapDouble()
        for k, v in sum_tses.items():
            std = np.std(v)
            if std == 0:
                std = 1e-60
            sum_variances[k] = std * std

        max_sig = -1e12
        scan_index = ev.scan_index
        for di in ev.det_info:
            if di.significance > max_sig:
                scan_index = di.scan_index
                max_sig = di.significance

        #now estimates the signficance from that
        amp_map = core.G3TimestreamMap()
        heavi_amp_map = core.G3TimestreamMap()

        ll_baseline_map = core.G3TimestreamMap()
        get_ts_poly_delta_heavi_ll(sum_tses, sum_variances,
                   fit_model_len, fit_poly_order,
                   True, False, #heavi, delta,
                   ll_baseline_map, amp_map, heavi_amp_map,
                   scan_index)

        ll_model_map = core.G3TimestreamMap()
        get_ts_poly_delta_heavi_ll(sum_tses, sum_variances,
                   fit_model_len, fit_poly_order,
                   True, True, #heavi, delta,
                   ll_model_map, amp_map, heavi_amp_map,
                   scan_index)

        #handles filling in the data
        for j, dinf in enumerate(ev.det_info):
            #handles device identification
            bid = dinf.bid

            wafer = wafer_template_group_inv[bid]

            dinf.variance = variances[bid]
            dinf.wafer_sig = ll_model_map[wafer][0] - ll_baseline_map[wafer][0]

            ev.det_info[j] = dinf
        frb_events[i] = ev        
    del frame[frb_events_key]
    frame[frb_events_key] = frb_events

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap',
                           observation_name = 'SourceName', 
                           #observation_number = 'ObservationNumber'
                           observation_number = 'ObservationID'
                           )
def AddExtraFrbInfo(frame, 
                    bs_ra_key, bs_dec_key,
                    frb_events_key,
                    ts_len_key = None,
                    scan_number_key = None,
                    observation_name = None, observation_number = None,
                    bolo_props = None, wiring_map = None ):
    frb_events = copy.copy(frame[frb_events_key])
    
    bs_ra = frame[bs_ra_key]
    bs_dec = frame[bs_dec_key]
    scan_num = frame[scan_number_key]

    ts_len = frame[ts_len_key]
    for i, ev in enumerate(frb_events):
        #calculates the detector pointing
        det_ra = core.G3MapDouble()                
        det_dec = core.G3MapDouble()
        dets = core.G3VectorString( [ di.bid for di in ev.det_info])
        mapmaker.bs_pointing_to_bolo_delta_alpha_one_samp(
            bs_ra[ev.scan_index],bs_dec[ev.scan_index],
            bolo_props, dets, det_ra, det_dec)
        
        ev.observation_name = observation_name
        ev.observation_number = observation_number
        ev.scan_number = scan_num

        #now estimates the signficance from that
        #handles filling in the data
        for j, dinf in enumerate(ev.det_info):
            #handles device identification
            bid = dinf.bid
            wm = wiring_map[bid]
            dinf.board_id = wm.board_serial #just needs to be unique for board
            dinf.module_id = wm.module

            dinf.ra = det_ra[bid]
            dinf.dec = det_dec[bid]

            dinf.was_injected = False
            ev.det_info[j] = dinf
        frb_events[i] = ev

    del frame[frb_events_key]
    frame[frb_events_key] = frb_events


def AddQFrbInfo(frame, 
                frb_events_key):
    if frame.type != core.G3FrameType.Scan:
        return
    frb_events = copy.copy(frame[frb_events_key])
    for i, ev in enumerate(frb_events):
        for j, dinf in enumerate(ev.det_info):
            #dinf.bid
            #dinf.scan_index
            dinf.q_ll_model = frame['Q_ModelLogLike'][dinf.bid][dinf.scan_index]
            dinf.q_ll_baseline = frame['Q_BaselineLogLike'][dinf.bid][dinf.scan_index]
            dinf.q_significance = dinf.q_ll_model - dinf.q_ll_baseline
            dinf.q_amplitude= frame['Q_FitDeltaAmplitude'][dinf.bid][dinf.scan_index]
            dinf.q_heavi_amp= frame['Q_FitHeaviAmplitude'][dinf.bid][dinf.scan_index]
            ev.det_info[j] = dinf
        frb_events[i] = ev
    del frame[frb_events_key]
    frame[frb_events_key] = frb_events
    return frame

def AddFrbEvTimeStamps(frame, 
                       ts_map_key,
                       frb_events_key):
    if frame.type != core.G3FrameType.Scan:
        return
    tsm = frame[ts_map_key]
    times = tsm[tsm.keys()[0]].times()
    frb_events = copy.copy(frame[frb_events_key])
    for i, ev in enumerate(frb_events):
        ev.event_time = times[ev.scan_index]
        frb_events[i] = ev
    del frame[frb_events_key]
    frame[frb_events_key] = frb_events
    return frame

@core.pipesegment
def LookForFrbEvents(pipe,
                     in_ts_map_key, filter_poly_order, filter_mhpf_cutoff, 
                     fit_poly_order, model_width, find_ll_thresh, other_ll_thresh,
                     do_rolling_mean = False, rolling_mean_width = 31,
                     
                     search_width = 2, #53
                     
                     min_peak_distance = 11, 

                     min_live_on_squid = 4,

                     clear_big_data = True,
                     scan_number_key = 'ScanNumber', 
                     sample_rate_key = 'SampleRate',
                     q_ts_map_key = None
                     ): 
    
    pipe.Add(FilterBumScans, ts_key = in_ts_map_key)
    #poly filter
    if do_rolling_mean:
        pipe.Add( todfilter.RollingMeanFilter,
                  ts_in_key = in_ts_map_key,
                  ts_out_key = 'PolyFilteredTs',
                  filter_width = rolling_mean_width
        )
    else:
        pipe.Add( todfilter.MaskedPolyHpf, 
                  in_ts_map_key = in_ts_map_key,
                  out_ts_map_key = 'PolyFilteredTs',
                  poly_order = filter_poly_order,
                  high_pass_freq_cutoff = filter_mhpf_cutoff,
                  is_masked = False,
                  mask_key = '')

    #est variance
    pipe.Add(AddSampleRate, sample_rate_key = 'SampleRate',
             ts_map_key = in_ts_map_key)

    VarianceModule = todfilter.VarianceAdder
    pipe.Add( VarianceModule,
              ts_key = 'PolyFilteredTs', 
              variance_output_key= 'TsVariance')

    ts_key_to_use = 'PolyFilteredTs'

    #now that we've flagged generate the sum timestreams
    pipe.Add(todfilter.polyutils.GenerateGroupTimestreamsWeighted, 
             in_ts_map_key = ts_key_to_use, var_key = 'TsVariance',
             out_ts_map_key = 'WaferTimestreams', 
             per_band = False, per_wafer=True, per_squid = False)

    pipe.Add(todfilter.polyutils.GenerateGroupTimestreamsAbs, 
             in_ts_map_key = ts_key_to_use, 
             out_ts_map_key = 'AbsBandTimestreams', 
             per_band = True, per_wafer = False, per_squid = False)

    pipe.Add(AddSquidLiveTimes, bolo_props_key = 'BolometerProperties', 
             wiring_map_key = 'WiringMap', ts_map_key = ts_key_to_use, 
             squid_liveness_key = 'SquidLivetime',
             bolo_squid_liveness_key = 'BoloSquidLivetime')

    #pipe.Add(core.InjectDebug, type=core.G3FrameType.Scan)
    #pipe.Add(core.Dump)

    pipe.Add(core.Delete, type=core.G3FrameType.Scan, keys=["PolyFlaggedTimestreams"])

    #fill in log likes
    pipe.Add( PolyLikelihoodFiller, 
              model_len = model_width,
              poly_order = fit_poly_order,
              include_heaviside = True,
              include_delta = False,
              ts_key = ts_key_to_use,
              variance_key = 'TsVariance',
              amp_map_output_key = 'FitDeltaAmplitude_baseline',
              hamp_map_output_key = 'FitHeaviAmplitude_baseline',#shouldn't be stored
              loglike_output_key = 'BaselineLogLike')

    pipe.Add( PolyLikelihoodFiller, 
              model_len = model_width,
              poly_order = fit_poly_order,
              include_heaviside = True,
              include_delta = True,
              ts_key = ts_key_to_use,
              variance_key = 'TsVariance',
              amp_map_output_key = 'FitDeltaAmplitude',
              hamp_map_output_key = 'FitHeaviAmplitude',            
              loglike_output_key = 'ModelLogLike')

    #find events
    pipe.Add( DeltaEventHunter,
              ll_model_key = 'ModelLogLike',
              ll_base_key = 'BaselineLogLike',
              trigger_thresh = find_ll_thresh,
              other_det_thresh = other_ll_thresh, 
              min_distance = min_peak_distance,
              output_event_key = 'FrbEvents',
              fit_amp_key = 'FitDeltaAmplitude',
              fit_hamp_key = 'FitHeaviAmplitude',
              search_width = search_width
          )
    #do the Q:
    if not q_ts_map_key is None:
        pipe.Add( todfilter.MaskedPolyHpf, 
                  in_ts_map_key = q_ts_map_key,
                  out_ts_map_key = 'Q_PolyFilteredTs',
                  poly_order = filter_poly_order,
                  high_pass_freq_cutoff = filter_mhpf_cutoff,
                  is_masked = False,
                  mask_key = '')
        pipe.Add( VarianceModule,
                  ts_key = 'Q_PolyFilteredTs', 
                  variance_output_key= 'Q_TsVariance')
        #pipe.Add(core.InjectDebug)

        #fill in log likes
        pipe.Add( PolyLikelihoodFiller, 
              model_len = model_width,
              poly_order = fit_poly_order,
              include_heaviside = True,
              include_delta = False,
              ts_key = 'Q_PolyFilteredTs',
              variance_key = 'Q_TsVariance',
              amp_map_output_key = 'Q_FitDeltaAmplitude_baseline',
              hamp_map_output_key = 'Q_FitHeaviAmplitude_baseline',#shouldn't be stored
              loglike_output_key = 'Q_BaselineLogLike')
        pipe.Add( PolyLikelihoodFiller, 
                  model_len = model_width,
                  poly_order = fit_poly_order,
                  include_heaviside = True,
                  include_delta = True,
                  ts_key = 'Q_PolyFilteredTs',
                  variance_key = 'Q_TsVariance',
                  amp_map_output_key = 'Q_FitDeltaAmplitude',
                  hamp_map_output_key = 'Q_FitHeaviAmplitude',            
                  loglike_output_key = 'Q_ModelLogLike')
        pipe.Add(AddQFrbInfo, frb_events_key = 'FrbEvents')
        pipe.Add(core.Delete, type=core.G3FrameType.Scan, keys=["Q_PolyFilteredTs"])
    #fill in the extra Frb information with our massive function
        #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
    pipe.Add(AddNSampsAboveSig, thresh = 5, output_key = 'NSampsAbove5Sig')
    pipe.Add(AddNSampsAboveSig, thresh = 6, output_key = 'NSampsAbove6Sig')
    pipe.Add(AddNSampsAboveSig, thresh = 7, output_key = 'NSampsAbove7Sig')

    pipe.Add(AddEventLiveTimeStats,
             id_ts_map_key = ts_key_to_use,
             det_id_lst_key = 'DetectorIdList',
             liveness_ts_map_key = 'TimestreamFlaggedFrb',
             ts_len_key = 'TimestreamLength',
             det_number_key = 'NumLiveDetectors',
             detector_livetime_stats_key = 'EffectiveDetectorLiveSamps'
    )
    
    pipe.Add(AddExtraFrbInfo,
             bolo_props = 'BolometerProperties', 
             wiring_map = 'WiringMap',
             #bs_ra_key = 'BoresightRa', 
             #bs_dec_key = 'BoresightDec',

             bs_ra_key = 'OnlineBoresightRa', 
             bs_dec_key = 'OnlineBoresightDec',

             scan_number_key = scan_number_key,
             frb_events_key = 'FrbEvents',
             ts_len_key = 'TimestreamLength'
    )

    pipe.Add(AddFrbWaferSig,
             ts_map_key = ts_key_to_use,
             var_key = 'TsVariance',
             frb_events_key = 'FrbEvents', 
             wafer_group_timestreams_key = 'WaferTimestreams',

             fit_poly_order = fit_poly_order, 
             fit_model_len = model_width)

    pipe.Add(AddFrbSquidSig,
             ll_model_key = 'ModelLogLike',
             ll_baseline_key = 'BaselineLogLike',
             frb_events_key='FrbEvents')

    #A bunch of processing for data characterization
    pipe.Add(HistogramLoglikelihoods, 
             model_key = 'ModelLogLike', baseline_key = 'BaselineLogLike',
             bins = get_histogram_bins(n_bins = 900, 
                                       min_bin = 0, max_bin = 90), 
             output_key = 'HistOfLoglikes',
             bin_key = 'HistOfLoglikesBins')

    # 8 sigma
    mean_sigma = 0.00121552045264 #useless magic number just for setting histogram scale
    pipe.Add(HistogramTimestreams, 
             ts_key = ts_key_to_use,
             bins = get_histogram_bins(n_bins = 1000, 
                                       min_bin = -20 * mean_sigma, 
                                       max_bin = 20 * mean_sigma), 
             output_key = 'HistOfTimestreams',
             bin_key = 'HistOfTimestreamsBins')
    
    if clear_big_data:
        pipe.Add(usefulmodules.FilterFramesOfTypes,
                 frame_type = core.G3FrameType.Scan,
                 keys_to_save = ['ModelLogLike', 'BaselineLogLike', 
                                 'WaferTimestreams', 'AbsBandTimestreams'],
                 types_to_remove = [ core.G3TimestreamMap,
                                     core.G3MapVectorDouble,
                                     core.G3MapVectorInt])
def get_frb_flags():
    ignore_ts_flags = ['too_many_glitches']
    invert_ts_flags = []
    ignore_bolo_flags = ['has_time_const', 'good_angle_fit', 'good_xpol_fit', 'has_polcal']
    invert_bolo_flags = ['has_pointing']
    enforce_partner_good = True
    return ignore_ts_flags, invert_ts_flags, ignore_bolo_flags, invert_bolo_flags, enforce_partner_good

def do_frb_search_idf(idf_fn, 
                  hwm_fn,
                  tes_pos_pkl_fn,
                  out_fn,

                  xtalk_constraint_file,

                  filter_poly_order = 17,
                  fit_poly_order = 1,
                  model_width = 17,
                  find_ll_thresh = 8.0,
                  other_ll_thresh = 8.0,
                  min_peak_distance = 11,

                  inject_fake_signal = False, 
                  
                  time_scale = 0.01,
                  curve_type = 3,
                  fluence = 200, 
                  ftsb_filename = None,

                  point_source_file = None,
                  min_live_on_squid = 4,
                  glitchy_prob_cutoff = 0.99,

                  simulate_timestreams = 0,
                  
                  skip_finding = False,
                  skip_filtering = False,
                  skip_writing_filter_data = False,
                  clear_big_data = True
              ):
    from spt3g.frbutils.frbfiltering import SptpolFrbFiltering

    pipe = core.G3Pipeline()

    if not skip_finding:
        xtalk_info = read_xtalk_file(xtalk_constraint_file)

        ignore_ts_flags, invert_ts_flags, ignore_bolo_flags, invert_bolo_flags, enforce_partner_good = get_frb_flags()
        #idf_names = map(lambda on: os.path.join(idf_key, on)+'.h5', obs_names)
        pipe.Add( DirectIdfReader, filename=idf_fn,
                  hwm_path = hwm_fn,
                  preload_data = False,
                  ignore_ts_flags=ignore_ts_flags,
                  invert_ts_flags=invert_ts_flags,
                  ignore_bolo_flags=ignore_bolo_flags,
                  invert_bolo_flags=invert_bolo_flags,
                  enforce_partner_good=enforce_partner_good,
                  #number_of_scans_to_read = 5
                  )

        #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
        pipe.Add(FilterBadScan, bad_scan_key = 'ScanIsBad')

        pipe.Add(RemoveFlaggedTimestreams,
                 input_ts_key = "CalTimestreams",
                 output_ts_key = "TimestreamFlaggedFrb",
                 input_flag_key = "Flags")

        if clear_big_data:
            pipe.Add(core.Delete, type = core.G3FrameType.Scan,
                     keys = ['CalTimestreams'])
            

        if inject_fake_signal:
            lookin_ts_key = 'InjTimestreams'
            ftsb = pickle.load(open(ftsb_filename))
            pipe.Add(InjectFrbSignal, ts_key = 'TimestreamFlaggedFrb',
                     out_ts_key = lookin_ts_key,
                     time_scale = time_scale, curve_type = curve_type,
                     fluence = fluence, fts = ftsb,
                     two_deg_s_speed = True)
        else:
            lookin_ts_key = 'TimestreamFlaggedFrb'

        if not tes_pos_pkl_fn is None:
            pipe.Add(AddTesPosition, pkl_file = tes_pos_pkl_fn)

            
        if simulate_timestreams:
            do_flat = False
            simulate_ts_map = True
            n_events_to_simulate = 50

            if simulate_timestreams == 2:
                do_flat = True
                simulate_ts_map = False
                n_events_to_simulate = 1000
            elif simulate_timestreams == 3:
                do_flat = False
                simulate_ts_map = False
                n_events_to_simulate = 1000
            elif simulate_timestreams == 4:
                do_flat = True
                simulate_ts_map = True
                n_events_to_simulate =0

            pipe.Add(CreateFakeFrbTimestream,
                     ts_map_key = lookin_ts_key,
                     xtalk_info = xtalk_info,
                     n_to_inject = n_events_to_simulate,
                     do_flat = do_flat,
                     simulate_ts_map = simulate_ts_map,
                     out_key = 'TimestreamTotallyFake')
            lookin_ts_key = 'TimestreamTotallyFake'

        pipe.Add( LookForFrbEvents,
                  in_ts_map_key = lookin_ts_key,
                  filter_poly_order = filter_poly_order,
                  filter_mhpf_cutoff = -1,
                  fit_poly_order = fit_poly_order,
                  model_width = model_width,
                  find_ll_thresh = find_ll_thresh,
                  other_ll_thresh = other_ll_thresh,
                  min_peak_distance = min_peak_distance,
                  clear_big_data = clear_big_data,
                  min_live_on_squid = min_live_on_squid,
                  
                  )

        if inject_fake_signal:        
            pipe.Add(LookForInjectedEvents,
                     frb_event_key = 'FrbEvents')
            #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
    else:
        pipe.Add(core.G3Reader, filename = search_fn)

    pipe.Add(AddBidCount, ev_key = 'FrbEvents')

    if not skip_filtering:
        pipe.Add(SptpolFrbFiltering,
                 point_source_file = point_source_file,
                 glitchy_prob_cutoff = glitchy_prob_cutoff,
                 filter_to_waf = False,
             )

    pipe.Add(usefulmodules.FilterFramesOfTypes,
             frame_type = core.G3FrameType.Scan,
             keys_to_save = ['BandTimestreams', 'AbsBandTimestreams', 'WaferTimestreams'],
             types_to_remove = [ core.G3TimestreamMap,
                                 core.G3MapVectorDouble,
                                 core.G3MapVectorInt,
                                 core.G3MapVectorString,
                                 core.G3Timestream,
                                 ])

    pipe.Add(core.G3Writer, filename = out_fn)
    #pipe.Add(core.Dump)            
    pipe.Run(profile=PROFILE_FRB_CODE)

def do_frb_search(fn_lst,

                  tes_pos_pkl_fn,
                  out_fn,

                  #hwm_fn = None,
                  xtalk_constraint_file = None,

                  filter_poly_order = 17,
                  fit_poly_order = 1,
                  model_width = 17,
                  find_ll_thresh = 8.0,
                  other_ll_thresh = 8.0,
                  min_peak_distance = 11,

                  inject_fake_signal = False, 
                  
                  time_scale = 0.01,
                  curve_type = 3,
                  fluence = 200, 
                  ftsb_filename = None,

                  point_source_file = None,
                  min_live_on_squid = 4,
                  glitchy_prob_cutoff = 0.99,

                  simulate_timestreams = 0,
                  
                  skip_finding = False,
                  skip_filtering = False,
                  skip_writing_filter_data = False,
                  clear_big_data = True
              ):
    from spt3g.frbutils.frbfiltering import SptpolFrbFiltering

    pipe = core.G3Pipeline()

    if not skip_finding:

        pipe.Add(core.G3Reader, filename = fn_lst)
        pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

        pipe.Add(AddPixelIdSptpol)
        pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

        pipe.Add(AddScanNumber)
        #pipe.Add(core.Dump)
        #calibrate
        pipe.Add(std_processing.CalibrateRawTimestreams,
                 q_output = 'Q_CalTimestreams')

        def RemoveKcmb(frame):
            if frame.type != core.G3FrameType.Scan:
                return
            i_data = copy.copy(frame['CalTimestreams'])
            q_data = copy.copy(frame['Q_CalTimestreams'])
            for k in i_data.keys():
                i_data[k] /= core.G3Units.kelvin
                q_data[k] /= core.G3Units.kelvin
            del frame['CalTimestreams']
            del frame['Q_CalTimestreams']
            frame['CalTimestreams'] = i_data
            frame['Q_CalTimestreams'] = q_data
        #flag
        pipe.Add(RemoveKcmb)
        pipe.Add(FlagNaNs, ts_key = 'CalTimestreams')
        pipe.Add(FlagTimestreamsWithoutProperties, ts_key = 'CalTimestreams')
        pipe.Add(FlagBadHousekeeping, ts_key = 'CalTimestreams')
        pipe.Add(FlagToBand, band = 150 * core.G3Units.GHz)
        pipe.Add(FlagMissing, key_ts_key = 'CalTimestreams', 
                 flagged_ts_key = 'Q_CalTimestreams')


        pipe.Add(todfilter.VarianceAdder, ts_key = 'CalTimestreams',
                 variance_output_key='FlaggingVariance')
        pipe.Add(FlagBadG3MapValue,m_key = 'FlaggingVariance',
                 min_val = 3e-5)

        pipe.Add(RemoveFlaggedTimestreams,
                 input_ts_key = "CalTimestreams",
                 output_ts_key = "TimestreamFlaggedFrb",
                 input_flag_key = "Flags")
        pipe.Add(RemoveFlaggedTimestreams,
                 input_ts_key = "Q_CalTimestreams",
                 output_ts_key = "Q_TimestreamFlaggedFrb",
                 input_flag_key = "Flags")

        
        if clear_big_data:
            pipe.Add(core.Delete, type = core.G3FrameType.Scan,
                     keys = ['CalTimestreams', 'Q_CalTimestreams'])
        if inject_fake_signal:
            lookin_ts_key = 'InjTimestreams'
            ftsb = pickle.load(open(ftsb_filename))
            pipe.Add(InjectFrbSignal, ts_key = 'TimestreamFlaggedFrb',
                     out_ts_key = lookin_ts_key,
                     time_scale = time_scale, curve_type = curve_type,
                     fluence = fluence, fts = ftsb)
        else:
            lookin_ts_key = 'TimestreamFlaggedFrb'

        if not tes_pos_pkl_fn is None:
            pipe.Add(AddTesPosition, pkl_file = tes_pos_pkl_fn)

            
        if simulate_timestreams:
            do_flat = False
            simulate_ts_map = True
            n_events_to_simulate = 50

            if simulate_timestreams == 2:
                do_flat = True
                simulate_ts_map = False
                n_events_to_simulate = 1000
            elif simulate_timestreams == 3:
                do_flat = False
                simulate_ts_map = False
                n_events_to_simulate = 1000
            elif simulate_timestreams == 4:
                do_flat = True
                simulate_ts_map = True
                n_events_to_simulate =0

            xtalk_info = read_xtalk_file(xtalk_constraint_file)
            pipe.Add(CreateFakeFrbTimestream,
                     ts_map_key = lookin_ts_key,
                     xtalk_info = xtalk_info,
                     n_to_inject = n_events_to_simulate,
                     do_flat = do_flat,
                     simulate_ts_map = simulate_ts_map,
                     out_key = 'TimestreamTotallyFake')
            lookin_ts_key = 'TimestreamTotallyFake'

        pipe.Add( LookForFrbEvents,
                  in_ts_map_key = lookin_ts_key,
                  filter_poly_order = filter_poly_order,
                  filter_mhpf_cutoff = -1,
                  fit_poly_order = fit_poly_order,
                  model_width = model_width,
                  find_ll_thresh = find_ll_thresh,
                  other_ll_thresh = other_ll_thresh,
                  min_peak_distance = min_peak_distance,
                  clear_big_data = clear_big_data,
                  min_live_on_squid = min_live_on_squid,
                  
                  q_ts_map_key = 'Q_TimestreamFlaggedFrb',
                  )

        if inject_fake_signal:        
            pipe.Add(LookForInjectedEvents,
                     frb_event_key = 'FrbEvents')
            #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
    else:
        pipe.Add(core.G3Reader, filename = search_fn)

    pipe.Add(AddBidCount, ev_key = 'FrbEvents')

    if not skip_filtering:
        pipe.Add(usefulmodules.FilterFramesOfTypes,
             frame_type = core.G3FrameType.Scan,
                              types_to_remove = [ 
                                 gcp.TrackerStatus,
                                 gcp.TrackerPointing,
                                 dfmux.DfMuxHousekeepingMap,
                                 gcp.ACUStatusVector])
        pipe.Add(SptpolFrbFiltering,
                 point_source_file = point_source_file,
                 glitchy_prob_cutoff = glitchy_prob_cutoff,
                 filter_to_waf = False,
             )

    pipe.Add(usefulmodules.FilterFramesOfTypes,
             frame_type = core.G3FrameType.Scan,
             keys_to_save = ['BandTimestreams', 'AbsBandTimestreams', 'WaferTimestreams'],
             types_to_remove = [ core.G3TimestreamMap,
                                 core.G3MapVectorDouble,
                                 core.G3MapVectorInt,
                                 core.G3MapVectorString,
                                 core.G3Timestream,
                                 gcp.TrackerStatus,
                                 gcp.TrackerPointing,
                                 dfmux.DfMuxHousekeepingMap,
                                 gcp.ACUStatusVector,
                                 ])

    pipe.Add(core.G3Writer, filename = out_fn)
    #pipe.Add(core.Dump)            
    pipe.Run(profile=PROFILE_FRB_CODE)

def sort_input_file_lst(fns):
    import os.path, string
    a = []
    b = []
    
    bns = map(os.path.basename, fns)
    for i, bn in enumerate(bns):
        if bn[:bn.rfind('.')].isdigit():
            b.append( fns[i] )
        else:
            a.append( fns[i] )
    return a + sorted(b)


if __name__ == '__main__':
    import argparse, glob
    
    idf_fn = '/spt/user/nlharr/idfs/data/ra23h30dec-55_idf_20120224_234758_150ghz.h5'
    hwm_fn = '/spt/user/nlharr/frb_side_products/hwm_20140318/'
    tes_pos_pkl = '/spt/user/nlharr/frb_side_products/TES_positions_150ghz.pkl'
    point_source_file = '/spt/user/nlharr/frb_side_products/point_source_map.fits'
    ftsb_file_name = '/spt/public/nlharr/frb_side_products/fast_transient_signal_conversion.pkl'
    xtalk_file = '/home/nlharr/frb_side_products/xtalk-from-venus-2013'
    

    fn_lst = (
        ['/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/-39387016.g3']+ 
        glob.glob('/spt/user/nwhitehorn/sptpol/converted/ra0hdec-57.5/-39387016/*.g3'))
    #print(fn_lst)

    parser = argparse.ArgumentParser()

    parser.add_argument("--fn_lst", type=str, 
                        default = fn_lst, nargs = '*')
    #parser.add_argument("--idf_fn", type=str, default = idf_fn)
    #parser.add_argument("--hwm_fn", type=str, default = hwm_fn)
    parser.add_argument("--tes_pos_pkl_fn", type=str, default = tes_pos_pkl)
    parser.add_argument("--out_fn", type=str, default = 'processed.g3')

    parser.add_argument("--xtalk_fn", type=str, default = xtalk_file)
    
    parser.add_argument("--filter_poly_order", type=int, default = 11)
    parser.add_argument("--fit_poly_order", type=int, default = 1)
    parser.add_argument("--model_width", type=int, default = 17)
    parser.add_argument("--find_ll_thresh", type=float, default = 6.0)
    parser.add_argument("--other_ll_thresh", type=float, default = 6.0)
    parser.add_argument("--min_peak_distance", type=int, default = 11)
    
    parser.add_argument("--inject_fake_signal", type = int, default = 0)
    parser.add_argument("--time_scale", type = float, default = 0.01, help='Is in units of seconds since we cant have access to g3units from the command line')
    parser.add_argument("--curve_type", type = int, default = 3)
    parser.add_argument("--fluence", type = float, default = 200, help = 'Jy ms')
    parser.add_argument("--ftsb_filename", type = str, default = ftsb_file_name)
    
    #parser.add_argument("--point_source_file", type = str, default = point_source_file)
    parser.add_argument("--point_source_file", type = str, default = None)
    parser.add_argument("--min_live_on_squid", type = int, default = 4)
    parser.add_argument("--glitchy_prob_cutoff", type = float, default = 0.97)

    parser.add_argument("--simulate_timestreams", type = int, default = 0)
    
    parser.add_argument("--skip_finding", type = int, default = 0)
    parser.add_argument("--skip_filtering", type = int, default = 0)
    parser.add_argument("--skip_writing_filter_data", type = int, default = 0)
    parser.add_argument("--clear_big_data", type = int, default = 1)

    parser.add_argument("--random_seed", type = int, default = 1)
    
    args = parser.parse_args()

    import random, time
    random.seed((args.random_seed * int(time.time()*100)) % 4294967295 )
    np.random.seed((args.random_seed * int(time.time()*100)) % 4294967295 )


    fn_lst = sort_input_file_lst(args.fn_lst)
    #print(fn_lst)
    do_frb_search(
        fn_lst = fn_lst,
        #idf_fn = args.idf_fn, 
        #hwm_fn = args.hwm_fn,
                  tes_pos_pkl_fn = args.tes_pos_pkl_fn,

                  out_fn = args.out_fn,

                  xtalk_constraint_file = args.xtalk_fn,
                  filter_poly_order = args.filter_poly_order,
                  fit_poly_order = args.fit_poly_order,
                  model_width = args.model_width,
                  find_ll_thresh = args.find_ll_thresh,
                  other_ll_thresh = args.other_ll_thresh,
                  min_peak_distance = args.min_peak_distance,
                  
                  inject_fake_signal = args.inject_fake_signal,
                  time_scale = args.time_scale,
                  curve_type = args.curve_type,
                  fluence = args.fluence,
                  ftsb_filename = args.ftsb_filename,
                  
                  point_source_file = args.point_source_file,

                  min_live_on_squid = args.min_live_on_squid,
                  
                  glitchy_prob_cutoff = args.glitchy_prob_cutoff,

                  simulate_timestreams = args.simulate_timestreams,
                  
                  skip_finding = args.skip_finding,
                  skip_filtering = args.skip_filtering,
                  skip_writing_filter_data = args.skip_writing_filter_data,
                  clear_big_data = args.clear_big_data
    )
    
#madness  410592006
#madness2 397647676
