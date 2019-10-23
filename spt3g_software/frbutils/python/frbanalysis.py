from __future__ import print_function

from spt3g import core, mapmaker, calibration, std_processing
from spt3g.frbutils import get_num_dets_dist_away, get_hist_bin, get_num_det_pairs
from spt3g.frbutils import get_exposure_hist, get_exposure_hist_even_bin
from spt3g.util.genericutils import uniquify_list
import spt3g.util.genericutils as GU
from spt3g.calibration import template_groups

from spt3g.frbutils.frbmontecarlo import apply_xtalk, kelvin_to_amp
from spt3g.frbutils.frbhunting import AddScanNumber
import numpy as np
import copy, random, math, pickle



from spt3g.sptpol.directidfreader import DirectIdfReader, MultipleIdfReader


def fill_in_timestream_map_scott(ev, xtalk_inv):
    file_name = str(ev.observation_number)
    file_name = '/spt/user/nlharr/idfs/data/ra23h30dec-55_idf_%s_150ghz.h5' % ( file_name[:-6] + '_' + file_name[-6:]  )
    print('loading file', file_name)
    idf_reader = DirectIdfReader( filename = file_name, hwm_path = None, preload_data = True)

    the_frame = None
    cal_frame = None
    while (1):
        print("grabbing scan")
        frames = idf_reader(None)
        if frames is None or frames == []:
            break
        for frame in frames:
            if frame.type == core.G3FrameType.Calibration:
                print(frame)
                cal_frame = frame
            if frame.type != core.G3FrameType.Scan:
                continue
            if frame['ScanNumber'] == ev.scan_number:
                the_frame = frame
                break
    if cal_frame is None:
        raise RuntimeError("Shit be wonky")

    ts_map = the_frame['CalTimestreams']
    amp_to_kcmb_map = cal_frame['AmpsToKcmb']

    ts_map = kelvin_to_amp(ts_map, amp_to_kcmb_map)
    ts_map = apply_xtalk(xtalk_inv, ts_map)
    ts_map = kelvin_to_amp(ts_map, amp_to_kcmb_map, go_backwards = True)

    info_str = '%d ' % ev.scan_len
    print('storing scan')

    o_ts_map = core.G3TimestreamMap()
    for bid in ev.active_bid_variance.keys():
        o_ts_map[bid] = ts_map[bid]
    ev.event_timestreams = o_ts_map
    return ev

def fit_exponential(ev, fit_width = 20):
    import scipy.optimize
    def fit_func(x, a, b, c, d):
        return np.abs(a) * np.exp( -np.abs(b) * x) + c + d * x
    ind = ev.scan_index
    tsm = ev.event_timestreams
    tses = [tsm[ev.det_info[0].bid], tsm[ev.det_info[1].bid]]

    start_ind = ind
    stop_ind = min(ind + fit_width, len(tses[0])-1)
    ps = []
    for ts in tses:
        x = range(stop_ind-start_ind)
        y = ts[start_ind : stop_ind]
        try:
            params = scipy.optimize.curve_fit(fit_func, x, y)[0]
            ps.append(params)
        except:
            ps.append( ( 0, 0, 0, 0) )
    return ps


def get_exp_coeffs(evs):
    ros = []
    i = 0
    for ev in evs:
        print(i)
        i += 1
        params = fit_exponential(ev)
        ros.append( ( params[0][1], params[1][1] ) )
    return ros

def fit_gaussian(tsm, ev, fit_width = 20):
    import scipy.optimize

    ind = ev.scan_index
    tses = [tsm[ev.det_info[0].bid], tsm[ev.det_info[1].bid]]

    start_ind = max(0, ind - fit_width) 
    stop_ind = min(ind + fit_width, len(tses[0])-1)
    ps = []
    for ts in tses:
        x = range(stop_ind-start_ind)
        y = ts[start_ind : stop_ind]
        def fit_func(x, a, b, c, d):
            return a * np.exp( -(x - (ind-start_ind) )**2 / (2 * b**2)) + c + d * x
        try:
            params = scipy.optimize.curve_fit(fit_func, x, y)[0]
            ps.append(params)
        except:
            ps.append( ( 0, 0, 0, 0, 0) )
    return ps

def find_zero_crossing(tsm, ev, fit_width = 10):
    import scipy.signal
    ind = ev.scan_index
    tses = [tsm[ev.det_info[0].bid], tsm[ev.det_info[1].bid]]
    start_ind = max(0, ind - fit_width)
    stop_ind = min(ind + fit_width, len(tses[0])-1)
    low_crossings = []
    high_crossings = []
    for ts in tses:
        local_ind = ind - start_ind
        local_ts = scipy.signal.detrend(ts[start_ind:stop_ind])
        local_ts *= np.sign(local_ts[ local_ind])
        for i in range(fit_width):
            if local_ts[i+local_ind] < 0:
                high_crossings.append(i)
        for i in range(fit_width):
            if local_ts[i-local_ind] < 0:
                low_crossings.append(i)
    return low_crossings, high_crossings

def grab_frb_evs_scott_idf(ev, outfile):
    import matplotlib.pyplot as plt
    import scipy.signal
    file_name = str(ev.observation_number)
    file_name = '/spt/user/nlharr/idfs/data/ra23h30dec-55_idf_%s_150ghz.h5' % ( file_name[:-6] + '_' + file_name[-6:]  )
    print('loading file', file_name)
    idf_reader = DirectIdfReader( filename = file_name, hwm_path = None, preload_data = True)

    out_frames = []
    while (1):
        print("grabbing scan")
        frames = idf_reader(None)
        if frames is None or frames == []:
            break
        for frame in frames:
            if frame.type != core.G3FrameType.Scan:
                out_frames.append(frame)
                continue
            if frame['ScanNumber'] == ev.scan_number:
                out_frames.append(frame)
                break
    writer = core.G3Writer(outfile)
    for frame in out_frames:
        writer(frame)
    writer(core.G3Frame(core.G3FrameType.EndProcessing))


def grab_frb_evs_scott_g3(ev, outfile, pic_file):
    import glob
    import matplotlib.pyplot as plt
    import scipy.signal

    obs_num = ev.observation_number

    cal = '/spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/%d.g3'%obs_num
    data = sorted(glob.glob('/spt/user/nwhitehorn/sptpol/converted/ra0hdec-57.5/%d/*.g3' % obs_num))

    filenames = [cal] + data

    g3_reader = core.G3Reader(filename = filenames)

    out_frames = []
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename = filenames)
    pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
    pipe.Add(AddScanNumber)
    pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan or fr['ScanNumber'] == ev.scan_number)
    pipe.Add(std_processing.CalibrateRawTimestreams, q_output = 'Q_CalTimestreams')
    pipe.Add(plot_frb_event_mod, ev = ev, output_file = pic_file)
    pipe.Add(core.Dump)
    pipe.Add(core.G3Writer, filename = outfile)
    pipe.Run()


def plot_frb_event_scott(ev, 
                         output_file = None,
                         plot_width = 20,
                         by_wafer = False
                         ):
    import matplotlib.pyplot as plt
    import scipy.signal
    #ev.scan_index
    #ev.scan_number
    file_name = str(ev.observation_number)
    file_name = '/spt/user/nlharr/idfs/data/ra23h30dec-55_idf_%s_150ghz.h5' % ( file_name[:-6] + '_' + file_name[-6:]  )
    print('loading file', file_name)
    idf_reader = DirectIdfReader( filename = file_name, hwm_path = None, preload_data = True)

    the_frame = None

    bprops = None
    wmap = None

    while (1):
        print("grabbing scan")
        frames = idf_reader(None)
        if frames is None or frames == []:
            break
        for frame in frames:
            if frame.type != core.G3FrameType.Scan:
                if frame.type == core.G3FrameType.Wiring:
                    wmap = frame['WiringMap']
                if frame.type == core.G3FrameType.Calibration:
                    bprops = frame['BolometerProperties']
                continue
            if frame['ScanNumber'] == ev.scan_number:
                the_frame = frame
                break
    ts_map = the_frame['CalTimestreams']
    plt.clf()


    if by_wafer:
        tgroups = template_groups.get_template_groups(bprops, wmap, per_wafer = True)
        for tg in tgroups:
            if ev.det_info[0].bid in tg:
                ts_lst = tg
                break
    else:
        ts_lst = ts_map.keys()


    info_str = ''
    print("plotting others")
    for bid in ts_lst:
        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[bid])-1)
        plt.plot(scipy.signal.detrend(ts_map[bid][start_ind:stop_ind]))

    print("plotting actual")
    for i, di in enumerate(ev.det_info):
        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[di.bid])-1)
        plt.plot(scipy.signal.detrend(ts_map[di.bid][start_ind:stop_ind]), linewidth = 4)
        info_str += '%d' % (i)

    plt.title(info_str)
    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    idf_reader.idf_file_.close()



def plot_frb_event_cached(ev, 
                          cached_file,
                          output_file = None,
                          plot_width = 20,
                          by_wafer = False
                          ):
    import matplotlib.pyplot as plt
    import scipy.signal
    #ev.scan_index
    #ev.scan_number

    g3_reader = core.G3Reader(filename = cached_file)

    the_frame = None

    bprops = None
    wmap = None

    while (1):
        #print("grabbing scan")
        frames = g3_reader(None)
        if frames is None or frames == []:
            break
        scan_found = False
        for frame in frames:
            if frame.type != core.G3FrameType.Scan:
                if frame.type == core.G3FrameType.Wiring:
                    wmap = frame['WiringMap']
                if frame.type == core.G3FrameType.Calibration:
                    bprops = frame['BolometerProperties']
                continue

            if frame.type == core.G3FrameType.Scan:
                the_frame = frame
                scan_found = True
                break
        if scan_found:
            break
    ts_map = the_frame['CalTimestreams']
    plt.clf()

    if by_wafer:
        tgroups = template_groups.get_template_groups(bprops, wmap, per_wafer = True)
        for tg in tgroups:
            if ev.det_info[0].bid in tg:
                ts_lst = tg
                break
    else:
        ts_lst = ts_map.keys()


    info_str = ''
    print("plotting others")
    for bid in ts_lst:
        if bid not in ts_map:
            continue
        if np.any(np.logical_not(np.isfinite(ts_map[bid]))):
            continue


        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[bid])-1)
        plt.plot(scipy.signal.detrend(ts_map[bid][start_ind:stop_ind]))

    print("plotting actual")

    min_vals = []
    max_vals = []
    for i, di in enumerate(ev.det_info):
        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[di.bid])-1)
        ts = scipy.signal.detrend(ts_map[di.bid][start_ind:stop_ind])
        plt.plot(ts, linewidth = 4)
        info_str += '%d' % (i)
        min_vals.append(min(ts))
        max_vals.append(max(ts))

        
    plt.title(info_str)
    plt.ylim( 2 * min(min_vals), 2*max(max_vals) )


    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)



def get_frb_width_cached(ev, cached_file, fit_width = 20 ):
    import matplotlib.pyplot as plt
    import scipy.signal

    g3_reader = core.G3Reader(filename = cached_file)

    the_frame = None

    bprops = None
    wmap = None

    while (1):
        #print("grabbing scan")
        frames = g3_reader(None)
        if frames is None or frames == []:
            break
        scan_found = False
        for frame in frames:
            if frame.type == core.G3FrameType.Scan:
                the_frame = frame
                scan_found = True
                break
        if scan_found:
            break
    ts_map = the_frame['CalTimestreams']

    #gauss_fit = fit_gaussian(ts_map, ev, fit_width)
    #fit_widths = [ gauss_fit[0][2], gauss_fit[1][2] ]

    #lowx, highx = find_zero_crossing(ts_map, ev, fit_width)
    #zero_widths = [highx[0]+lowx[0], highx[1]+lowx[1]]
    #return fit_widths, zero_widths
    lowx, highx = find_zero_crossing(ts_map, ev, fit_width)
    return sum([highx[0]+lowx[0], highx[1]+lowx[1]])/2.0


@core.scan_func_cache_data(bprops = 'BolometerProperties', wmap = 'WiringMap')
def plot_frb_event_mod(the_frame, 
                       ev, 
                       output_file = None,
                       plot_width = 20,
                       by_wafer = False,
                       wmap = None,
                       bprops = None
                       ):
    import matplotlib.pyplot as plt
    import scipy.signal
    #ev.scan_index
    #ev.scan_number
    ts_map = the_frame['CalTimestreams']
    plt.clf()

    if by_wafer:
        tgroups = template_groups.get_template_groups(bprops, wmap, per_wafer = True)
        for tg in tgroups:
            if ev.det_info[0].bid in tg:
                ts_lst = tg
                break
    else:
        ts_lst = ts_map.keys()


    info_str = ''
    print("plotting others")
    for bid in ts_lst:
        if bid not in ts_map:
            continue
        if np.any(np.logical_not(np.isfinite(ts_map[bid]))):
            continue


        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[bid])-1)
        plt.plot(scipy.signal.detrend(ts_map[bid][start_ind:stop_ind]))

    print("plotting actual")

    min_vals = []
    max_vals = []
    for i, di in enumerate(ev.det_info):
        start_ind = max(0, ev.scan_index-plot_width)
        stop_ind = min( ev.scan_index+plot_width, len(ts_map[di.bid])-1)
        ts = scipy.signal.detrend(ts_map[di.bid][start_ind:stop_ind])
        plt.plot(ts, linewidth = 4)
        info_str += '%d' % (i)
        min_vals.append(min(ts))
        max_vals.append(max(ts))

        
    plt.title(info_str)
    plt.ylim( 2 * min(min_vals), 2*max(max_vals) )


    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)




def calculate_poisson_ll( expected_rates, actual_rates):
    return np.sum(actual_rates * np.log( expected_rates ) -  expected_rates)


def get_percent_found_threshes( frame, low_cutoffs, high_cutoffs, 
                                ll_0_key, ll_1_key, unity_key ):
    lls_0 = np.asarray(frame[ll_0_key])
    lls_1 = np.asarray(frame[ll_1_key])
    n_scans = frame[unity_key]
    percent_found = []
    for i in range(len(low_cutoffs)):
        #import pdb; pdb.set_trace()
        nfound = np.size(np.where( np.logical_and ( np.logical_and( lls_0 > low_cutoffs[i], lls_0 < high_cutoffs[i]),
                                                    np.logical_and( lls_1 > low_cutoffs[i], lls_1 < high_cutoffs[i]))))

        percent_found.append( float(nfound) / float(n_scans))
    return percent_found


class CosmicRayLdf(object):
    def __init__(self, 
                 pos_hist_bins, 
                 frb_events_key,
                 bolo_names_key,
                 bolo_props_key = 'BolometerProperties',
                 bolo_x_pos_key = 'TesPosX', 
                 bolo_y_pos_key = 'TesPosY',
                 
                 ts_len_key = 'TimestreamLength',

                 get_pf = False,

                 do_by_squid = False,
                 do_by_wafer = False,

                 force_not_shared_squid = False,
                 force_not_shared_board = False,

                 do_geometry = False,
                 skip_pixel_evs = True,
                 pop_shared = True,
                 cache_all_evs = False
    ):
        self.force_not_shared_squid = force_not_shared_squid
        self.force_not_shared_board = force_not_shared_board

        #key detritus
        self.frb_events_key = frb_events_key
        
        self.bolo_names_key = bolo_names_key

        self.pixel_to_bolo_map = None
        self.bolo_to_pixel_map = None
        

        self.ts_len_key = ts_len_key
        
        self.bolo_x_pos_key = bolo_x_pos_key
        self.bolo_y_pos_key = bolo_y_pos_key
        self.bolo_x_pos = None
        self.bolo_y_pos = None
        
        self.bolo_props_key = bolo_props_key
        self.pop_shared = pop_shared

        self.hist_bins = core.G3VectorDouble(pos_hist_bins) # [low, high]
        self.hist_exposure = np.zeros(len(pos_hist_bins)) # [low, high, ovflow]
        self.hist = np.zeros(len(pos_hist_bins))

        self.inverse_exp_hist = np.zeros(len(pos_hist_bins))
        
        self.n_bolos_in_evs = []

        self.do_by_squid = do_by_squid
        self.do_by_wafer = do_by_wafer

        self.wiring_map = None
        self.bias_props = None

        self.skip_pixel_evs = skip_pixel_evs
        self.pixel_evs = []


        self.num_evs_in_scan = []

        self.n_events = 0

        self.do_geometry = do_geometry

        self.live_time = 0

        self.get_pf = get_pf
        self.n_found = 0
        self.n_scans = 0


        self.cache_all_evs = cache_all_evs
        self.all_evs = []

    def degrade(self, deg_fac):
        old_hist_bins = list(self.hist_bins)
        new_hist_bins = list(self.hist_bins)
        new_hist_bins = new_hist_bins[:-1:deg_fac] + [new_hist_bins[-1]]
        
        new_hist = np.zeros(len(new_hist_bins))
        new_hist_exposure = np.zeros(len(new_hist_bins))
        for i in range(len(new_hist)-1):
            low_ind = old_hist_bins.index( new_hist_bins[i] )
            high_ind = old_hist_bins.index( new_hist_bins[i+1] )
            
            new_hist[i] = np.sum(self.hist[low_ind:high_ind])
            new_hist_exposure[i] = np.sum(self.hist_exposure[low_ind:high_ind])

        self.hist_bins = core.G3VectorDouble(new_hist_bins)
        self.hist_exposure = new_hist_exposure
        self.hist = new_hist

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wiring_map = frame['WiringMap']

        if frame.type == core.G3FrameType.Calibration:
            self.bias_props = copy.copy(frame[self.bolo_props_key])

            bolo_x_pos_tmp = frame[self.bolo_x_pos_key]
            bolo_y_pos_tmp = frame[self.bolo_y_pos_key]
            self.bolo_x_pos = core.G3MapDouble()
            self.bolo_y_pos = core.G3MapDouble()
            for k in self.bias_props.keys():
                bid = self.bias_props[k].physical_name
                if bid in bolo_x_pos_tmp:
                    self.bolo_x_pos[k] = bolo_x_pos_tmp[bid]
                    self.bolo_y_pos[k] = bolo_y_pos_tmp[bid]


        if ((frame.type == core.G3FrameType.Wiring or 
             frame.type == core.G3FrameType.Calibration) and
            (not self.bias_props is None) and 
            (not self.wiring_map is None) ):
            self.pixel_to_bolo_map = template_groups.get_template_groups(
                self.bias_props, per_band = False, per_pixel = True, include_keys = True)
            self.bolo_to_pixel_map = template_groups.get_template_groups_inverse(
                self.bias_props, per_band = False, per_pixel = True)

            self.squid_to_bolo_map = template_groups.get_template_groups(
                self.bias_props, self.wiring_map, per_band = False, per_squid = True, 
                include_keys = True)
            self.bolo_to_squid_map = template_groups.get_template_groups_inverse(
                self.bias_props, self.wiring_map, per_band = False, per_squid = True)

            self.wafer_to_bolo_map = template_groups.get_template_groups(
                self.bias_props, self.wiring_map, per_band = False, per_wafer = True, 
                include_keys = True)
            self.bolo_to_wafer_map = template_groups.get_template_groups_inverse(
                self.bias_props, self.wiring_map, per_band = False, per_wafer = True)

            self.board_to_bolo_map = template_groups.get_template_groups(
                self.bias_props, self.wiring_map, per_band = False, per_board = True, 
                include_keys = True)
            self.bolo_to_board_map = template_groups.get_template_groups_inverse(
                self.bias_props, self.wiring_map, per_band = False, per_board = True)

            #template_groups.validate_template_groups(self.bias_props, self.wiring_map)

        if frame.type != core.G3FrameType.Scan:
            return

        self.n_scans += 1
        self.live_time += len(frame[self.bolo_names_key]) * frame[self.ts_len_key]

        bid_to_pix = self.bolo_to_pixel_map
        pix_to_bid = self.pixel_to_bolo_map

        x_pos = self.bolo_x_pos
        y_pos = self.bolo_y_pos

        frb_evs = frame[self.frb_events_key]

        loop_i = 0
        n_evs = 0

        for ev in frb_evs:
            bids = copy.copy(frame[self.bolo_names_key])

            loop_i += 1

            assert(len(ev.det_info) == 2)
            lls = [abs(ev.det_info[i].significance  ) for i in range(len(ev.det_info))]
            big_ind = 0 if lls[0] > lls[1] else 1
            small_ind = 1 if lls[0] > lls[1] else 0

            if self.do_by_squid:
                squid = self.bolo_to_squid_map[ev.det_info[big_ind].bid]
                if ( squid != self.bolo_to_squid_map[ev.det_info[small_ind].bid] ):
                    continue
                new_bids = []
                for bid in self.squid_to_bolo_map[squid]:
                    if bid in bids:
                        new_bids.append(bid)
                bids = new_bids

            if self.do_by_wafer:
                wafer = self.bolo_to_wafer_map[ev.det_info[big_ind].bid]
                if ( wafer != self.bolo_to_wafer_map[ev.det_info[small_ind].bid] ):
                    continue
                new_bids = []
                for bid in self.wafer_to_bolo_map[wafer]:
                    if bid in bids:
                        new_bids.append(bid)
                bids = new_bids

            if self.force_not_shared_squid:
                assert(False)
                bids = list(bids)
                squid = self.bolo_to_squid_map[ev.det_info[big_ind].bid]

                if ( squid == self.bolo_to_squid_map[ev.det_info[small_ind].bid] ):
                    continue
                for bid in self.squid_to_bolo_map[squid]:
                    if ((bid in bids) and (bid != ev.det_info[big_ind].bid) ):
                        bids.remove(bid)
                bids = list(bids)

            if self.force_not_shared_board:
                assert(False)
                bids = list(bids)
                board = self.bolo_to_board_map[ev.det_info[big_ind].bid]
                if ( board == self.bolo_to_board_map[ev.det_info[small_ind].bid] ):
                    continue
                for bid in self.board_to_bolo_map[board]:
                    if ((bid in bids) and (bid != ev.det_info[big_ind].bid) ):
                        bids.remove(bid)
                bids = list(bids)
            if ( (not ev.det_info[0].bid in bids) or 
                 (not ev.det_info[1].bid in bids)):
                continue

            n_evs += 1

            if self.get_pf:
                print("PF")
                if (ev.det_info[0].was_injected and
                    ev.det_info[1].was_injected):
                    print("FOUND")
                    self.n_found += 1
                return

            pix = self.bolo_to_pixel_map[ev.det_info[big_ind].bid]


            assert(big_ind + small_ind == 1)
            assert(big_ind * small_ind == 0)
            
            if ( pix  == 
                 self.bolo_to_pixel_map[ev.det_info[small_ind].bid]):
                self.n_events += 1
                if not self.skip_pixel_evs:
                    self.pixel_evs.append(ev)
                if self.pop_shared:
                    continue
            
            if self.pop_shared:
                bids = set(bids)
                #pop the detectors in the pixel
                for pix_bid in self.pixel_to_bolo_map[pix]:
                    if pix_bid in bids:
                        bids.remove(pix_bid)
                bids = list(bids)

            if self.cache_all_evs:
                self.all_evs.append(ev)


            #print("histing")
            #import pdb; pdb.set_trace()
            #get the histogram bin
            hist_bin = get_hist_bin(ev.det_info[big_ind].bid,
                                    ev.det_info[small_ind].bid,
                                    x_pos,
                                    y_pos,
                                    self.hist_bins)
            assert(hist_bin >= 0)
            self.hist[hist_bin] += 1
            self.n_bolos_in_evs.append( len( bids) )

            #get the exposure histogram
            if self.do_geometry:
                exposure_hist = core.G3VectorDouble()
                get_exposure_hist(ev.det_info[big_ind].bid,
                                  bids,
                                  x_pos, y_pos,
                                  self.hist_bins,
                                  exposure_hist)
                self.hist_exposure += exposure_hist
                self.inverse_exp_hist[hist_bin] += 1.0/exposure_hist[hist_bin]

        self.num_evs_in_scan = n_evs

    def __iadd__(self, other):
        self.hist_exposure += other.hist_exposure
        self.hist += other.hist
        self.n_bolos_in_evs += other.n_bolos_in_evs
        self.n_events += other.n_events
        self.inverse_exp_hist += other.inverse_exp_hist

        self.live_time += other.live_time
        
        self.n_found += other.n_found
        self.n_scans += other.n_scans

        return self

    def __isub__(self, other):
        self.hist_exposure -= other.hist_exposure
        self.hist -= other.hist
        self.n_events -= other.n_events
        self.inverse_exp_hist -= other.inverse_exp_hist
        return self

def make_150_pairs(tes_pos_pkl):
    return map(lambda p: (p, p[:-1]+'Y'), filter(lambda s: s[-1] == 'X', zip(*(pickle.load(open(tes_pos_pkl))))[0]))

def estimate_detector_rate_from_fit(ldf, fitting_function, pair_lst, x_pos, y_pos):
    dists = []
    for p in pair_lst:
        dist = ( ( x_pos[p[0]] - x_pos[p[1]])**2 + 
                 ( y_pos[p[0]] - y_pos[p[1]])**2 )**0.5
        dists.append(dist)
    dists  = np.asarray(filter(lambda x: x, sorted(np.nan_to_num(np.array(dists)))))
    
    dist = np.mean(dist)
    print(dist)

    fit_x = (np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0
    res, dof = fit_ldf(ldf)
    #print(get_saturated_ll(ldf))

    fit_vals =  np.nan_to_num(ldf.hist / ldf.hist_exposure * np.sum(ldf.hist) )[:-1]
    fit_y = fitting_func_real(fit_x, res['x']) * np.sum(ldf.hist)
    return np.interp(dist, fit_x, fit_y)


def print_det_rate_vs_estimate(ldf_file, pair_count_file, 
                               tes_pkl_file = '/home/nlharr/frb_side_products/TES_positions_150ghz.pkl',
                               x_pos_key = 'TesPosX', 
                               y_pos_key = 'TesPosY'):
    fitting_function = fitting_func_real
    pixel_pairs = get_full_pixel_pair_list(tes_pkl_file)
    ldf = pickle.load(open(ldf_file))
    pair_count = pickle.load(open(pair_count_file))
    frame = pair_count.frame_data[0]
    x_pos = frame[x_pos_key]
    y_pos = frame[y_pos_key]

    actual_count = sum(pair_count.pair_counts.values())
    estimated = estimate_detector_rate_from_fit(ldf, fitting_function, pixel_pairs, x_pos, y_pos)
    print(actual_count, estimated)

def estimate_pair_sum_from_ldf(ldf, pair_lst, 
                               x_pos, y_pos):
    n_bolos = np.max(np.asarray(ldf.n_bolos_in_evs))
    print('n-bolos stats', np.mean(np.asarray(ldf.n_bolos_in_evs)), np.max(np.asarray(ldf.n_bolos_in_evs)), np.mean(np.asarray(ldf.n_bolos_in_evs))/np.max(np.asarray(ldf.n_bolos_in_evs)))
    print("stats", np.mean(np.asarray(ldf.n_bolos_in_evs)), np.median(np.asarray(ldf.n_bolos_in_evs)))
    print('n_bolos', n_bolos)

    ldf_det_rate_r = np.asarray(ldf.hist_bins)
    ldf_det_rate_r = (ldf_det_rate_r[1:] + ldf_det_rate_r[:-1]) / 2.0
    ldf_det_rate  = np.asarray(2 * np.nan_to_num( ldf.hist / ldf.hist_exposure * np.sum(ldf.hist) / n_bolos ))[:-1]
    inds = np.where(ldf_det_rate > 0)
    ldf_det_rate_r = ldf_det_rate_r[inds]
    ldf_det_rate = ldf_det_rate[inds]
    #print list(ldf_det_rate)
    dists = []
    for p in pair_lst:
        dist = ( ( x_pos[p[0]] - x_pos[p[1]])**2 + 
                 ( y_pos[p[0]] - y_pos[p[1]])**2 )**0.5
        dists.append(dist)
    dists  = np.asarray(filter(lambda x: x, sorted(np.nan_to_num(np.array(dists)))))
    return np.sum(np.nan_to_num(np.interp(dists, ldf_det_rate_r, ldf_det_rate)))

'''
if __name__ == '__main__':
    from glob import glob
    histogram_bins = range(165)
    filenames = glob('/cosmo/spt/nlharr/processed_data/mw13_fllt12_ollt10_injfs0_injsig0/*_frb_first_pass_filter.g3.gz')
    do_collective_analysis(filenames, 'mw13_fllt12_ollt10_injfs0_injsig0_ldf_better.pkl', CosmicRayLdf, histogram_bins, 'GroupFrbFilter')
'''



def make_pretty_actual_vs_expected(start_ts, stop_ts, expected_counts, actual_counts, min_ll = 7):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.path as path

    good_inds = np.where(np.array(stop_ts) < 10000)
    start_ts = np.array(start_ts)[good_inds]
    stop_ts = np.array(stop_ts)[good_inds]
    expected_counts = np.array(expected_counts)[good_inds]
    actual_counts = np.array(actual_counts)[good_inds]
    
    fig, ax = plt.subplots()

    # histogram our data with numpy
    data = np.random.randn(1000)
    n, bins = np.histogram(data, 50)
    
    # get the corners of the rectangles for the histogram
    left = np.array(start_ts)
    right = np.array(stop_ts)
    bottom = np.zeros(len(left))
    top = bottom + actual_counts


    # we need a (numrects x numsides x 2) numpy array for the path helper
    # function to build a compound path
    XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T
    
    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)
    
    # make a patch out of it
    patch = patches.PathPatch(
        barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)

    # update the view limits
    ax.set_xlim(left[0], right[-1])
    ax.set_ylim(bottom.min(), top.max() + 2 * top.max()**.5)

    y = expected_counts
    #ax = axs[1,0]
    er = plt.errorbar((left + right)/2.0, y, yerr=[y**.5, y**.5], linestyle = 'o', color = 'r', linewidth = 2)
    plt.legend([patch, er[2]],['Pixel Pair Counts', 'Predicted Rate ($\pm \sqrt{n}$)'] )
    plt.xlabel("Test Statistic For Most Significant TES")
    plt.ylabel("Counts")
    plt.title("Background Rate Estimation, Minimum Test Statistic of %d" % min_ll)
    #plt.ion()
    #import pdb, rlcompleter
    #pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
    #pdb.set_trace()
    plt.show()


if __name__ == '__main__':
    import pickle, glob
    import pylab as pl
    
    #fns = glob.glob('dp2/ldf_7_7_start_looking_*.pkl')
    fns = glob.glob('dp4/ldf_7_7_start_looking_*_9.000000_final.pkl')

    fns = sorted(fns, cmp = GU.str_cmp_with_numbers_sorted)
    
    start_ts = []
    stop_ts = []
    expected_counts = []
    actual_counts = []

    legend_plot = []
    legend_label = []
    for fn in fns:
        fnsplit = fn.split('_')
        high_ts = fnsplit[-3]
        low_ts = fnsplit[-4]

        #high_ts = fnsplit[-2]
        #low_ts = fnsplit[-3]
        start_ts.append(float(low_ts))
        stop_ts.append(float(high_ts))

        mfunc = lambda s: s[:s.find('.')+2] if s != '1000000000000.000000' else '$\infty$'
        low_ts_label = mfunc(low_ts)
        high_ts_label = mfunc(high_ts)
        
        pair_pkl_file = fn.replace('ldf','pair_count')
        try:
            pair_count = pickle.load(open(pair_pkl_file))
        except:
            continue
        n_pairs = sum(pair_count.pair_counts.values())
        #n_pairs = 0
        
        split_distance = 2.67272912956
        
        pl.clf()
        out_fn = fn[:fn.rfind('.')]+'.png'
        ldf = pickle.load(open(fn))
        fit_x = (np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0
        res, dof = fit_ldf(ldf)
        print('Sat fit delta', (res['fun'] - get_saturated_ll(ldf)) / dof, dof)

        '''
        import pdb, rlcompleter
        pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
        pdb.set_trace()
        '''

        #'''
        pl.annotate('Actual Counts: %d' % n_pairs, xy=(0.7, 0.9), xycoords='axes fraction', 
                    bbox=dict(facecolor='none', edgecolor='black', pad=10.0))

        #'''
        #pl.legend(['test'])

        expected_counts.append(np.interp(split_distance, fit_x, fitting_func_real(fit_x, res['x']) * np.sum(ldf.hist)))
        actual_counts.append(n_pairs)
        
        if high_ts_label == '$\infty$':
            continue
        leg, = pl.plot(fit_x, np.nan_to_num(ldf.hist / ldf.hist_exposure * np.sum(ldf.hist) )[:-1])
        legend_plot.append(leg)
        legend_label.append(' %s < Test Stat. < %s' % (low_ts_label, high_ts_label))
        pl.plot(fit_x, fitting_func_real(fit_x, res['x']) * np.sum(ldf.hist))
        pl.xlabel('Detector Pair Distance (vertical line is actual)')
        pl.ylabel('Expected Background Rate')
        pl.title('Background Rates, %s < Test Stat. < %s' % (low_ts_label, high_ts_label))
        pl.savefig(out_fn)
        #pl.show()
    pl.ylabel('Lateral Distirbution Function Normalized to Expected Count Rate')
    pl.xlabel('Detector Pair Distance (vertical line is actual)')
    pl.axvline(x = split_distance,  linestyle = '--', color = 'r')
    pl.legend(legend_plot, legend_label)
    pl.title("Lateral Distribution Function (Minimum Small Test Stat. of 7)")
    #pl.show()
    pl.clf()
    start_ts, stop_ts, expected_counts, actual_counts = GU.sort_n_lists(start_ts, stop_ts, expected_counts, actual_counts)
    print('ll of estimation', len(np.array(expected_counts)))
    ll =  calculate_poisson_ll(np.array(expected_counts), np.array(actual_counts))
    sat_ll = calculate_poisson_ll(np.array(actual_counts), np.array(actual_counts))
    print(ll, sat_ll, (ll-sat_ll)/len(expected_counts))
    
    make_pretty_actual_vs_expected(start_ts, stop_ts, expected_counts, actual_counts, min_ll = 7)
