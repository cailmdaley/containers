from spt3g import core, xtalk
from spt3g.frbutils.impulseinjections import get_signal_response, add_fake_signal_to_timestream
from spt3g import todfilter

import numpy as np
import copy

#min sptpol pixel separation is 2.4324016938
# for sptpol to get

#monte carlo event simulation
def get_r_exp(exp_fac = -0.0025, max_sep = 163):
    min_sep = 1.25**2
    max_sep = max_sep**2
    A = 1.0 / ( np.exp( exp_fac * min_sep ) - np.exp( exp_fac * max_sep ) )
    v = A / exp_fac * np.exp( exp_fac * max_sep)
    tmp_val = (np.random.rand() + exp_fac * v) / A
    out_r = 1.0/(exp_fac) * np.log( tmp_val )
    return out_r**0.5

def get_r_flat(min_sep = 1, max_sep = 200):
    norm_fac = max_sep**2 - min_sep**2
    v = min_sep**2 / norm_fac
    return ( (np.random.rand() + v) * norm_fac ) ** 0.5

def get_theta():
    return np.random.rand() * 2 * np.pi

def get_random_amp(plaw = -2, a_min = 0.01, a_max = 0.5):
    assert(plaw < 0 and plaw != -1)
    norm_fac = a_min**(1. + plaw) - a_max**(1. + plaw)
    #print norm_fac
    return (np.random.rand() * norm_fac + a_max**(1. + plaw) )**(1.0 /(1. + plaw) )

# detector hitting  (slow and dumb)
def x_y_to_det_id( x,y, det_xs, det_ys, half_min_sep = 1.1):
    r_sq = half_min_sep**2.0
    for k in det_xs.keys():
        if (det_xs[k] - x)**2.0 + (det_ys[k] - y)**2.0 < r_sq:
            return k
    return None

def simulate_event(det_xs, det_ys, exp_fac = -0.05, amp_pow_law = -2, do_flat = False):
    amp_0 = get_random_amp(plaw = amp_pow_law)
    amp_1 = get_random_amp(plaw = amp_pow_law)
    theta = get_theta()

    if do_flat:
        r = get_r_flat()
    else:
        r = get_r_exp(exp_fac = exp_fac)

    ks = det_xs.keys()
    k_det = ks[np.random.randint(0, len(ks))]
    x = det_xs[k_det]
    y = det_ys[k_det]
    
    x += r * np.cos(theta)
    y += r * np.sin(theta)

    o_det = x_y_to_det_id(x,y, det_xs, det_ys)
    return k_det, o_det, amp_0, amp_1

def read_xtalk_file(fn):
    f = open(fn)
    #[id id xtalk]
    xtalk_info = []
    for line in f:
        ls =  line.strip().split()
        if len(ls) != 3:
            continue
        xtalk_info.append((ls[0], ls[1], float(ls[2])))
    f.close()
    return xtalk_info
    
def apply_xtalk(xtalk_info, in_ts_map):
    out_ts_map = copy.copy(in_ts_map)
    for xti in xtalk_info:
        if (not xti[0] in in_ts_map) or (not xti[1] in in_ts_map):
            continue
        if xti[0] == xti[1]:
            out_ts_map[ xti[0] ] += (xti[2] - 1) * in_ts_map[ xti[1] ]
        else:
            out_ts_map[ xti[0] ] += ( xti[2] ) * in_ts_map[ xti[1] ]
    return out_ts_map

def get_inverse_xtalk(xtalk_info):
    xtalk_mapping = core.G3MapMapDouble()
    for xti in xtalk_info:
        if not xti[0] in xtalk_mapping:
            xtalk_mapping[ xti[0] ] = core.G3MapDouble()
        xtalk_mapping[xti[0]][xti[1]] = xti[2]
    inv_xtalk = xtalk.xtalkinvert.invert_xtalk_matrix(xtalk_mapping)
    xtalk_inv_elems = []
    for k in inv_xtalk.keys():
        for l in inv_xtalk[k].keys():
            xtalk_inv_elems.append( (k, l, inv_xtalk[k][l]))
    return xtalk_inv_elems

def kelvin_to_amp( in_ts_map, amps_to_kcmb, go_backwards = False ):
    out_ts_map = copy.copy(in_ts_map)
    for k in out_ts_map.keys():
        if amps_to_kcmb[k] is None or amps_to_kcmb[k] == 0:
            out_ts_map[k] *= 0
            continue
        if go_backwards:
            out_ts_map[k] *= amps_to_kcmb[k]
        else:
            out_ts_map[k] /= amps_to_kcmb[k]
    return out_ts_map
        




@core.scan_func_cache_data(x_pos = 'TesPosX', y_pos = 'TesPosY',
                           amp_to_kcmb_map='AmpsToKcmb')
def CreateFakeFrbTimestream(
        frame, ts_map_key, xtalk_info,
        out_key,
        sim_params = ((-0.0025, -3 ),),
        curve_type = 3, 
        n_to_inject = 50, 
        time_scale = 1e-3, fts = None, 
        sptpol_amp_scaling = 6.66,
        simulate_ts_map = True,
        do_flat = False,
        amp_to_kcmb_map= None,
        x_pos = None, y_pos = None
        ):

    if simulate_ts_map:
        out_tsm = todfilter.polyutils.poly_filter_g3_timestream_map(frame[ts_map_key], 
                                                                    poly_order = 11)
        for k in out_tsm.keys():
            std = np.std(out_tsm[k])
            ts = copy.copy(out_tsm[k])
            ts[:] = np.random.normal(loc = 0, scale = std, size = len(ts))
            #import pdb; pdb.set_trace()
            #print('std', std, np.std(out_tsm[k]))
            out_tsm[k] = ts
    else:
        out_tsm = copy.copy(frame[ts_map_key])

    ts_len = len(out_tsm[out_tsm.keys()[0]])
        #injects our frbs
    for i in range(np.random.poisson(n_to_inject)):
        sp = sim_params[np.random.randint(len(sim_params))]
        k_det, o_det, amp_0, amp_1 = simulate_event(
            x_pos, y_pos, exp_fac = sp[0], amp_pow_law = sp[1], do_flat = do_flat)        

        amp_0 *= sptpol_amp_scaling
        amp_1 *= sptpol_amp_scaling

        s_index = np.random.randint(0, ts_len)    
        if (not o_det is None) and o_det in out_tsm:
            ts_1 = out_tsm[o_det]
        else:
            ts_1 = None

        if k_det in out_tsm:
            ts_0 = out_tsm[k_det]
        else:
            ts_0 = None

        if ts_0 is None and ts_1 is None:
            continue
        if ts_0 is None:
            ts_0 = ts_1
            ts_1 = None

        add_fake_signal_to_timestream( ts_0, None, time_scale, curve_type,
                                       fluence = amp_0, fts=None, inject_index = s_index)
        if not ts_1 is None:
            add_fake_signal_to_timestream( ts_1, None, time_scale, curve_type,
                                           fluence = amp_1, fts=None, inject_index = s_index)
            
    if simulate_ts_map:
        out_tsm = kelvin_to_amp(out_tsm, amp_to_kcmb_map)
        out_tsm = apply_xtalk(xtalk_info, out_tsm)
        out_tsm = kelvin_to_amp(out_tsm, amp_to_kcmb_map, go_backwards = True)

    frame[out_key] = out_tsm
    


