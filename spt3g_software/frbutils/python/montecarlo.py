import random, copy
import cPickle as pickle
import numpy as np
import pylab as pl
from spt3g import core
from spt3g.core.genericutils import uniquify_list
from spt3g import frbutils
from spt3g.frbutils import FrbEventInfoVector
from spt3g.frbutils.frbanalysis import make_fake_scan_frame_from_frb_ev_lst, do_inputless_collective_analysis

def one_over_rish_lde():
    rval = random.random()
    #print rval
    if rval < 1e-12:
        rval == 1e-12
    return 1.0/rval
'''
'''
def generate_event(
        lde_func = one_over_rish_lde, #lde func doesn't include r normalization, so scale by r if you care. generates a random r value
        x_bounds = (-90, 90),
        y_bounds = (-90, 90)):
    while 1:
        ev_0 = (random.random() * (x_bounds[1] - x_bounds[0]) + x_bounds[0],
                random.random() * (y_bounds[1] - y_bounds[0]) + y_bounds[0])
        ev_1_r_off = lde_func()
        ev_1_theta_off = random.random() * 2 * np.pi
        ev_1 = (ev_0[0] + ev_1_r_off * np.cos(ev_1_theta_off),
                ev_0[1] + ev_1_r_off * np.sin(ev_1_theta_off))
        yield ev_0, ev_1, ev_1_r_off

def get_grid_pnt(grid_des, x,y):
    min_x = grid_des[0]
    max_x = grid_des[1]
    min_y = grid_des[2]
    max_y = grid_des[3]
    ngridpnts = grid_des[4]
    gp = ( int(round(np.nan_to_num((x - min_x) / ( max_x - min_x ) * ngridpnts))),
           int(round(np.nan_to_num((y - min_y) / ( max_y - min_y ) * ngridpnts))))
    return gp

def search_grid(ev, grid_des, grid  ):
    gp_base = get_grid_pnt(grid_des, ev[0], ev[1])
    gps = []
    for i in range(-1,2):
        for j in range(-1,2):
            gps.append(( gp_base[0]+ i, gp_base[1] + j))
    dlst = []
    for gp in gps:
        if gp in grid:
            dlst += grid[gp]
    out_dets = []
    for d in dlst:
        if (ev[0]-d[1])**2 + (ev[1]-d[2])**2 < 1:
            out_dets.append(d[0])
    return out_dets

def generate_grid_map(ns, xs, ys, ngridpnts = 160):
    min_x =  min(xs)
    max_x = max(xs)

    min_y =  min(ys)
    max_y =  max(ys)

    grid_des = (min_x, max_x, min_y, max_y, ngridpnts)
    grid = {}
    for i in range(len(xs)):
        gp = get_grid_pnt(grid_des, xs[i], ys[i])
        if not gp in grid:
            grid[gp] = []
        grid[gp].append( (ns[i], xs[i], ys[i]))
    return grid, grid_des


def get_dets_hit(ev, ns, xs, ys, radius = 1):
    #ev is one event
    out_dets = []
    for i in range(len(xs)):
        if (xs[i] - ev[0])**2.0 + (ys[i] - ev[1])**2.0  < 1:
            out_dets.append(ns[i])
    return out_dets


class MonteCarloWrapper(object):
    def __init__(self, dic):
        self.dic = dic
    def __iadd__(self, other):
        self.dic['lde'] += other.dic['lde']
        self.dic['frames'] += other.dic['frames']
        return self

def filter_pixel_pairs_monte(bad_dets, ns):
    pix_loc = {}
    for i, n in enumerate(ns):
        pix_loc[n] = i
    out_bds = copy.copy(bad_dets)
    for i in bad_dets:
        bd = ns[i] 
        if bd[-1] == 'X' or bd[-1] == 'Y':
            opix = 'X' if bd[-1] == 'Y' else 'Y'
            other_bad = pix_loc[bd[:-1]+opix]
            out_bds.append(other_bad)
    return out_bds

def generate_monte_carlo(n_good_evs, max_bad_frac):
    print('running', n_good_evs, max_bad_frac)
    from spt3g import core, frbutils
    from spt3g.frbutils import frbanalysis
    import pickle, random

    #d = zip(*pickle.load(open('/home/nlharr/tmp/frb_side_products/TES_positions_150ghz.pkl')))
    d = zip(*pickle.load(open('/homes/nlharr/frb_side_products/TES_positions_150ghz.pkl')))

    n_dets = len(d[0])
    max_bad = int(n_dets * max_bad_frac)
    n_bad = random.randint(0,max_bad)
    bad_dets = [random.randint(0,n_dets - 1) for i in range(n_bad)]
    print('before', len(bad_dets))
    bad_dets = filter_pixel_pairs_monte(bad_dets, list(d[0]))
    print('after', len(bad_dets))
    good_inds = filter(lambda i: not i in bad_dets, range(n_dets))

    #print('good inds', good_inds)
    ns = map(lambda i: d[0][i], good_inds)
    xs = map(lambda i: float(d[1][i]), good_inds)
    ys = map(lambda i: float(d[2][i]), good_inds)

    grid, grid_des = generate_grid_map(ns, xs, ys, ngridpnts = 40)
    good_evs = frbutils.G3VectorFrbEventInfo()
    rs = []
    print('generating ev')
    i = 0
    for ev in generate_event():
        rs.append(ev[2])
        det_lg = uniquify_list( search_grid(ev[0], grid_des, grid) + search_grid(ev[1], grid_des, grid))
        if len(det_lg) == 2:
            #print det_lg
            evinfo = frbutils.FrbEventInfo()
            evinfo.triggers = det_lg
            evinfo.trigger_ll_model = [1,1]
            evinfo.trigger_ll_baseline = [.1,.1]
            good_evs.append(evinfo)
        i += 1
        if i >= n_good_evs:
            break
    #for ev in good_evs:
    #    print len(ev.triggers)
    print('bundling evs', len(good_evs))
    frame = make_fake_scan_frame_from_frb_ev_lst(good_evs, 'MonteCarlo', 
                                                 hwm_folder = '/homes/nlharr/frb_side_products/hwm_20140318/',
                                                 pkl_file = '/homes/nlharr/frb_side_products/TES_positions_150ghz.pkl')
    frame['DetectorIdList'] = core.G3VectorString(ns)

    sf = 20
    histogram_bins = map(lambda x: x/float(sf), range(165*sf))
    #lde = frbanalysis.LatDist(histogram_bins, 'MonteCarlo', 'DetectorIdList')
    print('LDFing')
    lde = frbanalysis.CosmicRayLdf(histogram_bins, 'MonteCarlo')
    lde(frame)
    
    dumped_info = {}
    dumped_info['lde'] = lde
    dumped_info['frames'] = [frame]
    print('returning')
    return MonteCarloWrapper(dumped_info)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("outputfile")
    args = parser.parse_args()
    do_inputless_collective_analysis(20, args.outputfile, generate_monte_carlo, (1e6, 0.0))
