import numpy as np
import matplotlib.pyplot as plt
import pickle


from spt3g import core, todfilter, calibration, coordinateutils, mapmaker
from spt3g import todfilter, frbutils, util

pos_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/damn_ldfs_S3p0W4p0H0p5Nv0p0Xv1p0E-04P4p0E-03Ng0Gr2Co6p0.pkl'))[1][0]

#loc
if 1:
    f = open('/home/nlharr/frb_side_products/ptsrc_config_ra0hdec-57p5_both_50mJy.txt')
    for i in range(5):
        f.readline()
    pnt_ra = []
    pnt_dec = []
    for line in f:
        if line.strip()[0] == '#':
            continue
        ls = line.split()
        pnt_ra.append( float(ls[1]))
        pnt_dec.append( float(ls[2]))



    
    evs = filter( lambda ev: min(ev[1].det_info[0].significance, ev[1].det_info[1].significance)  > 9, enumerate(pos_ldf.pixel_evs))

    inds, evs = zip(*evs)
    ras = map(lambda ev: ev.det_info[0].ra / core.G3Units.arcmin, evs)
    decs = map(lambda ev: ev.det_info[0].dec / core.G3Units.arcmin, evs)

    pids = map(lambda ev: ev.det_info[0].bid, evs)
    snums = map(lambda ev: ev.scan_number, evs)
    sinds = map(lambda ev: ev.scan_index, evs)
    
    obs_num = map(lambda ev: ev.observation_number, evs)

    pos_amps = map(lambda ev: abs(ev.det_info[0].amplitude + ev.det_info[1].amplitude)/2.0, evs)
    pos_angs = map(lambda ev: np.mean(np.abs((np.arctan2(ev.det_info[0].q_amplitude, 
                                                 ev.det_info[0].amplitude),
                                              np.arctan2(ev.det_info[0].q_amplitude, 
                                                         ev.det_info[0].amplitude)))), evs)
    pos_sigs = map(lambda ev: abs(ev.det_info[0].significance + ev.det_info[1].significance)/2.0, evs)
    
    import copy
    dummy_, obs_num, snums, sinds, ras, decs, pids, pos_amps, pos_angs, pos_sigs, inds = util.genericutils.sort_n_lists(copy.copy(obs_num), obs_num, snums, sinds, ras, decs, pids, pos_amps, pos_angs, pos_sigs, inds)
    
    for on, ra, dec, pid, sind, snum, pamp, pang, psig, ind in zip(obs_num, ras, decs, pids, sinds, snums, pos_amps, pos_angs, pos_sigs, inds):
        print(on, snum, sind, ra, dec, pid, ' ::: ', 'amp', pamp, 'ang', pang, 'sig', psig, 'ind', ind)

