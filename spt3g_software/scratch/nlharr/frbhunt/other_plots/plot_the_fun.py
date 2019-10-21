import numpy as np
import matplotlib.pyplot as plt
import pickle

from collections import Counter

from spt3g import core, todfilter, calibration, coordinateutils, mapmaker
from spt3g import todfilter, frbutils, util

#pos_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/damn_time_ldfs_S3p0W4p0H0p5Nv0p0Xv1p0E-04P4p0E-03Ng0Gr2Co6p0.pkl'))[1][0]
#neg_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/damn_ldfs_S3p0W4p0H0p5Nv0p0Xv1p0E-04P4p0E-03Ng1Gr2Co6p0.pkl'))[1][0]

pos_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H1p0Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1.pkl'))[1][0]
neg_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H1p0Nv0p0Xv1p5E-04P1p1E-03Ng1Gr2Co6p5Fs1.pkl'))[1][0]



#sig hist
if 0:
    pos_sig = map(lambda ev: (ev.det_info[0].significance + ev.det_info[1].significance)/2.0, pos_ldf.pixel_evs)
    neg_sig = map(lambda ev: (ev.det_info[0].significance + ev.det_info[1].significance)/2.0, neg_ldf.pixel_evs)    
    b = plt.hist( pos_sig, bins = 1000, alpha = 0.75, color = 'r')
    plt.hist( neg_sig, bins = b[1], alpha = 0.75, color = 'b')

    plt.title("Significance Histogram, Positive Fluctuations vs Negative Fluctuations")

    plt.show()

#amp
if 0:
    pos_amp = map(lambda ev: abs(ev.det_info[0].amplitude + ev.det_info[1].amplitude)/2.0, pos_ldf.pixel_evs)
    neg_amp = map(lambda ev: abs(ev.det_info[0].amplitude + ev.det_info[1].amplitude)/2.0, neg_ldf.pixel_evs)    

    b = plt.hist( pos_amp, bins = 100, alpha = 0.75, color = 'r')
    plt.hist( neg_amp, bins = b[1], alpha = 0.75, color = 'b')
    plt.title("Amplitude Histogram, Positive Fluctuations vs Negative Fluctuations")
    plt.show()


if 0:
    pos_amp = map(lambda ev: abs(ev.det_info[0].amplitude + ev.det_info[1].amplitude)/2.0, pos_ldf.pixel_evs)
    neg_amp = map(lambda ev: abs(ev.det_info[0].amplitude + ev.det_info[1].amplitude)/2.0, neg_ldf.pixel_evs)    

    pos_sig = map(lambda ev: abs(ev.det_info[0].significance + ev.det_info[1].significance)/2.0, pos_ldf.pixel_evs)
    neg_sig = map(lambda ev: abs(ev.det_info[0].significance + ev.det_info[1].significance)/2.0, neg_ldf.pixel_evs)    

    plt.plot(pos_amp, pos_sig, 'ro')
    plt.plot(neg_amp, neg_sig, 'b*')
    plt.show()

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

    evs = filter( lambda ev: min(ev.det_info[0].significance, ev.det_info[1].significance)  > 6, pos_ldf.pixel_evs)

    snums = map(lambda ev: ev.scan_number, evs)
    obs_num = map(lambda ev: ev.observation_number, evs)
    c = Counter(zip(obs_num, snums))

    evs = filter( lambda ev: min(ev.det_info[0].significance, ev.det_info[1].significance)  > 6, evs)
    evs = filter( lambda ev: c[(ev.observation_number, ev.scan_number)] > 2, evs)

    print('DA THING', len(evs))
    ras = map(lambda ev: ev.det_info[0].ra / core.G3Units.rad, evs)
    decs = map(lambda ev: ev.det_info[0].dec / core.G3Units.rad, evs)
    pids = map(lambda ev: ev.det_info[0].bid, evs)
    snums = map(lambda ev: ev.scan_number, evs)
    sinds = map(lambda ev: ev.scan_index, evs)
    obs_num = map(lambda ev: ev.observation_number, evs)
    ev_times =     map(lambda ev: ev.event_time, evs)

    obs_num, snums, sinds, ras, decs, pids, ev_times = util.genericutils.sort_n_lists(obs_num, snums, sinds, ras, decs, pids, ev_times)
    
    colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']
    types = ['.',  'o', '+', 'x', 'h', 's', 'D']
    col_index = 0
    
    n_vals = len(colors) * len(types)

    for on, ra, dec, pid, sind, snum, et in zip(obs_num, ras, decs, pids, sinds, snums, ev_times):
        #print(on, on/24./3600., snum, sind, ra if ra > 0 else ra + 2 * np.pi, dec, pid, et.isoformat())
        print( ra if ra > 0 else ra + 2 * np.pi, dec, et.isoformat())
        t_ind = on % n_vals
        plt.plot( np.array(ra), dec, types[t_ind//len(colors)] + colors[t_ind % len(colors)], markersize=10)

    plt.title("Position of Events")
    plt.xlabel("Ra (Radian)")
    plt.ylabel("Dec (Radian)")
    plt.show()

#obs num
if 0:
    num = map(lambda ev: ev.observation_number, pos_ldf.pixel_evs)
    plt.hist(num, bins = 100)
    plt.show()

#phase
if 0:
    neg_angs = map(lambda ev: np.arctan2(-ev.det_info[0].q_amplitude, 
                                         -ev.det_info[0].amplitude),
                   neg_ldf.pixel_evs)
    pos_angs = map(lambda ev: np.arctan2(ev.det_info[0].q_amplitude, 
                                         ev.det_info[0].amplitude),
                   pos_ldf.pixel_evs)
    b = plt.hist( pos_angs, bins = 30, alpha = 0.75, color = 'r')
    plt.hist( neg_angs, bins = b[1], alpha = 0.75, color = 'b')
    plt.title("Events Phase Angle")
    plt.xlabel("I/Q Phase of Response (Rad)")
    plt.show()
