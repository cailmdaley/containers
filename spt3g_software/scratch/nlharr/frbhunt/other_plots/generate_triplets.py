import numpy as np
import matplotlib.pyplot as plt
import pickle, glob, os.path

from collections import Counter

from spt3g import core, todfilter, calibration, coordinateutils, mapmaker
from spt3g import todfilter, frbutils, util

#pos_ldf = pickle.load(open('/home/nlharr/tmp/ldfs/damn_time_ldfs_S3p0W4p0H0p5Nv0p0Xv1p0E-04P4p0E-03Ng0Gr2Co6p0.pkl'))[1][0]


#ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_good_counted_ldfs_*pkl')
#ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p5E-03Ng0Gr2Co6p5Fs1.pkl')
#ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p3E-03Ng0Gr2Co6p5Fs1.pkl')
#ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_NEW_GLTCH_G99_p75_ldfs_*')
ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_phinal_v2_G98_*')
o_nam = map(lambda s: os.path.basename(s)[:-4]+'_trip.txt', ldf_fns)
#loc

print ldf_fns, o_nam
for ldf_fn, o_fn in zip(ldf_fns, o_nam):
    out_str = ''
    pos_ldf = pickle.load(open(ldf_fn))[1][0]
    evs = filter( lambda ev: min(ev.det_info[0].significance, ev.det_info[1].significance)  > 6, pos_ldf.pixel_evs)

    print( ldf_fn, len(pos_ldf.pixel_evs))

    snums = map(lambda ev: ev.scan_number, evs)
    obs_num = map(lambda ev: ev.observation_number, evs)
    c = Counter(zip(obs_num, snums))

    evs = filter( lambda ev: min(ev.det_info[0].significance, ev.det_info[1].significance)  > 7, evs)
    #evs = filter( lambda ev: c[(ev.observation_number, ev.scan_number)] == 1, evs)

    ras = map(lambda ev: ev.det_info[0].ra / core.G3Units.rad, evs)
    decs = map(lambda ev: ev.det_info[0].dec / core.G3Units.rad, evs)
    pids = map(lambda ev: ev.det_info[0].bid, evs)
    snums = map(lambda ev: ev.scan_number, evs)
    sinds = map(lambda ev: ev.scan_index, evs)
    obs_num = map(lambda ev: ev.observation_number, evs)
    ev_times =     map(lambda ev: ev.event_time, evs)

    obs_num, snums, sinds, ras, decs, pids, ev_times = util.genericutils.sort_n_lists(obs_num, snums, sinds, ras, decs, pids, ev_times)
    
    of = open(o_fn, 'w')
    for on, ra, dec, pid, sind, snum, et in zip(obs_num, ras, decs, pids, sinds, snums, ev_times):
        out_str += str((sind, ra if ra > 0 else ra + 2 * np.pi, dec, et.isoformat()))+'\n'
        of.write(out_str)
    of.close()


