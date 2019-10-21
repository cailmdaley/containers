import numpy as np
import argparse, pickle, copy
from spt3g import core, todfilter, calibration, coordinateutils, mapmaker, todfilter, frbutils, util
import random, time
import numpy as np
import matplotlib.pyplot as plt

from spt3g.frbutils.frbbackgrounds import *

def get_bounds(c, nsig = 1):
    if nsig == 2:
        if c == 0:
            return (1e-9,1.29)
        elif c == 1:
            return (0.37,2.75)
        elif c == 2:
            return (0.74, 4.25 )
        else:
            return ( c - 2 * c**0.5, 
                     c + 2 * c**0.5 )        
    elif nsig == 1:
        if c == 0:
            return (1e-9,3.09)
        elif c == 1:
            return (0.05, 5.14)
        elif c == 2:
            return (0.36, 6.72 )
        elif c == 3:
            return (1.10, 5.3 )
        else:
            return ( c - c**0.5, 
                     c + c**0.5 )

def get_bounds_arr(bs):
    lbs = []
    ubs = []
    for b in bs:
        lb, ub = get_bounds(b)
        lbs.append(b-lb)
        ubs.append(ub-b)
    return np.array( [lbs, ubs] )

PAIR_SPACING = 2.67272912956

ldf_set = 'a'

deg_fac = 100
exposure_ldf_index = 0
fluence_adjustment_fac = 0.88
sky_area = 9 * core.G3Units.arcmin * core.G3Units.arcmin

significances, pf_fluences, sat_filt_counts, pf_per_sig_dic, exposure_ldf, real_ldfs_s = pickle.load(open('plotting_bundle.pkl'))

live_time_samps = exposure_ldf.live_time
live_time = 0.005820766091346741 * core.G3Units.seconds * live_time_samps
A_correction = 1.0/((live_time/(86400 * core.G3Units.seconds)) * 
                    (sky_area/(4*np.pi*core.G3Units.rad*core.G3Units.rad)))
significances = sorted(significances)
pf_fluences = [flu * fluence_adjustment_fac for flu in pf_fluences]

real_ldfs = []
for sig in significances:
    real_ldfs.append(real_ldfs_s[sig])

exposure_ldf.degrade(deg_fac)
for i in range(len(real_ldfs)):
    real_ldfs[i].degrade(deg_fac)

pair_bin = get_pair_bin(exposure_ldf)
bgs = []
counts = []
exposure_fac = exposure_ldf.hist_exposure[pair_bin] / np.sum(exposure_ldf.hist[:-1])

counts = []
for i, s in enumerate(sorted(significances)):
    counts.append(sat_filt_counts[s])
for i in range(len(real_ldfs)):
    bgs.append(real_ldfs[i].hist[pair_bin])


pos_counts = counts
pos_bgs = bgs
pos_exposure_fac = exposure_fac

significances, pf_fluences, sat_filt_counts, pf_per_sig_dic, exposure_ldf, real_ldfs_s = pickle.load(open('plotting_bundle_neg.pkl'))

live_time_samps = exposure_ldf.live_time
live_time = 0.005820766091346741 * core.G3Units.seconds * live_time_samps
A_correction = 1.0/((live_time/(86400 * core.G3Units.seconds)) * 
                    (sky_area/(4*np.pi*core.G3Units.rad*core.G3Units.rad)))
significances = sorted(significances)
pf_fluences = [flu * fluence_adjustment_fac for flu in pf_fluences]

real_ldfs = []
for sig in significances:
    real_ldfs.append(real_ldfs_s[sig])

exposure_ldf.degrade(deg_fac)
for i in range(len(real_ldfs)):
    real_ldfs[i].degrade(deg_fac)

pair_bin = get_pair_bin(exposure_ldf)
bgs = []
counts = []
exposure_fac = exposure_ldf.hist_exposure[pair_bin] / np.sum(exposure_ldf.hist[:-1])

counts = []
for i, s in enumerate(sorted(significances)):
    counts.append(sat_filt_counts[s])
for i in range(len(real_ldfs)):
    bgs.append(real_ldfs[i].hist[pair_bin])
neg_counts = counts
neg_bgs = bgs
neg_exposure_fac = exposure_fac

print('counts',pos_counts, neg_counts)
print('bgs',pos_bgs, neg_bgs)
print('exp', pos_exposure_fac, neg_exposure_fac)

significances = np.array(significances)

c1 = (0.1,0.1,1.0)
c2 = (0.1,1.0,0.1)
linewidth = 2

plt.errorbar(significances-.06, pos_bgs/pos_exposure_fac, 
             get_bounds_arr(pos_bgs)/pos_exposure_fac, fmt='b.', color=c1, linewidth=linewidth)
plt.errorbar(significances-.03, pos_counts, get_bounds_arr(pos_counts), fmt = 'c*', linewidth=linewidth)

plt.errorbar(significances+.03, neg_bgs/neg_exposure_fac, 
             get_bounds_arr(neg_bgs)/neg_exposure_fac, fmt='r.', linewidth=linewidth)
plt.errorbar(significances+.06, neg_counts, get_bounds_arr(neg_counts), fmt = 'm*', linewidth=linewidth)


plt.legend( ('Pos BG Rate', 'Pos Events', 
             'Neg BG Rate', 'Neg Events'),
            numpoints=1)
plt.xlim([6.5, 11.5])
plt.semilogy()

plt.ylabel("Number of Counts")
plt.title("Number of Events and Background Rate Estimates vs Significance Cutoff")
plt.xlabel("Significance Cutoff")

plt.ylim(ymin=1e0)
plt.show()
