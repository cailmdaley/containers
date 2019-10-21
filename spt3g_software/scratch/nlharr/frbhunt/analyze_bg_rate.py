import argparse, pickle, copy
from spt3g import core, todfilter, calibration, coordinateutils, mapmaker, todfilter, frbutils, util
import random, time
import numpy as np
import matplotlib.pyplot as plt

from spt3g.frbutils.frbbackgrounds import *

PAIR_SPACING = 2.67272912956

ldf_set = 'a'

deg_fac = 100
exposure_ldf_index = 0

fluence_adjustment_fac = 0.88

sky_area = 9 * core.G3Units.arcmin * core.G3Units.arcmin

significances, pf_fluences, sat_filt_counts, pf_per_sig_dic, exposure_ldf, real_ldfs_s = pickle.load(open('bg_analysis_bundle.pkl'))
#significances, pf_fluences, sat_filt_counts, pf_per_sig_dic, exposure_ldf, real_ldfs_s = pickle.load(open('bg_analysis_bundle_5ms.pkl'))

live_time_samps = exposure_ldf.live_time
live_time = 0.005820766091346741 * core.G3Units.seconds * live_time_samps

#A_correction = 1.0/((live_time/(31536000 * core.G3Units.seconds)) * 
#                    (sky_area/(4*np.pi*core.G3Units.rad*core.G3Units.rad)))

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

if 1:
    print('sig',significances)
    print('counts',counts)
    #print(bgs)
    print('bgs',bgs)
    print('bgs/exp',bgs/exposure_fac)

if 0:
    counts = map(lambda c: c, counts)

    As = np.arange(100, 1e3, 1e1)
    phis = np.arange(-1., -2, -.2)
    lls = np.zeros( (len(As), len(phis)))
    for i, A in enumerate(As):
        print(i)
        for j, phi in enumerate(phis):
            lls[i,j] = ( -1 *
                         get_stellar_rate_likelihood_counts(
                             A, phi, pf_per_sig_dic, pf_fluences, significances[:],
                             counts, bgs, exposure_fac) )

    saturated_ll = get_saturated_ll_counts(counts)

    As = A_correction * As    

    X, Y = np.meshgrid(phis, As[::-1])

    lls = lls[::-1,:]

    dof = len(counts) - 2
    lls = (lls - np.max(lls))
    #lls = 2 * (s + lls)
    import scipy.stats
    plt.imshow(lls, interpolation='none', vmax = np.max(lls), vmin = np.max(lls)-12,
               extent = (np.max(phis), np.min(phis), np.min(As), np.max(As)),
               aspect = 'auto')
    plt.contour(X,Y,lls, levels = [-scipy.stats.chi2.ppf(0.68, dof), 
                                   -scipy.stats.chi2.ppf(0.95, dof)],
                                   colors='k')
    plt.title("Loglikelihood Difference")
    plt.ylabel("Sky rate in Events per (sky day)")
    plt.xlabel("$\lambda$ where # above fluence F can be written as $F^\lambda$")
    plt.show()

if 1:
    #1,082,340
    phi = -1.5

    #confs = [0.68, 0.95]
    #confs = [0.68, 0.95]
    confs = [0.9]
    #counts = map(lambda x: int(x*0.8), counts) #test for low ll bias correction
    if False:
        counts = map(lambda x: 0, counts) #test for low ll bias correction
        bgs = map(lambda x: 1, bgs) #test for low ll bias correction
        exposure_fac = 1.0
    for c in confs:
        def h(A):
            return rate_in_fc_band(A, phi, pf_per_sig_dic, pf_fluences, significances, 
                                   counts, bgs, exposure_fac, 
                                   confidences = [c],
                                   n_sims = 50000, F0=10)[0]
                                   #n_sims = 150000, F0=15)[0]
        A = binary_search_single_param(h, 100, 1)
        print(c, A, A*A_correction)


if 0:
    import pylab as pl
    counts = map(lambda c: c, counts)

    confs = [0.68, 0.9, 0.95]
    n_confs = len(confs)
    As = []
    for i in range(n_confs):
        As.append([])
    phis = np.arange(-1.1,-2.0,-0.2)
    #phis = [-1.5]
    #confs = [0.68, .9]
    for i, c in enumerate(confs):
        for phi in phis:
            print(phi)
            def h(A):
                print('calling ha')
                return rate_in_fc_band(A, phi, pf_per_sig_dic, pf_fluences, significances, 
                                       counts, bgs, exposure_fac, 
                                       confidences = [c],
                                       F0=15,
                                       n_sims = 500000)[0]
            A = binary_search_single_param(h, 200, 1)
            As[i].append(A)
    for i in range(len(As)):
        As[i] = A_correction * np.array(As[i])
    print(phis)
    print(As)
    cs = [ (30/256., 191/256., 216/256.), 
           (102/256., 74/256., 158/256.), 
           (206/256., 26/256., 170./256)]
    lbs = phis*0
    for i,c,a in zip(range(n_confs), cs, As[::]):
        pl.fill_between(phis, As[i-1] if i >0 else lbs, a, facecolor=c,
                        linewidth=0, alpha = 0.7)
        pl.plot(phis, As[i-1] if i > 0 else lbs, color = c)
    pl.legend(['0.68 Confidence', '0.9 Confidence', 
               '0.95 Confidence'])

    pl.xlim([min(phis), max(phis)])
    pl.xlabel(r'$\alpha$')
    pl.ylabel('Rate (Transient Events $sky^{-1} day^{-1})$ Above 15 Jy ms')
    pl.title("Confidence Intervals for the On Sky Rate of 1 ms 150 GHz Transients")
    pl.show()

