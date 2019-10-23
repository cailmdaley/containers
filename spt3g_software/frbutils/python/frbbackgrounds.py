from __future__ import print_function
from spt3g.frbutils.frbanalysis import *
from spt3g import util 

import numpy as np
import math
import scipy.optimize
from spt3g.util import genericutils as GU

PAIR_SPACING = 2.67272912956

############################################################################
############################################################################
############################################################################

def ll_poisson(counts, rate):
    counts = int(counts)
    if rate < 1e-5:
        rate = 1e-5
    return -1 * (counts * np.log(rate) - rate)

def linear_fit_func(x, ps, fixed_bg_rate = None, pair_spacing = PAIR_SPACING):
    if fixed_bg_rate is None:
        return ps[1] + (x - pair_spacing) * ps[0]# + ps[2] * (x - pair_spacing)**2.0
    else:
        return fixed_bg_rate + (x - pair_spacing) * ps[0]

def get_model_ll(theory_rates, hist, exposure):
    #returns the inverse for minimization
    #rate is theory rate * exposure
    rate = theory_rates * exposure
    inds = np.where(rate == 0)
    rate[inds] = 1e-12
    return -1 * np.sum( (hist * np.log( rate )  -  rate))    

def get_info_from_ldfs(ldf, ldf_to_get_exposure, max_x = 1e6, min_x = 0):
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    hist = ldf.hist
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])
    hist = hist[:-1]
    exposure = exposure[:-1]

    inds = np.logical_and( x <= max_x, x >= min_x)
    x = x[inds]
    exposure = exposure[inds]
    hist = hist[inds]
    return ldf.n_events, x, hist, exposure

def get_ldf_fitting_function(x, hist, exposure, 
                            fixed_bg_rate = None):
    inds = np.where(exposure > 0)
    x = x[inds]
    hist = hist[inds]
    exposure = exposure[inds]
    def help_f(ps):
        theory_rates = linear_fit_func(x, ps, fixed_bg_rate = fixed_bg_rate)
        return get_model_ll(theory_rates, hist, exposure)
    return help_f

def fit_ldf( x, hist, exposure, counts, 
             fixed_bg_rate = None, pair_spacing = PAIR_SPACING):
    from scipy.optimize import minimize
    #fit_type = 'Powell'
    fit_type = 'nelder-mead'
    ffunc = get_ldf_fitting_function(x, hist, exposure, fixed_bg_rate)
    base_params = [0.0, np.mean( np.nan_to_num(hist/exposure)), 0.0]
    res = minimize(ffunc, base_params, method=fit_type,
                   tol = 1e-12,
                   options={'xtol': 1e-12, 'disp': False} )
    fit_params = res['x']
    fit_bg_rate = linear_fit_func(pair_spacing, fit_params)
    fit_rate = max([0, counts - fit_bg_rate])
    fit_ll = res['fun']
    #saturated_ll = get_saturated_ll(hist, exposure)
    #n_pnts = len(hist[inds])
    return fit_ll, fit_rate, fit_bg_rate, fit_params

def do_fluence_integral( phi, pf_percents, pf_fluences, F0 = 10):
    low_int = scipy.integrate.simps(  
        (-phi) * pf_percents * (pf_fluences)**(phi - 1)/F0**phi,
        pf_fluences)
    high_int = pf_percents[-1] * pf_fluences[-1]**phi/F0**phi
    #high_int = 0
    return low_int + high_int

def get_stellar_rate_likelihood_counts( A, phi, pfs, pf_fluences, 
                                        sigs, counts, bgs, exposure_fac):
    ll_sum = 0
    for i, sig in enumerate(sigs):
        pf = np.array(pfs[sig])
        #import pdb; pdb.set_trace()
        rate = A * do_fluence_integral(phi, pf, np.asarray(pf_fluences))
        #bg_rate = bgs[i]/exposure_fac
        bg_rate = profile_bg(rate,counts[i],bgs[i],exposure_fac)
        ll_sum += ll_poisson(counts[i], rate + bg_rate)
    return ll_sum

def profile_bg(rate, counts, bg_counts, exposure_fac, rate_scaling = 0.912 ):
    def ll(bg_rate):
        return ( ll_poisson(counts, rate_scaling * (rate+bg_rate)) + 
                 ll_poisson(bg_counts, bg_rate*exposure_fac) )
    from scipy.optimize import minimize
    fit_type = 'nelder-mead'
    base_params = [bg_counts/exposure_fac]
    res = minimize(ll, base_params, method=fit_type,tol = 1e-12,
                   options={'xtol': 1e-12, 'disp': False} )
    rate = res['x'][0]
    return abs(rate)

def get_feldmann_cousins_band_multi(counts, rates, bg_rates, confidences, 
                                    n_sims = 10000, rate_scaling = 0.912):
    assert(len(counts) == len(rates))
    assert(len(counts) == len(bg_rates))
    counts = tuple(counts)
    sim_rates = (np.asarray(rates) + np.asarray(bg_rates)) * rate_scaling
    sim_counts = []
    sim_ranks = []
    for i in range(n_sims):
        scounts = []
        for r in sim_rates:
            if r < 0:
                print( rates, bg_rates, rate_scaling)
            scounts.append(np.random.poisson(r))
        sim_counts.append(tuple(scounts))
        this_ll = sum([ll_poisson(s,r) for s,r in zip(scounts, sim_rates)])
        best_ll = sum([ll_poisson(s,s) for s in scounts])
        sim_ranks.append(this_ll - best_ll)
    sim_ranks, sim_counts = GU.sort_two_lists(sim_ranks, sim_counts)
    rate_in_confidence = []
    if not counts in sim_counts:
        core.log_warn("Counts Not Found in Sims, maybe you need more iterations")
    for con in confidences:
        rate_in_confidence.append(counts in sim_counts[:int(np.ceil(n_sims * con))])
    return rate_in_confidence

def rate_in_fc_band(A, phi, pfs, pf_fluences, sigs, counts, bgs, exposure_fac, 
                    confidences, n_sims = 10000, rate_scaling = 0.912, F0 = 10):
    #estimate rates
    rates = []
    bg_rates = []
    for i, sig in enumerate(sigs):
        pf = np.array(pfs[sig])
        rate = A * do_fluence_integral(phi, pf, np.asarray(pf_fluences), F0=F0)

        rates.append(rate)
        bg_rate = profile_bg(rate,counts[i],bgs[i],exposure_fac, rate_scaling = rate_scaling)
        bg_rates.append(bg_rate)
        #print(rate)
        #print(bg_rate, bgs[i]/exposure_fac, rate, counts[i])
    #print('info')
    #print('r', rates)
    #print( 'b',bg_rates)
    #print( 'c',counts)
    #print('x', np.array(bgs)/exposure_fac)
    return get_feldmann_cousins_band_multi(counts, rates, bg_rates, confidences, n_sims, 
                                           rate_scaling = rate_scaling)


def get_stellar_rate_profile_likelihood( A, phi, ldfs, pfs, pf_fluences, sigs,
                                         exposure_index = 0, bg_corr_fac = 1.0,
                                         _cached_values = {}):
    ll_sum = 0
    for sig, ldf in zip(sigs, ldfs):
        pf = np.array(pfs[sig])
        #import pdb; pdb.set_trace()
        rate = A * do_fluence_integral(phi, pf, np.asarray(pf_fluences))
        if not sig in _cached_values:
            counts, x, hist, exposure = get_info_from_ldfs(ldf, ldfs[exposure_index], 
                                                           max_x = 20, min_x = 0)
            fit_ll, fit_rate, fit_bg_rate, fit_params = fit_ldf( x, hist, exposure, counts)
            bg_rate = fit_bg_rate * bg_corr_fac
            _cached_values[sig] = (bg_rate, counts)
        else:
            bg_rate, counts = _cached_values[sig]
        #print(counts, rate, bg_rate)
        ll_sum += ll_poisson(counts, rate + bg_rate)
    return ll_sum





#Estimated Rate, Counts, Background counts.
# need belt on Sky Rate and BG Rates.

# n rates, n counts, n bg counts, n expsosure factors



def get_pair_bin(ldf_to_get_exposure, pair_spacing = PAIR_SPACING):
    ldf = ldf_to_get_exposure
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])
    dist = np.abs(x - pair_spacing) + 1e12 * (1 - np.sign(exposure[:-1]))
    min_ind = int(np.where(dist == np.min(dist))[0])
    return min_ind

def find_best_a(phi, ldfs, pfs, pf_fluences, sigs,
                exposure_index = 0, bg_corr_fac = 1.0):
    def helpf_(A):
        return  get_stellar_rate_profile_likelihood(A, phi, ldfs, pfs, pf_fluences, sigs,
                                                        exposure_index, bg_corr_fac)
    from scipy.optimize import minimize

    res = minimize(helpf_, [1e6], method='nelder-mead',
                   tol = 1e-12,
                   options={'xtol': 1e-12, 'disp': False} )
    return res['x'], res['fun']


def find_best_a_counts(phi, pfs, pf_fluences, 
                       sigs, counts, bgs, exposure_fac):
    def helpf_(A):
        return get_stellar_rate_likelihood_counts( A, phi, pfs, pf_fluences,
                                                   sigs, counts, bgs, exposure_fac)
    from scipy.optimize import minimize

    res = minimize(helpf_, [1e6], method='nelder-mead',
                   tol = 1e-12,
                   options={'xtol': 1e-12, 'disp': False} )
    return res['x'], res['fun']

def get_saturated_ll_counts(counts):
    ll_sum = 0
    for c in counts:
        ll_sum += ll_poisson(c,c)
    return ll_sum
    
############################################################################
############################################################################
############################################################################


def fitting_func_real(x, ps):
    #tval = np.abs(ps[0]) * np.exp(-1 * abs(ps[1] * x) ) + np.abs(ps[2])# + x * ps[3] 
    #tval = np.abs(ps[0]) * np.exp(-1 * abs(ps[1] * x) ) + np.abs(ps[2])# + x * ps[3] + x**2 *ps[4] + x**3 * ps[5]
    #tval = np.abs(ps[0]) * np.exp(-1 * abs(ps[1] * x) ) + np.abs(ps[2])# + x * ps[3] + x**2 *ps[4] + x**3 * ps[5]
    tval = ps[0] + x * ps[1]# + x**2 *ps[2] #+ x**3 * ps[3] + x**4 * ps[4] + x**5 * ps[5]

    #tval = np.abs(ps[0]) + ps[1] * x
    #tval = np.abs(ps[0]) + ps[1] * x
    return tval

def fitting_func_low_sig(x, ps, exp_val = 0.45):
    tval = np.abs(ps[0]) * np.exp(-1 * abs(exp_val * x) ) 
    return tval

def get_saturated_ll( hist, exposure ):
    theory_rates = np.nan_to_num(hist / exposure)
    return get_model_ll(theory_rates, hist, exposure)

def get_model_fitting_function(x, hist, exposure, counts, 
                               fixed_rate = None, 
                               skip_count_ll = False,
                               fitting_func = fitting_func_real,
                               max_x = 20.0,
                               pair_spacing = PAIR_SPACING):

    inds = np.where(np.logical_and(np.logical_and(exposure > 0, x < max_x), x > 0.0))
    #inds = np.where(np.logical_and(exposure > 0, x < max_x))
    x = x[inds]
    hist = hist[inds]
    #print(x)
    exposure = exposure[inds]
    def fit_val(ps):
        theory_rates = fitting_func(x, ps)
        ll = get_model_ll(theory_rates, hist, exposure)

        if skip_count_ll:
            return ll
        if fixed_rate is None:
            ll += ll_poisson(counts, counts)
        else:
            bg_rate = fitting_func(pair_spacing, ps)
            ll += ll_poisson( counts, fixed_rate + bg_rate)
        return ll
    return fit_val, inds


def convert_pf_to_usable_format(pf, flux_correction = 0.89):
    fluences, sigs = zip(*pf.keys())
    fluences = list(set(fluences))
    fluences_vals = map(lambda s: flux_correction * float(s.replace('p','.')), fluences)
    fluence_vals, fluences = util.genericutils.sort_two_lists(fluences_vals, fluences)

    sigs = sorted(list(set(sigs)))
    od = {}
    for sig in sigs:
        od[sig] = []
        for flu in fluences:
            od[sig].append(pf[(flu, sig)])
    return fluence_vals, od


def fit_ldf_lazy(ldf, ldf_to_get_exposure, min_x = 0, max_x = 20.0, pair_spacing = PAIR_SPACING):
    counts, x, hist, exposure =  get_info_from_ldfs(ldf, ldf_to_get_exposure, max_x, min_x = min_x)
    return fit_ldf(x, hist, exposure, counts = counts)

def closest_bin_bg_est(ldf, ldf_to_get_exposure, pair_spacing = PAIR_SPACING):
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    hist = ldf.hist
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])
    hist = hist[:-1]
    exposure = exposure[:-1]

    good_inds = np.where(exposure != 0)
    x = x[good_inds]
    hist = hist[good_inds]
    exposure = exposure[good_inds]
    
    tmp_arr = np.abs(x - pair_spacing)
    min_ind = int(np.where(tmp_arr == np.min(tmp_arr))[0])
    return (hist/exposure)[min_ind]

def plot_fit(x, hist, exposure, ps, fitting_func = fitting_func_real):
    import matplotlib.pyplot as plt
    plt.plot( x, np.nan_to_num(hist / exposure))
    plt.plot( x, fitting_func(x, ps))
    #plt.show()

def create_fake_ldf(x, exposure, ps, fitting_func = fitting_func_real):
    rate = exposure * fitting_func(x, ps)
    return np.random.poisson(rate)

def test_fake_ldf(x, hist, exposure, counts, n_plot = 10):
    import matplotlib.pyplot as plt
    fit_ll, fit_rate, fit_bg_rate, fit_params = fit_ldf( x, hist, exposure, counts)
    for i in range(n_plot):
        new_hist = create_fake_ldf(x, exposure, fit_params)
        plt.plot(new_hist)
    plt.plot(hist)
    plt.show()

def get_bg_rate_distribution(ldf, ldf_to_get_exposure,
                             n_iters, degradation_factor = 100,
                             fitting_func = fitting_func_real, 
                             max_x = None, 
                             pause_on_low_val = False, 
                             pause_low_val = 0.3):
    ldf = copy.copy(ldf)
    ldf_to_get_exposure = copy.copy(ldf_to_get_exposure)
    ldf.degrade(degradation_factor)
    ldf_to_get_exposure.degrade(degradation_factor)
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    hist = ldf.hist
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])

    hist = hist[:-1]
    exposure = exposure[:-1]

    if max_x != None:
        inds = np.where( x > max_x)
        mind = np.min(inds)
        x = x[:mind]
        hist = hist[:mind]
        exposure = exposure[:mind]

    dummy, dummy, dummy, ldf_fit_params = fit_ldf(
        x, hist, exposure, 0,  skip_count_ll = True, fitting_func = fitting_func)

    sim_bg_rates = []

    for i in range(n_iters):
        sim_hist = create_fake_ldf(x, exposure, ldf_fit_params)
        sim_best_fit_ll, sim_rate, sim_bg_rate, sim_fit_params = fit_ldf(
            x, sim_hist, exposure, 0, fitting_func=fitting_func, skip_count_ll = True)
        if pause_on_low_val and (sim_bg_rate < 0.1):
            import pdb; pdb.set_trace()
        sim_bg_rates.append(sim_bg_rate)

    return sim_bg_rates


def wilks_theorem_constraints(x, hist, exposure, counts, fixed_rates, 
                              p_cutoff = 0.95):
    '''
    Shitty wilks theorem constraints, should use a root finder to find the cutoffs, but, fuck it,
    this is just here to test the feldman cousins fit
    '''
    optimum_ll, fit_rate, fit_bg_rate, fit_params = fit_ldf(x, hist, exposure, counts)
    lls = []
    for fr in fixed_rates:
        lls.append(fit_ldf(x, hist, exposure, counts, fixed_rate = fr)[0])
    llrat_cutoff = scipy.stats.chi2.ppf(p_cutoff, 1)
    llrat = 2.0 * np.abs(optimum_ll - np.array(lls))
    inds = np.where( llrat < llrat_cutoff)

    print(llrat)
    print("fit rate", fit_rate)

    return fixed_rates[inds]

def feldmann_cousins_constraints(x, hist, exposure, counts, fixed_rates, n_iterations, 
                                 confidence_percent, fitting_func = fitting_func_real):
    '''
    ordering param is: P(x|m) / P(x|m_best)
    '''
    best_fit_ll, best_fit_rate, best_fit_bg_rate, best_fit_params = fit_ldf(
        x, hist, exposure, counts, fitting_func = fitting_func)

    dummy, dummy, dummy, ldf_fit_params = fit_ldf(
        x, hist, exposure, counts,  skip_count_ll = True, fitting_func = fitting_func)

    rates_dic = {}
    ordering_params_dic = {}

    bg_rates_dic = {}
    counts_dic = {}

    belts = []
    cf_counts = []
    for i, fr in enumerate(fixed_rates):
        ordering_params = []
        rates = []
        bg_rates = []
        sim_counts_lst = []
        for i in range(n_iterations):
            sim_hist = create_fake_ldf(x, exposure, ldf_fit_params)
            sim_counts = np.random.poisson(fr + best_fit_bg_rate)
            sim_best_fit_ll, sim_rate, sim_bg_rate, sim_fit_params = fit_ldf(
                x, sim_hist, exposure, sim_counts, fitting_func = fitting_func)
            sim_true_ll, sim_rate_fixed, sim_true_bg_rate, sim_true_params = fit_ldf(
                x, sim_hist, exposure, sim_counts, fixed_rate = fr, 
                fitting_func = fitting_func)
            '''
            print("sc", sim_counts, sim_true_bg_rate + fr)
            print("lls", ll_poisson(sim_counts, sim_counts),
                  ll_poisson(sim_counts, sim_true_bg_rate + fr))
            print('true bg, best bg, counts', sim_true_bg_rate, sim_bg_rate, sim_counts)
            print('bf', sim_best_fit_ll, sim_true_ll)
            '''
            #fit rate is 
            if(sim_best_fit_ll - sim_true_ll > 0.1):
                print('err', sim_best_fit_ll, sim_true_ll)
                continue
            ordering_param = abs(sim_best_fit_ll - sim_true_ll)
            ordering_params.append( ordering_param )
            rates.append( sim_rate )
            bg_rates.append( sim_bg_rate )
            sim_counts_lst.append(sim_counts)

        #print(ordering_params)
        ordering_params, sim_counts_lst, rates, bg_rates  = GU.sort_n_lists(ordering_params, 
                                                                            sim_counts_lst, rates, bg_rates)
        ordering_params_dic[fr] = ordering_params
        rates_dic[fr] = rates
        bg_rates_dic[fr] = bg_rates
        counts_dic[fr] = sim_counts_lst

        n_valid_samps = int( np.ceil( len(sim_counts_lst) * confidence_percent ) )
        valid_counts = sim_counts_lst[:n_valid_samps]
        belt = ( np.min(valid_counts), np.max(valid_counts) )
        belts.append(belt)
        if ( ( belt[0] <= counts ) and 
             ( counts <= belt[1]) ):
            cf_counts.append(fr)
    return best_fit_rate, cf_counts, belts, ordering_params_dic, rates_dic, bg_rates_dic, counts_dic

def get_fc_rates_and_orders(ldf, ldf_to_get_exposure, 
                            counts, fixed_rates, n_iterations, 
                            confidence_percent, degradation_factor = 100,
                            just_give_me_best_fit = False,
                            fitting_func = fitting_func_real,
                            max_x = None
                            ):
    ldf = copy.copy(ldf)
    ldf_to_get_exposure = copy.copy(ldf_to_get_exposure)
    ldf.degrade(degradation_factor)
    ldf_to_get_exposure.degrade(degradation_factor)
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    hist = ldf.hist
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])

    hist = hist[:-1]
    exposure = exposure[:-1]

    if max_x != None:
        inds = np.where( x > max_x)
        mind = np.min(inds)
        x = x[:mind]
        hist = hist[:mind]
        exposure = exposure[:mind]

    best_fit_rate, cf_counts, belts, ordering_params_dic, rates_dic, bg_rates_dic, counts_dic = feldmann_cousins_constraints(
        x, hist, exposure, counts, fixed_rates, n_iterations,
        confidence_percent = confidence_percent,  fitting_func = fitting_func  )
    if just_give_me_best_fit:
        return best_fit_rate
    else:
        return rates_dic, ordering_params_dic, bg_rates_dic, counts_dic



def get_some_constraints( ldf, ldf_to_get_exposure, counts, 
                          confidence, fixed_rates, n_iterations, 
                          wilks = True,
                          degradation_factor = 100):

    ldf = copy.copy(ldf)
    ldf_to_get_exposure = copy.copy(ldf_to_get_exposure)
    ldf.degrade(degradation_factor)
    ldf_to_get_exposure.degrade(degradation_factor)
    x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
    hist = ldf.hist
    exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])

    hist = hist[:-1]
    exposure = exposure[:-1]
    if wilks:
        return  wilks_theorem_constraints(x, hist, exposure, counts, fixed_rates,
                                                      p_cutoff = confidence)
    else:
        return feldmann_cousins_constraints(x, hist, exposure, counts, 
                                            fixed_rates, n_iterations,
                                            confidence_percent = confidence)


#used for searching fc intervals.  Because it matches the data
#assumes that f(0) = true
def binary_search_single_param(f, p_start, accuracy):
    p = p_start
    while f(p):
        if p == 0:
            p = 8
        p *= 2.0
    p_step = -p/2.0
    while 2*abs(p_step) >= accuracy:
        #print(p, p_step)
        p += p_step
        if f(p):
            p_step = abs(p_step)/2.0
        else:
            p_step = -abs(p_step)/2.0
    return p


if __name__ == '__main__':
    import pickle
    ldf_x_file = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldf_test_8.0_7.pkl'
    ldf_h_file = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldf_test_12.0_7.pkl'
    d, ldf_to_get_exposure, __ = pickle.load(open(ldf_x_file))
    counts, ldf, __ = pickle.load(open(ldf_h_file))

    #plot the fitting algorithm, for funsies.
    if False:
        print("counts", counts)
        ldf_to_get_exposure.degrade(100)
        ldf.degrade(100)
        
        x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
        hist = ldf.hist
        exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])
        
        hist = hist[:-1]
        exposure = exposure[:-1]

        fit_ll, fit_rate, fit_bg_rate, fit_params = fit_ldf( x, hist, exposure, counts)
        plot_fit(x, hist, exposure, fit_params)

    #get wilks theorem constraints
    if False:
        wvr = get_some_constraints( ldf, ldf_to_get_exposure, counts,
                                    confidence = 0.95,
                                           fixed_rates = np.arange(0,10,0.2), 
                                    n_iterations = 100, wilks = True)
        print(wvr)
        print( min(wvr),max(wvr))
    #feld cous
    if True:
        bf, fc_valid_counts, belts, ordering_params_dic, rates_dic, bg_rates_dic, counts_dic = get_some_constraints( ldf, ldf_to_get_exposure, counts,
                                                                                                                     confidence = 0.9,
                                                                                                                     fixed_rates = range(10),
                                                                                                                     n_iterations = 100, wilks = False)
        print('belts', belts)
        print('valid_counts', fc_valid_counts, counts)
        print('bf', bf)
    # fake ldf construction
    if False:
        ldf_to_get_exposure.degrade(100)
        ldf.degrade(100)
        x = ((np.asarray(ldf.hist_bins)[1:] + np.asarray(ldf.hist_bins)[:-1]) / 2.0)
        hist = ldf.hist
        exposure = ldf_to_get_exposure.hist_exposure / np.sum(ldf_to_get_exposure.hist[:-1])        
        hist = hist[:-1]
        exposure = exposure[:-1]
        test_fake_ldf(x, hist, exposure, counts)
