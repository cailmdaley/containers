from spt3g.frbutils import frbbackgrounds
import spt3g.core
U = spt3g.core.G3Units

from spt3g.util import genericutils as GU
from spt3g import core, todfilter, mapmaker, frbutils


import scipy.misc
import scipy.integrate
import scipy.stats
from scipy import interpolate

import numpy as np
import pylab as pl

import pickle

#sum p = 1
'''
aint = a0 / ( 1 + alpha) F_{start}^{alpha + 1}

a0: fluence PL normalization

alpha: fluence PL index

gamma:  spectrum PL index
'''



def do_pf_int( alpha, fluences, p_founds, fluence_co):
    low_int = scipy.integrate.simps(  (-alpha) * p_founds * (fluences)**(alpha - 1) / fluence_co**alpha,
                         fluences )
    high_int = p_founds[-1] * fluences[-1]**alpha/fluence_co**alpha
    return low_int + high_int

def poisson_likelihood(x, rates):
    return rates**x * np.exp( -1 * rates) / scipy.misc.factorial(x)

def get_estimated_rate( anorm, fnorm, alpha, gamma,
                        perc_found_fluence, perc_found,
                        detector_sky_time, detector_sky_area, 
                        nu_ours, nu_theirs,
                        fluences_in_integral = np.arange(0.1, 128, 0.1)):
    fluence_conversion = (nu_ours / nu_theirs)**gamma
    
    modified_fluences = fluences_in_integral * fluence_conversion
    
    pfound_vals = np.interp(modified_fluences,  perc_found_fluence, perc_found, 
                            left = 0, right = np.max(perc_found))
    low_f_cr = ( (anorm) / (fnorm**alpha) ) * alpha * fluences_in_integral**(alpha - 1)

    rate_per_sky_day = -np.trapz(low_f_cr * pfound_vals, fluences_in_integral) 
    upper_bound_fluence = np.max(fluences_in_integral)
    pf = np.max(perc_found)
    upper_integral_value = anorm * (upper_bound_fluence / fnorm)**alpha * pf
                        
    rate_per_sky_day *= upper_integral_value
    
    det_rate = rate_per_sky_day * (detector_sky_area / (4.0 * np.pi * U.radians**2.0)) * (detector_sky_time/(U.day))
    return det_rate

def get_anorm_prior(anorm_values):
    '''
    Rare Event Statistics Applied to Fast Radio Bursts
    95% confidence: 2 Jy ms 2397(998, 5754), 4 Jy ms 1065(308, 3681)
    at 1.4 GHz
    rough approximation right now


    returns: (pdf values, fluence_of_anorm, frequency_of_anorm)
    '''
    #bullshit, bad, wrong
    return scipy.stats.norm(loc = 2397, scale = 200).pdf( anorm_values)

def get_alpha_prior(alpha_values):
    '''
    0.52 to 1.0 at 90% confidence  (remember that for all my calculations I haven't been assuming negative alpha)
    rough approximation right now

    returns: (pdf value)
    '''
    return scipy.stats.norm(loc = -0.74, scale = .1).pdf( alpha_values)

def get_rate_likelihoods(xs, rates, bg_rate_dist, n_bg_bins):
    '''
    rate = r + bg_r
     p(x | r) = int Dbg  p(x | bg, r) * p( bg | r)
     this code assumes p(bg | r) = p(bg)
    returns 2d numpy of likelihood for counts (xs) and rates for a given bg_rate distribution
    '''
    counts, vs = np.histogram(bg_rate_dist, bins = n_bg_bins)
    centers = (vs[1:] + vs[:-1])/2.0
    integ = np.trapz(counts, centers)
    counts = counts / integ

    #centers is the actual bg rate, counts is the probability of that
    print("counts, centers", counts, centers)
    out_likelihoods = np.zeros( (len(rates), len(xs)) )
    for j, r in enumerate(rates):
        for i, x in enumerate(xs):
            out_likelihoods[j, i] = np.trapz( poisson_likelihood(x, r + centers) * counts, centers)
    return out_likelihoods

def construct_rate_to_likelihood_function( x_values, interp_rates, 
                                           bg_rate_dist, n_bg_bins):
    #x is fast changing
    # y 
    #x, y = np.meshgrid(x_values, interp_rates, indexing = 'xy')
    rate_likelihoods = get_rate_likelihoods(x_values, interp_rates, bg_rate_dist, n_bg_bins)
    print("interpolation constructionJ")
    return interpolate.interp2d(x_values, interp_rates, rate_likelihoods, 
                                kind='linear', bounds_error = True)

def map_out_rate_space(alpha_values, anorm_values, gamma_values,
                       fnorm, 
                       perc_found_fluence, perc_found,
                       detector_sky_time, detector_sky_area, 
                       nu_ours, nu_theirs,
                       fluences_in_integral = np.arange(0.1, 128, 0.1)):
    n_alpha = len(alpha_values)
    n_anorm = len(anorm_values)
    n_gamma = len(gamma_values)
    rate_matrix = np.zeros( (n_alpha, n_anorm, n_gamma))
    for i, alpha in enumerate(alpha_values):
        for j, anorm in enumerate(anorm_values):
            for k, gamma in enumerate(gamma_values):
                est_rate = get_estimated_rate( anorm, fnorm, alpha, gamma,
                                               perc_found_fluence, perc_found,
                                               detector_sky_time, detector_sky_area, 
                                               nu_ours, nu_theirs,
                                               fluences_in_integral)
                rate_matrix[i,j,k] = est_rate
    return rate_matrix

def construct_likelihood_space( rate_matrix, 
                                alphas, anorms, x_values, gammas, 
                                rate_likelihood_function, ignore_prior = False ):
    '''
    returns a 3d volume element for the 4d matrix since x is quantized


    returns p(x, alpha, anorm | gamma)
    return likelihood matrix, volume element
    '''
    
    likelihoods    = np.zeros( (len(alphas), len(anorms), len(x_values), len(gammas)) )
    
    alpha_prior = get_alpha_prior(alphas)
    anorm_prior = get_anorm_prior(anorms)
    
    for i, alpha in enumerate(alphas):
        print('like', i)
        for j, anorm in enumerate(anorms):
            for k, x in enumerate(x_values):
                for l, gamma in enumerate(gammas):
                    if ignore_prior:
                        likelihoods[i,j,k,l] = (rate_likelihood_function(x, rate_matrix[i,j,l])[0])
                    else:
                        likelihoods[i,j,k,l] = (rate_likelihood_function(x, rate_matrix[i,j,l])[0]
                                                * alpha_prior[i] * anorm_prior[j] )
    return likelihoods


def construct_ordering_space( likelihoods, 
                              alphas, anorms, x_values, gammas):
    ordering_param    = np.zeros( (len(alphas), len(anorms), len(x_values), len(gammas)) )
    for i, alpha in enumerate(alphas):
        print('ord', i)
        for j, anorm in enumerate(anorms):
            for k, x in enumerate(x_values):
                for l, gamma in enumerate(gammas):
                    #ordering_param[i,j,k,l] = likelihoods[i,j,k,l]/np.max(likelihoods[i,j,k,:])
                    ordering_param[i,j,k,l] = likelihoods[i,j,k,l]/np.max(likelihoods[:,:,k,:])
    return ordering_param


def create_x_belt_restricted(likelihoods, ordering_param, 
                             alphas, anorms, x_values, gammas,
                             alpha_index, anorm_index, gamma_index, 
                             confidence):
    dummy, dummy, x_values = np.meshgrid( alphas, anorms, x_values, indexing = 'ij')
    x_values = x_values.flatten()

    likes = likelihoods[alpha_index, anorm_index, :, gamma_index].flatten()
    ordering_param = ordering_param[alpha_index, anorm_index, :, gamma_index].flatten()

    #import pdb; pdb.set_trace()
    ordering_param, likes, x_values = GU.sort_n_lists( ordering_param, likes, x_values )

    ordering_param = ordering_param[::-1]
    likes = likes[::-1]
    x_values = x_values[::-1]
    sum_like = 0
    xs_in_belt = []
    for i, like in enumerate(likes):
        sum_like += like
        xs_in_belt.append( x_values[i] )
        if sum_like > confidence:
            break
    return xs_in_belt, ordering_param


def create_x_belt(likelihoods, ordering_param, 
                  alphas, anorms, x_values, gammas,
                  gamma_index, confidence):
    dummy, dummy, x_values = np.meshgrid( alphas, anorms, x_values, indexing = 'ij')
    x_values = x_values.flatten()

    likes = likelihoods[:,:,:, gamma_index].flatten()
    ordering_param = ordering_param[:,:,:, gamma_index].flatten()

    #import pdb; pdb.set_trace()
    ordering_param, likes, x_values = GU.sort_n_lists( ordering_param, likes, x_values )

    ordering_param = ordering_param[::-1]
    likes = likes[::-1]
    x_values = x_values[::-1]

    vol_element = (alphas[1] - alphas[0]) * (anorms[1] - anorms[0])
    #import pdb; pdb.set_trace()
    sum_like = 0
    xs_in_belt = []
    for i, like in enumerate(likes):
        sum_like += like * vol_element
        xs_in_belt.append( x_values[i])

        if sum_like > confidence:
            break
    return xs_in_belt



if __name__ == '__main__':
    from spt3g.frbutils.impulseinjections import *
    import pickle
    psf_fn = '/home/nlharr/frb_side_products/fast_transient_signal_conversion.pkl'
    search_fn = '/home/nlharr/tmp/combined_88_search_v4.g3'
    
    ldf_exposure_fn = '/home/nlharr/tmp/ldfs/ldf_test_8.0_7.pkl'
    
    ll = 14
    ldf_fn = '/home/nlharr/tmp/ldfs/ldf_test_%d.0_7.pkl'%ll
    
    percent_found_file = '/home/nlharr/tmp/percent_found_dic_v2.pkl'
    
    #gets our bg_rate_distribution
    exposure_counts, exposure_ldf = pickle.load(open(ldf_exposure_fn))
    ev_counts, ldf = pickle.load(open(ldf_fn))
    
    sim_bg_rates = frbbackgrounds.get_bg_rate_distribution(ldf, exposure_ldf, 1000)


    #gets the effective live time
    live_time = list(core.G3File(search_fn))[2]['EffectiveDetectorLiveSamps'] / (25e6 / 2**17. * U.Hz)
    det_sky_area = pickle.load(open(psf_fn)).sky_coverage
    
    #percent found to lists
    pfound_dic = pickle.load(open( percent_found_file ))
    fluences = []
    pfounds = []
    for k in sorted(pfound_dic.keys()):
        if k[2] == ll:
            fluences.append(k[0])
            pfounds.append(pfound_dic[k])

    #estimate the rate likelihoods for our background distribution convolved 
    xs = np.arange(0,80)
    rates = np.arange(0,1000,1.0)
    alphas = np.arange(-2.0,-0.1,0.1)
    anorms = np.arange(100,4000,100)
    gammas = np.arange(-2.0,2.0,0.1)


    #rate_likelihoods = get_rate_likelihoods(xs, rates, bg_rate_dist=sim_bg_rates, n_bg_bins = 10)
    rate_likelihood_func = construct_rate_to_likelihood_function(
        xs, rates, bg_rate_dist=sim_bg_rates, n_bg_bins = 10)
    print("done interpolation")

    #[alpha, anorm, gamma]
    rate_cube = map_out_rate_space(alpha_values = alphas, 
                                   anorm_values = anorms, 
                                   gamma_values = gammas,
                                   fnorm = 2.0, 
                                   perc_found_fluence = fluences, 
                                   perc_found = pfounds,
                                   detector_sky_time = live_time, 
                                   detector_sky_area = det_sky_area, 
                                   nu_ours = 150, nu_theirs = 1.4,
                                   fluences_in_integral = np.arange(0.1, 1024, 0.1))
    
    likelihood_hyper_cube = construct_likelihood_space( rate_cube, 
                                                        alphas, anorms, xs, gammas, 
                                                        rate_likelihood_func )
    ord_param_hyper_cube = construct_ordering_space( likelihood_hyper_cube, 
                                                     alphas, anorms, xs, gammas)
    
    belt_x = create_x_belt(likelihood_hyper_cube,
                       ord_param_hyper_cube,
                       alphas, anorms, xs, gammas,
                       gamma_index = 5, confidence = 0.95)
    print( min(belt_x), max(belt_x) )

