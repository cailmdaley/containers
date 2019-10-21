'''
average_cls.py

This script reads in cross spectra, such as those produce by
cross_spectra_multisubmit.py, and combines them.
'''
import os
import numpy as np
from glob import glob
from spt3g import core
from spt3g.util import files


def average_cls(data_dir,
                spectra=None,
                make_covariance=False,
                weights=None,
                verbose=False,
                rebin=None):
    """
    This function combines cross-spectra which have been produced by
    map_analysis.calculate_powerspectra() (or a wrapper thereof, such as
    cross_spectra_multisubmit.py) and saved as .pkl files. It computes the
    average, standard deviation, and optionally the full covariance matrix.

    Parameters
    ----------
    data_dir: The directory containing .pkls of cross-spectra
    spectra [None]: (list) A list of the spectra which you want to average.
        Ex. ['TT', 'TE', 'EE']. If None, all spectra found will be used.
    make_covariance [False]: (bool) Compute the noise covariance matrix for the
        specified spectra using calculate_covariance_auto(). The output will
        contain additonal 'error_XX' and 'covariance_XX' keys for each
        XX in spectra.
    weights [None]: (list or array) If provided, a weight for each set of Cls
        in the average. If not provided, we'll equally weight all Cls.

    Returns
    -------
    Dictionary keyed like the input Cl pkl files and containing the average
    power spectra, 'sigma_XX' keys storing errors on the mean, and 'n_spectra'
    recording the number of input cross-spectra.
    If make_covariance = True, will also include 'covariance_XX' keys with full
    covariance matrices, and 'error_XX" keys with the sqrt(abs(diagonals)) of
    the covariance matrices.

    Notes
    -----
    Taken from averageCls in sptpol_software authored by
        Stephen Hoover and Abby Crites.
    Barely modified by DPD for SPT-3G, April 2019
        'sigma_ps' is now the error on the mean, so no need to divide by root(N)
        Weights are now properly accounted for in this calculation.
    """
    xspectra = sorted(glob(os.path.join(data_dir, '*.pkl')))
    if len(xspectra)==0:
        raise FileNotFoundError(
            'No pkls in %s' % data_dir)

    if spectra is not None:
        if not isinstance(spectra, list):
            spectra = [spectra]

    # Check the "weights" input, and make sure it's an array.
    if weights is None:
        weights = np.ones(len(xspectra), dtype=float)
    else:
        weights = np.asarray(weights)

    final_ps = dict()
    mapspec1d = None
    for i, pth in enumerate(xspectra):
        if verbose:
            print('Loading pkl %0d out of %0d' % (i+1, len(xspectra)))
        now = files.load_pickle(pth)
        # For each spectrum, create a 2d-array with shape
        # (N_cross_spectra, N_ell_bins)
        if spectra is None:
            spectra = now.keys()
        for ps in spectra:
            if rebin is not None:
                now[ps] = now[ps].rebin(rebin)
            if mapspec1d is None:
                # Use this to return spectra as MapSpectrum1D
                mapspec1d = now[ps].copy()
            if ps not in final_ps.keys():
                final_ps[ps] = np.zeros((len(xspectra), len(now[ps])))
            final_ps[ps][i,:] = now[ps]

    for ps in list(final_ps.keys()):
        # Weighted average
        avg_cls = np.average(final_ps[ps], axis=0, weights=weights)
        avg_cls_2d = np.resize(avg_cls, np.shape(final_ps[ps]))
        delta_cls = final_ps[ps] - avg_cls_2d
        if spectra is not None and ps not in spectra:
            continue
        # Calculate the (weighted) variance of the Cls
        # The prefactor is equal to 1 / (N-1) if the weights are all 1.
        prefactor = (np.sum(weights) /
                     (np.sum(weights)**2 - np.sum(weights**2)) )
        var = prefactor*np.dot(weights, delta_cls**2)
        # Each Cl is an average, so use the variance of the
        # (weighted) mean. If all weights are 1, this is dividing by N
        var /=  (np.sum(weights)**2 / np.sum(weights**2))
        final_ps['sigma_'+ps] = np.sqrt(var)
        if make_covariance:
            if verbose:
                print('Calculating covariance matrix for %s ...' % ps)
            covariance = calculate_covariance_auto(
                final_ps[ps],
                weights=weights,
                total_cov_only=False,
                indep_cov_only=False)
            final_ps['covariance_'+ps] = covariance
            final_ps['error_'+ps] = np.array(
                [np.sqrt(np.abs(covariance['total'][ii,ii])) \
                 for ii in range(covariance['total'].shape[0])] )
        np.asarray(mapspec1d)[:] = avg_cls
        final_ps[ps] = mapspec1d.copy()

    final_ps['n_spectra'] = len(xspectra)
    return final_ps


def calculate_covariance_auto(cls,
                              maps_used=None,
                              total_cov_only=False,
                              weights=None,
                              indep_cov_only=False):
    """
    This function takes many cross-spectra for a single polarization spectrum
    and calculates the covariance due to instrument / atmospheric noise.
    This covariance does not include any sample variance contribution.
    
    This function assumes that we're crossing all maps with all other maps.
    It is an implementation of equation (A15) from Lueker et al. 2009,
    arXiv:0912.4317v1. Lueker's treatment of the covariance is based on
    Tristram et al., MNRAS 358, 833-842 (2005), but is written with analysis
    of SPT data in mind.
    
    Parameters:
    -----------
    cls: (2-d array) An array of shape (N_cross-spectra, N_ell-bins).

    maps_used [None]: (2d array or list) An (N_cross-spectra, 2) (or transpose)
        iterable describing which maps were used for each cross-spectrum.
        maps_used[0] should be a 2-element iterable where the two elements are
        indices of the maps crossed to make the cross-spectrum at index 0.
        If we don't get an input, create maps_used with
            n_maps = np.int(np.ceil(np.sqrt(2*cls.shape[0])))
            maps_used =
                [[x,y] for x in range(n_maps-1) for y in range(x+1,n_maps)]

    total_cov_only [False]: (bool) If True, then return only the total
        covariance matrix, without the separate independent and interactions
        portions. The latter portion is unreliable for noise-dominated spectra

    weights [None]: (list or array) If provided, a weight for each 
        cross-spectrum in the covariance. If not provided, we'll equally
        weight all cross-spectra.

    indep_cov_only [False]: (bool) If True, return only the independent
        covariance terms, i.e the raw variance among your spectra, without
        calculating interaction terms.  This comes in handy when your
        spectra are too noise to efficiently measure those extra terms.
    
    Returns:
    --------
    A dictionary of three 2D numpy.ndarrays, each (N_ell-bins, N_ell-bins).

    The array keyed as "total" is the estimated noise covariance matrix for 
    the set of input cross-spectra. This is the uncertainty on the mean 
    of the cross-spectra.

    The array keyed as "independent" is the portion of the noise covariance 
    which comes from the diagonal of the cross-spectrum / cross-spectrum
    covariance matrix. This is equal to the covariance of the mean of the 
    cross-spectra on the assumption that each is completely independent 
    (which isn't true).

    The array keyed as "interactions" is the off-diagonal portion of the 
    cross-spectrum / cross-spectrum covariance matrix (summed in each ell bin).
    This represents the extra covariance from the maps shared between
    cross-spectra.

    If total_cov_only is True, output is only the 2D total covariance array. 
    
    AUTHOR
    ------
    Stephen Hoover, 10 January 2014
    Barely modified by DPD for SPT-3G, April 2019
        To more closely agree with how covariance matrices are commonly
        computed, we now only use one factor each of "prefactor" and weights.
        The variance is adjusted by np.sum(weights)**2 / np.sum(weights**2).
        For uniform weights the difference should be <1%
    """
    # For convenience, pull out the number of cross-spectra here.
    n_spectra = cls.shape[0]
    
    # If we weren't given weights, weight everything equally.
    # If we are given weights, make sure they're formatted as an array.
    if weights is None:
        weights = np.ones(n_spectra, dtype=float)
    else:
        weights = np.asarray(weights)

    # Reconstruct the number of maps used in the cross-spectrum analysis, 
    # assuming n_spectra = n_maps*(n_maps-1)/2
    n_maps = np.int(np.ceil(np.sqrt(2*n_spectra)))

    # If we weren't told which map was in which cross-spectrum, reconstruct
    # our best guess. The covariance will be wrong if we're wrong about the
    # ordering of maps!  (Unless we ignore interaction terms -- JTS)
    if not indep_cov_only:
        if maps_used is None:
            core.log_warn(
                "You didn't tell me which cross-spectra use which maps. "
                "I'll need to make assumptions.")
            maps_used = [
                [x,y] for x in range(n_maps-1) for y in range(x+1,n_maps)]

        else:
            # If we got input maps_used, we can check that we're correct about 
            # how many maps we have. Assume that the number of maps is equal to
            # the maximum index. ... The rest of this function doesn't assume
            # that the map index is anything in particular, it only cares about
            # matching map indices used in making cross-spectra. 
            if np.max(maps_used)+1 != n_maps:
                core.log_warn(
                    "I think there should be %d maps that went into "%n_maps +
                    "this cross-spectrum, but the maps_used indicate %d. "%(
                        np.max(maps_used)+1) + "Problem?")

        # Standardize the "maps_used" input to be a (n_spectra, 2) array.
        maps_used = np.asarray( maps_used )
        if maps_used.ndim != 2:
            raise ValueError("The maps_used input must have two dimensions!")

        # Make sure that the "maps_used" array has N_spectra as the slow index.
        if maps_used.shape[0]==2:
            maps_used = np.transpose(maps_used)

        if not (maps_used.ndim==2 and maps_used.shape==(n_spectra,2)): 
            raise ValueError("The maps_used input doesn't look right! "
                             "It should have shape (%d,2)." % n_spectra)

    ###########################################################################
    # Divide the covariance into two parts. The first part is the noise
    # contribution assuming that each set of C_ells is independent (which is
    # obviously not true). The second part is the increase in the covariance
    # resulting from the fact that two cross-spectrum measurements which share
    # an observation are correlated. Assume that no cross-spectra in the input
    # C_ells are repeated (that is, the cross-spectrum of (map 1)x(map 2)
    # doesn't reappear as (map 2)x(map 1).
    ###########################################################################

    # First calculate the difference of each value of C_ells from the
    # (weighted) average.
    avg_cls = np.average(cls, axis=0, weights=weights)
    avg_cls_2d = np.resize( avg_cls, cls.shape )
    delta_cls = cls-avg_cls_2d

    # Multiply by this prefactor to normalize the sum.
    # The prefactor is equal to the usual (1 / (N-1)) if the weights are all 1.
    prefactor = np.sum(weights) / (np.sum(weights)**2 - np.sum(weights**2))
    # Reshape weights so we can multiply the delta_cls.
    weights = weights.reshape([len(weights),1])

    # Calculate the covariance under the assumption of independence.
    cov_independent = prefactor*np.dot(np.transpose(delta_cls),
                                       weights*delta_cls)
    # Return error on the (weighted) mean.
    # This is dividing by N if the weights are all 1.
    cov_independent /= (np.sum(weights)**2 / np.sum(weights**2))
    if indep_cov_only:
        return {'independent':cov_independent}

    ###########################################################################
    # Now calculate the second part, the increase in variance which comes from
    # correlation between cross-spectra. 
    #
    # To start, construct an (N_spectra x N_spectra) matrix.
    # Each row corresponds to a cross-spectrum made from a pair of maps, and
    # each column corresponds to a cross-spectrum from another pair of maps.
    # The matrix elements are 1 for (row, column) pairs where the
    # cross-spectra share one map. Since we never include the auto-spectrum,
    # the only time spectra share two maps is when we compare a cross-spectrum
    # to itself. This is the diagonal of this matrix. We've separated out 
    # this term as "cov_independent", so we skip the diagonal here.
    #
    # What we're doing is summing together the cross-spectrum / cross-spectrum
    # covariance terms. We assume that noise in cross-spectra is uncorrelated
    # (zero covariance) between cross-spectra which do not share any data.
    map_overlaps = np.zeros( [n_spectra, n_spectra] )

    for i_row, x in enumerate(maps_used):
        v = maps_used[i_row + 1:]
        i_col = np.where(np.isin(v, x).reshape(v.shape))[0] + i_row + 1
        map_overlaps[i_row, i_col] = map_overlaps[i_col, i_row] = 1

    # Now calculate the portion of the covariance coming from the correlations 
    # between different cross-spectra. This is the second term in eq. A15.
    # Note that instead of 
    #    2*Sum_{\beta \ne \lambda, \alpha} \Delta D^{\lambda \beta}, 
    # I use
    #    Sum_{\beta \ne \lambda, \alpha} \Delta D^{\lambda \beta} +
    #    Sum_{\beta \ne \lambda, \alpha} \Delta D^{\beta \alpha}.
    # I've also assumed that the cross-spectrum made from maps (A,B)
    # is the same as the cross-spectrum made from maps (B,A), and that
    # only one of those occurs in the list.
    #
    # This is the matrix multiplication
    # (Delta_Cls^T x map_overlaps x Delta_Cls), with shape
    # (n_bins, n_spectra) x (n_spectra, n_spectra) x (n_spectra, n_bins).
    # The "map_overlaps" term averages the delta_cls for cross-spectra with
    # map overlaps.
    cov_interactions = np.dot(np.transpose(delta_cls),
                              np.dot(map_overlaps, weights*delta_cls))
    cov_interactions *= prefactor
    cov_interactions /= (np.sum(weights)**2 / np.sum(weights**2))
   
    #########################################################################
    # Combine the two parts of the covariance to find the total covariance. 
    # If you want to compare this with IDL's unbiased_multiset_pspec.pro,
    # note that here I leave the diagonal elements of map_overlaps as zero. 
    # The IDL version includes them, and this is why it finishes
    # with cov = 2*cov2 - cov1: The diagonal elements are equal to cov1.
    cov_total = cov_independent + cov_interactions

    # Return the resulting covariance matrix.
    if total_cov_only:
        return cov_total
    else:
        return {'total':cov_total, 'independent':cov_independent,
                'interactions':cov_interactions}
