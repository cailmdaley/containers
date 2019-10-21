# null_tests.py
#
# For each type of null test, need to 
# 1) find each ObsIds value for that test (if applicable)
# 2) rank the bundles according to the ObsIds they contain
# 3) Difference the bundles, producing null maps
# 4) Take cross-spectra, etc.

import os
import numpy as np
from spt3g import core, std_processing, mapmaker
from spt3g.mapspectra.map_analysis import subtract_two_maps
from spt3g.util import files
from glob import glob
import astropy
import scipy.stats
import argparse

def calculate_cl_chi_square(input_cls,
                            spectrum,
                            use_sigma_errors=True, 
                            expected_cls=None,
                            ell_range=[299, 3001]):

    """
    Calculate the chi^2 and PTEs for an input set of C_ells, compared to an
    identically binned set of expected C_ells. This was made with null tests
    in mind, but could be used for anything.
    
    Parameters
    ----------
    input_cls : array
        C_ells, as output from average_cls.
        Must have uncertainties, preferably as an "error_XX" field.

    spectrum : str
        Which spectrum should we look at (e.g. "EE", "TT")?

    use_sigma_errors : bool
        If True, use std(cross_spectra)/sqrt(n_spectra) as the errors.
        In the limit where the input_cls are noise dominated (e.g. for null tests),
        this is the same, in the mean, as the true diagonal covariance matrix elements,
        and is much less noisy.

    expected_cls : array
        C_ells, containing at least the 'ell' and spectrum field.
        Must be binned identically to input_cls.
        If None, assume the expected values are zero.

    ell_range : 2-element list
        Use only ell bins with central values in this range for chi^2 calculation.

    Returns
    -------
    A 3-tuple of (chi_sq, ndf, pte).

    Notes
    -----
    Stephen Hoover, 27 January 2014 (calculateJackknifeChiSquare)
    Barely modified by DPD for SPT-3G, April 2019
    """
    # Figure out the uncertainties on the input spectrum.
    if 'error_'+spectrum in input_cls and not use_sigma_errors:
        input_error = input_cls['error_'+spectrum]
    elif 'sigma_'+spectrum in input_cls:
        if not use_sigma_errors:
            print("[calculateClChiSquare] "
                  "Using sigma errors!")
        input_error = input_cls['sigma_'+spectrum]
    else:
        raise ValueError("I can't determine the uncertainty on the input Cls.")

    # If we didn't get any expected values, then set them to zero.
    if expected_cls is None:
        expected_cls = {'ell' : input_cls[spectrum].bin_centers,
                        spectrum : np.zeros_like(input_cls[spectrum])}
    elif not (expected_cls['ell'][spectrum]==input_cls[spectrum].bin_centers).all():
        raise ValueError("The expected C_ells are binned differently "
                         "than the input spectrum!")

    # Create a masking array for the ell bins we want to use.
    ell = input_cls[spectrum].bin_centers
    ell_to_use = (ell >= ell_range[0]) & (ell <= ell_range[1])

    # Calculate the chi^2.
    chi_sq = np.sum((
        (input_cls[spectrum][ell_to_use] - expected_cls[spectrum][ell_to_use]) /
                      input_error[ell_to_use])**2 )

    # Calculate the probability to exceed the calculated chi^2.
    ndf = np.sum(ell_to_use)
    pte = 1 - scipy.stats.chi2.cdf(chi_sq, ndf)

    return chi_sq, ndf, pte


def make_all_null_maps(bundle_dir='/spt/user/ddutcher/bundles/150GHz/total',
                       bundle_def='/spt/user/ddutcher/bundles/bundle_obsids.pkl',
                       az_bundle_def='/spt/user/ddutcher/bundles/az_bundle_obsids.pkl',
                       tests_to_do=[],
                       output_dir='/spt/user/ddutcher/null_maps/150GHz',
                       verbose=False):
    """
    Note: assumes bundles are labelled like 'bun00_...', 'bun01_...', etc.


    Parameters:
    ----------
    bundle_dir : str
        Parent directory containing bundles out of which to make null maps.
        Contains subdirectories for each type of bundle
        ('total', 'left_right', 'high_low', 'azimuth')
    bundle_def : str or dict
        Dictionary mapping bundle index to list of observation ids that it contains.
        The keys of this dict are used to map a bundle path to its bundle index.
    tests_to_do : list
        List of strings of null tests for which to make maps.
        Current options are:
            'left_right': left-going scans vs right-going
            '1_2'       : First-half vs second-half
            'azimuth'   : noisy ground vs cleaner ground
            'sun'       : Sun up vs Sun down
            'moon'      : Moon up vs Moon down
            'saturation': # of saturated bolos
            'wafer'     : Low vs High responsivity wafers
        Note 'azimuth', 'left-right', and 'wafer' tests use specialized bundles.
    output_dir : str
        Parent directory to store null map directroies. Each null test has
        its own subdirectory created.
    """
    if not isinstance(tests_to_do, list):
        tests_to_do = [tests_to_do]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def _find_bundles(bundle_dir, filepth='*.g3*', bundle_def=None): 
        # Define the bundles we're going to difference
        bundles = sorted(glob(os.path.join(bundle_dir, filepth)))
        assert len(bundles) == 30
        if len(bundles) == 0:
            raise FileNotFoundError(os.path.join(bundle_dir, filepth))
        if isinstance(bundle_def, str):
            bundle_def = files.load_pickle(bundle_def)
        if not isinstance(bundle_def, dict):
            raise TypeError('bundle_def should be dictionary')

        label_pth = dict()
        for bundle_label in bundle_def.keys():
            for bun_pth in bundles:
                if 'bun%02d'%int(bundle_label)+'_' in bun_pth:
                    label_pth[bundle_label] = bun_pth
        return label_pth

    null_sets = dict()
    for test in tests_to_do:
        null_sets[test] = []
        if test == 'left_right':
            print('%s test requires special bundles; hope you made them!'%test)
            left = sorted(glob(os.path.join(
                bundle_dir, 'left_right', '*left*.g3*')))
            right = sorted(glob(os.path.join(
                bundle_dir, 'left_right', '*right*.g3*')))
            if len(left) != len(right):
                raise AssertionError("Unequal number of left and right maps!")
            for i in np.arange(0, len(left)):
                null_sets[test].append([left[i], right[i]])

        elif test == '1_2':
            # bundling is chronological, so skip the metric business
            # and just use bundle label.
            label_pth = _find_bundles(
                bundle_dir, 'total/*.g3*', bundle_def=bundle_def)
            set1 = np.sort(np.array(list(label_pth.keys())))
            set2 = np.sort(np.array(list(label_pth.keys())))[::-1]
            for i in np.arange(0, int(len(set1)/2)):
                null_sets[test].append([label_pth[set1[i]], label_pth[set2[i]]])

        elif test == 'azimuth':
            print('%s test requires special bundles; hope you made them!'%test)
            label_pth = _find_bundles(
                bundle_dir, 'azimuth/*.g3*', bundle_def=az_bundle_def)
            set1 = np.sort(np.array(list(label_pth.keys())))
            set2 = np.sort(np.array(list(label_pth.keys())))[::-1]
            for i in np.arange(0, int(len(set1)/2)):
                null_sets[test].append([label_pth[set1[i]], label_pth[set2[i]]])

        elif test == 'sun':
            label_pth = _find_bundles(
                bundle_dir, 'total/*.g3*', bundle_def=bundle_def)
            metric = get_generic_metric(
                bundle_def, 
                config_file = '/spt/user/ddutcher/null_configs/sun_up.pkl')
            ordering = np.array(list(metric.values())).argsort()
            set1 = np.array(list(metric.keys()))[ordering]
            set2 = np.array(list(metric.keys()))[ordering][::-1]
            for i in np.arange(0, int(len(set1)/2)):
                null_sets[test].append([label_pth[set1[i]], label_pth[set2[i]]])

        elif test == 'moon':
            label_pth = _find_bundles(
                bundle_dir, 'total/*.g3*', bundle_def=bundle_def)
            metric = get_generic_metric(
                bundle_def, 
                config_file = '/spt/user/ddutcher/null_configs/moon_up.pkl')
            ordering = np.array(list(metric.values())).argsort()
            set1 = np.array(list(metric.keys()))[ordering]
            set2 = np.array(list(metric.keys()))[ordering][::-1]
            for i in np.arange(0, int(len(set1)/2)):
                null_sets[test].append([label_pth[set1[i]], label_pth[set2[i]]])

        elif test == 'saturation':
            label_pth = _find_bundles(
                bundle_dir, 'total/*.g3*', bundle_def=bundle_def)
            metric = get_generic_metric(
                bundle_def, 
                config_file = '/spt/user/ddutcher/null_configs/saturation.pkl')
            ordering = np.array(list(metric.values())).argsort()
            set1 = np.array(list(metric.keys()))[ordering]
            set2 = np.array(list(metric.keys()))[ordering][::-1]
            for i in np.arange(0, int(len(set1)/2)):
                null_sets[test].append([label_pth[set1[i]], label_pth[set2[i]]])

        elif test == 'wafer':
            print('%s test requires special bundles; hope you made them!'%test)
            low = sorted(glob(os.path.join(
                bundle_dir, 'high_low', '*low*.g3*')))                
            high = sorted(glob(os.path.join(
                bundle_dir, 'high_low', '*high*.g3*')))
            if len(low) != len(high):
                raise AssertionError("Unequal number of wafer maps!")
            for i in np.arange(0, len(low)):
                null_sets[test].append([low[i], high[i]])

        else:
            print('Unrecognized test %s, skipping.'%test)

    # Make null maps
    for test, sets in null_sets.items():
        if verbose:
            print('Making %s null maps'%test)
        if len(sets) == 0:
            continue
        if not os.path.exists(os.path.join(output_dir, test)):
            os.mkdir(os.path.join(output_dir, test))
        for ind, (pth1, pth2) in enumerate(sets):
            if verbose:
                print('Making %s null map %s'%(test, ind))
            label1 = pth1.split('/')[-1].partition('.')[0]
            label2 = pth2.split('/')[-1].partition('.')[0]
            difference_map = subtract_two_maps(pth1, pth2, in_place=True)
            difference_map['Id'] = '%s-%s'%(label1, label2)
            filename = os.path.join(output_dir,test,
                                    'null_'+test+'_%02d'%ind+'.g3.gz')
            core.G3Writer(filename)(difference_map)

# =============================================================================
# Getting split metrics: assign one number to each bundle based on its obs.
# =============================================================================

def get_1_2_metric(bundle_def):
    """
    Figures out which observations were in the first half and which
    were in the second half, chronologically.
    Then assigns 0 to bundles purely in the first half, 1 to bundles
    purely in the second half, and a number in between for others.
    """
    # This is a little redundant as the bundles were
    # made chronologically.
    chronologic_metric = dict()
    if isinstance(bundle_def, str):
        bundle_def = files.load_pickle(bundle_def)
    all_obs = []
    for bundle, obsids in bundle_def.items():
        for obsid in obsids:
            all_obs.append(int(obsid))
    all_obs = np.sort(all_obs)
    midpoint = all_obs[int(len(all_obs)/2)]
    for bundle, obsids in bundle_def.items():
        chronologic_metric[bundle] = 0
        for obsid in obsids:
            if int(obsid) >= midpoint:
                chronologic_metric[bundle] += 1
        chronologic_metric[bundle] /= len(obsids)
    return chronologic_metric


def get_azimuth_metric(bundle_def, config_file):
    """
    Note: I did not end up using this. Instead I re-bundled
        based on azimuth, and used those az_bundles.

    Azimuth positions are defined as good(0) or bad(1)
    based on the known locations of buildings at the South Pole.
    Ground maps are too noisy to reliably extract actual
    information about ground contamination.
    """
    bad_az_range = [95, 275]
    if isinstance(bundle_def, str):
        bundle_def = files.load_pickle(bundle_def)
    config = files.load_pickle(config_file)
    metric = dict()
    for bundle, obsids in bundle_def.items():
        metric[bundle] = 0
        for obsid in obsids:
            dist =  np.abs(153-config[int(obsid)])
            if dist > 180:
                dist = 180 - (dist % 180)
            metric[bundle] += dist
        metric[bundle] /= len(obsids)
    return metric


def get_generic_metric(bundle_def, config_file):
    """
    Passed a dictionary mapping bundle id to the list of
    observations it contains, and another dictionary mapping each
    observation to a quantity, this will return the average quantity
    in each bundle.
    """
    if isinstance(bundle_def, str):
        bundle_def = files.load_pickle(bundle_def)
    config = files.load_pickle(config_file)
    metric = dict()
    for bundle, obsids in bundle_def.items():
        metric[bundle] = 0
        for obsid in obsids:
            metric[bundle] += config[int(obsid)]
        metric[bundle] /= len(obsids)
    return metric


# =============================================================================
# Config files: getting observation-specific quantities
# =============================================================================

def make_moon_config(filename='/spt/user/ddutcher/null_configs/moon_up.pkl',
                     source='ra0hdec-??.??',
                     verbose=True):
    """
    For every 1500d observation, this will figure out if the moon was up (1)
    or down (0).
    """
    if os.path.isfile(filename):
        out = files.load_pickle(filename)
    else:
        out = dict()
    obs = sorted(glob('/spt/data/bolodata/downsampled/'+source+'/*'))
    for pth in obs:
        obsid = int(pth.split('/')[-1])
        if obsid in out.keys():
            continue
        if verbose:
            print(obsid)
        tt = std_processing.obsid_to_g3time(int(obsid))
        moon_pos = std_processing.sourcemaps.get_source_ra_dec('moon',
                                                               at_time=tt)
        if moon_pos[1] < 0:
            out[obsid] = 1 # 'up'
        else:
            out[obsid] = 0 #' down'
    files.save_pickle(out, filename)


def make_sun_config(filename='/spt/user/ddutcher/null_configs/sun_up.pkl',
                    source='ra0hdec-??.??',
                    verbose=True):
    """
    For every 1500d observation, this will figure out if the sun was up (1)
    or down (0).
    """
    if os.path.isfile(filename):
        out = files.load_pickle(filename)
    else:
        out = dict()
    obs = sorted(glob('/spt/data/bolodata/downsampled/'+source+'/*'))
    for pth in obs:
        obsid = int(pth.split('/')[-1])
        if obsid in out.keys():
            continue
        if verbose:
            print(obsid)
        tt = std_processing.obsid_to_g3time(int(obsid))
        sun_pos = std_processing.sourcemaps.get_source_ra_dec('sun',
                                                              at_time=tt)
        if sun_pos[1] < 0:
            out[obsid] = 1 # 'up'
        else:
            out[obsid] = 0 #' down'
    files.save_pickle(out, filename)


def make_saturated_config(filename='/spt/user/ddutcher/null_configs/saturation.pkl'):
    """
    This stores the mean number of bolometers flagged as Overbiased per scan
    for each observation for which a map exists.
    Right now uses pre-computed statistics.
    """
    out = dict()
    ob_flags = files.load_pickle(
        '/home/ddutcher/data/2018/flags/ob_vs_scan.pkl')
    for obsid, arr in ob_flags.items():
        out[int(obsid)] = np.mean(arr)
    files.save_pickle(out, filename)


def make_azimuth_config(filename='/spt/user/ddutcher/null_configs/azimuth.pkl',
                        source='ra0hdec-??.??',
                        verbose=True):
    """
    Return the average azimuth for each observation, in range [0,360).
    """
    spt = astropy.coordinates.EarthLocation(lat=-89.991066*astropy.units.deg,
                                            lon=-44.65*astropy.units.deg,
                                            height=2835.0*astropy.units.meter)
    if os.path.isfile(filename):
        out = files.load_pickle(filename)
    else:
        out = dict()
    obs = sorted(glob('/spt/data/bolodata/downsampled/'+source+'/*/0000.g3'))
    # Extracting BoresightAz for every scan of every observation is insane.
    # Instead, get source Az at the start and stop of observation.
    for pth in obs:
        obsid = int(pth.split('/')[-2])
        if obsid in out.keys():
            continue
        if verbose:
            print(obsid)
        dat = core.G3File(pth)
        for frame in dat:
            if frame.type == core.G3FrameType.Observation:
                break
        try:
            obs_src = frame['SourceName']
            times = [frame['ObservationStart'], frame['ObservationStop']]
        except KeyError:
            continue
            
        ra0, dec0 = std_processing.sourcemaps.get_source_ra_dec(
            obs_src, at_time = times[0])
        ra1, dec1 = std_processing.sourcemaps.get_source_ra_dec(
            obs_src, at_time = times[1])
        ra = [ra0, ra1]
        dec = [dec0, dec1]

        # The below borrows from coordinateutils.azel.convert_radec_to_azel
        t = astropy.time.Time(np.asarray([i.mjd for i in times]), format='mjd')
        k = astropy.coordinates.FK5(
            ra=np.asarray(ra)/core.G3Units.deg*astropy.units.deg,
            dec=np.asarray(dec)/core.G3Units.deg*astropy.units.deg)
        kt = k.transform_to(astropy.coordinates.AltAz(
            obstime=t, location=spt, pressure=0))
        az = kt.az/astropy.units.deg
        if np.abs(az[0] - az[1]) > 180.:
            for i, aa in enumerate(az):
                if aa > 180:
                    az[i] = aa - 360.
        avg = (az[0]+az[1])/2.
        if avg < 0:
            avg += 360.
        elif avg > 360.:
            avg -= 360.
        out[obsid] = np.float(avg)
    files.save_pickle(out, filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bundle_dir')
    parser.add_argument('--output', '-o')
    args = parser.parse_args()

    for band in ['90GHz', '150GHz', '220GHz']:
        make_all_null_maps(bundle_dir=os.path.join(args.bundle_dir, band),
                           bundle_def='/spt/user/ddutcher/bundles/bundle_obsids.pkl',
                           az_bundle_def='/spt/user/ddutcher/bundles/az_bundle_obsids.pkl',
                           output_dir=os.path.join(args.output, band),
                           tests_to_do=[
                               '1_2', 'azimuth', 'left_right', 'moon', 'saturation', 'wafer'],
                           verbose=True)
