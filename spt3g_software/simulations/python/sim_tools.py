'''
simulations/python/sim_tools.py

Functions for combining and saving sims.
'''
import os
import numpy as np
import healpy as hp
from spt3g import core
from spt3g.util import healpix_tools as hpt


def save_sims(
    maps=None,
    alms=None,
    verbose=False,
    map_out='sim_map.fits',
    alm_out='sim_alm.fits',
    store_alm=True,
    store_map=False,
    lmax=None,
    nside=None,
    pol=True,
):
    """
    Save simulated healpix maps or alms.
    
    Can convert alms to maps or maps to alms if required.
    """
    if maps is None and alms is None:
        raise ValueError("Must specify either maps or alms to save.")
    map_out = os.path.abspath(map_out)
    alm_out = os.path.abspath(alm_out)
    if isinstance(alms, np.ndarray) and len(alms.shape) > 1:
        alms = list(alms)
    if maps is not None and alms is not None:
        if store_map:
            hp.write_map(map_out, maps)
        if store_alm:
            hp.write_alm(alm_out, alms)
    elif maps is not None:
        if store_map:
            hp.write_map(map_out, maps)
        if store_alm:
            # get the alms
            curralm = hp.map2alm(maps, lmax=lmax)
            hp.write_alm(alm_out, curralm)
    elif alms is not None:
        if store_alm:
            hp.write_alm(alm_out, alms)
        if store_map:
            # get the maps
            currmap = hp.alm2map(
                alms, nside, lmax=lmax, pol=pol, inplace=True, verbose=verbose
            )
            hp.write_map(map_out, currmap)


def combine_alms_into_map(
    alms_to_add,
    nside=2048,
    lmax=2e4,
    pol=True,
    freq=150,
    add_beam=True,
    beamfile=None,
    fwhm_90=1.7,
    fwhm_150=1.2,
    fwhm_220=1.0,
    verbose=False,
):
    """
    Combine all the alms made into a single map

    Parameters
    ----------
    alms_to_add : list of strings
        List of filepaths corresponding to maps to be added
    pol : bool
        make polarised skies. default is True.
    freq : int
        Frequency channel of alms. Used to get proper beam.
    add_beam : bool
        If True add beam. It invokes get_beams to get the desired beam
    beamfile : str
        Filepath to .txt file containing Bl. Overrides `fwhm_*`
    fwhm_ : float
        The beam FWHM, in arcminutes, at the corresponding frequency band.
        Used to compute Bl if `add_beam`=True but `beamfile` unspecified.
    lmax : int
        for alms if beam needs to be applied
    nside : int
        for generating final map

    Returns
    -------
    combined_map : array or sequence of three arrays
        Map containing simulated CMB (+ FG + noise).
        If `pol`=True, the arrays will corresond to T,Q,U.

    See Also
    --------
    get_beams
    """
    combined_alm = None
    for curralmname in alms_to_add:
        if pol:
            curralm = [0, 1, 2]
        else:
            curralm = [0]
        for teb in range(len(curralm)):
            curralm[teb] = hp.read_alm(curralmname, hdu=int(teb + 1))
        curralm = np.asarray(curralm)
        if combined_alm is not None:
            combined_alm += curralm
        else:
            combined_alm = np.copy(curralm)

    if add_beam:
        core.log_info('Applying beam to combined alms')
        from .instrument import get_beams

        bls = get_beams(
            fwhm_90=fwhm_90,
            fwhm_150=fwhm_150,
            fwhm_220=fwhm_220,
            beamfile=beamfile,
            lmax=lmax,
        )
        bls = bls[freq]
        for tqu in range(len(combined_alm)):
            combined_alm[tqu] = hp.almxfl(combined_alm[tqu], bls, inplace=True)

    # Get the map back with the same nside as the inputs
    combined_map = hp.alm2map(list(combined_alm), nside, lmax=lmax, verbose=verbose)
    combined_map = np.asarray(combined_map)

    return combined_map


def save_healpix_as_spt3g_map(
    hp_map, filename, maskfile=None, verbose=False, units=core.G3Units.uK
):
    '''
    Save map as .fits that can be read by coordinateutils.fitsio.load_skymap_fits

    Parameters
    ----------
    hp_map : str, array, or sequence of 3 arrays
        The healpix map to be stored
    filename : str
        The fits file name
    maskfile : str
        SPT-3G survey area mask, as fits file
    verbose : bool
        print or not print healpy terminal outputs
    units: float
        The corresponding G3Units value of the physical units of
        the input map. Used to put output map in G3Units system.

    See Also
    --------
    util.healpix_tools.write_map
    '''
    if isinstance(hp_map, str):
        t = hp.read_map(hp_map, field=0, verbose=verbose)
        try:
            q = hp.read_map(hp_map, field=1, verbose=verbose)
            u = hp.read_map(hp_map, field=2, verbose=verbose)
            column_names = ['T', 'Q', 'U']
            map_arr = np.asarray([t, q, u])
        except:
            map_arr = t
            column_names = ['T']
    else:
        map_arr = np.asarray(hp_map)
        if len(np.shape(hp_map)) == 1:
            column_names = ['T']
        elif len(hp_map) == 3:
            column_names = ['T', 'Q', 'U']
        else:
            raise ValueError("Input healpix map has unrecognized shape.")

    # Convert output to G3Units
    map_arr *= units

    nside = hp.npix2nside(np.shape(map_arr)[-1])
    if maskfile not in [None, '']:
        mask = hp.ud_grade(hp.read_map(maskfile, verbose=verbose, dtype=bool), nside)
    else:
        mask = None

    if len(column_names) == 3:
        map_arr = list(map_arr)

    hpt.write_map(
        filename,
        map_arr,
        coord='C',
        mask=mask,
        column_names=column_names,
        extra_header={'WEIGHTED': False, 'UNITS': 'Tcmb'},
        partial=True,
    )
