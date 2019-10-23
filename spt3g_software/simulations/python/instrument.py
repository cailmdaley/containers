"""
simulations/python/instrument.py

Properties of the instrument for sky simulations, including
beam response, white-noise level, and atmospheric noise.
"""
import numpy as np
import healpy as hp


def get_noise_sim_spectra(
    freq, component='instrument', lmax=3000, delta_t=None, lknee_t=None, lknee_p=None
):
    """
    Get instrumental or atmospheric noise power spectra.

    If left as None, the spectrum parameters will assume defaults
    based on the 5-year SPT-3G projected levels for frequency `freq`.
    Default output units are uK^2

    Parameters
    ----------
    freq : int
        Frequency band. Must be one of 90, 150, or 220
    component : str
        Either 'instrument', or 'atmosphere'
    lmax : int
        max multipole for Nls
    delta_t : float
        Temperature map depth in uK-arcmin
    lknee_t : float
        Where temperature 1/f flattens out to white noise
    lknee_p : float
        Where polarization 1/f flattens out to white nosie

    Returns
    -------
    nls_t : array
        Temperature noise angular power spectrum
    nls_p : array
        Polarization noise angular power spectrum

    See Also
    --------
    get_nl
    """
    assert component == 'instrument' or component == 'atmosphere'

    # default values are SPT-3G 5-year noise values
    if freq == 90:
        if delta_t is None:
            delta_t = 3.0
        if lknee_t is None:
            lknee_t = 1200
        if lknee_p is None:
            lknee_p = 300
    elif freq == 150:
        if delta_t is None:
            delta_t = 2.2
        if lknee_t is None:
            lknee_t = 2200
        if lknee_p is None:
            lknee_p = 300
    elif freq == 220:
        if delta_t is None:
            delta_t = 8.8
        if lknee_t is None:
            lknee_t = 2300
        if lknee_p is None:
            lknee_p = 300
    else:
        raise ValueError('freq must be one of 90, 150, 220')

    # white noise component from instrument
    nls_t = get_nl(delta_t, lmax)
    nls_p = get_nl(delta_t * np.sqrt(2.0), lmax)

    if component == 'instrument':
        return nls_t, nls_p

    if component == 'atmosphere':
        # atmosphere 1/f component
        ells = np.arange(len(nls_t))
        nls_t_atm = nls_t * (1.0 + (lknee_t / ells))
        nls_p_atm = nls_p * (1.0 + (lknee_p / ells))

        nls_t_atm[ells > lknee_t] = 0.0
        nls_p_atm[ells > lknee_p] = 0.0

        nls_t_atm[np.invert(np.isfinite(nls_t_atm))] = 0
        nls_p_atm[np.invert(np.isfinite(nls_p_atm))] = 0
        return nls_t_atm, nls_p_atm


def get_nl(delta_t, lmax, ells=None, output_units='uk'):
    """
    Returns white-noise power spectrum up to lmax

    Parameters
    ----------
    delta_t : float
        White-noise level in uK-arcmin
    lmax : int
        Maximum ell up to which to compute the spectrum
    ells : array_like, optional
        If specified, only return values at these ells
    output_units : str
        If set to 'uk', Nls will be in uK^2. Otherwise K^2

    Returns
    -------
    nl : array
        White-noise power spectrum
    """
    if output_units.lower() == 'k':
        delta_t /= 1e6

    delta_t_radians = delta_t * np.radians(1.0 / 60.0)
    nl = np.full(int(lmax) + 1, delta_t_radians ** 2.0)

    if ells is not None:
        nl = nl[ells.astype(int)]

    return nl


def get_beams(fwhm_90=1.7, fwhm_150=1.2, fwhm_220=1.0, beamfile=None, lmax=2e4):
    """
    Get the Bl for each frequency

    Uses beamfile if supplied, otherwise makes gaussian beams

    Parameters
    ----------
    fwhm_ : float
        The beam FWHM, in arcminutes, at the corresponding frequency band.
    beamfile : str
        Filepath to .txt file containing beams
    lmax : int
        maximum ell up to which to compute the beam

    Returns
    -------
    bls : dict
        Dictionary of Bl keyed by observing frequency.
    """
    bls = {}
    if beamfile not in ['', None]:
        ell, b90, b150, b220 = np.loadtxt(beamfile, unpack=True)
        bls[90] = b90[ell < lmax + 1]
        bls[150] = b150[ell < lmax + 1]
        bls[220] = b220[ell < lmax + 1]
    else:
        fwhms = {90: fwhm_90, 150: fwhm_150, 220: fwhm_220}
        for freq, fwhm in fwhms.items():
            fwhm_rad = np.radians(fwhm / 60.0)
            bls[freq] = hp.gauss_beam(fwhm_rad, lmax)

    return bls
