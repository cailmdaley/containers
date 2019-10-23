"""
A collection of tools for handling healpix maps.

 * Wrappers around healpy functions with some extended functionality (mostly masking)
 * Python wrappers for the compiled Fortran executables
 * Functions for creating masks and outlines
 * Functions for rotating T/Q/U maps and pointing matrix maps between coordinate systems
   in a way that properly preserves polarization power
 * Functions for stacking hot and cold spots
"""
import numpy as np
import healpy as hp
import shutil
import tempfile as tf
import subprocess as sp
import os
from spt3g import core

__all__ = ['synfast', 'write_alm', 'write_map', 'read_map', 'rotate_map',
           'rotate_proj', 'alteralm', 'rotate_map_alm', 'rotate_qu',
           'outline_map', 'latlon_mask', 'circle_mask', 'oval_mask', 'mask_bad',
           'mask_good', 'smoothalm', 'smoothing', 'anafast', 'fill_outline',
           'hotspot', 'stack_map', 'rotate_qu_radial', 'rotate_qu_stack',
           'ud_grade', 'spice']

synfast_exe = 'synfast'
"""
default synfast executable.  Maybe do some smart hostname lookup later.
"""

alteralm_exe = 'alteralm'
"""
default alteralm executable.  Maybe do some smart hostname lookup later.
"""

hotspot_exe = 'hotspot'
"""
default hotspot executable.  Maybe do some smart hostname lookup later.
"""

spice_exe = 'spice'
"""
default spice executable.  Maybe do some smart hostname lookup later.
"""

@core.usefulfunc
def mask_bad(m, badval=hp.UNSEEN, rtol=1.e-5, atol=1.e-8,
             badnan=True, badinf=True):
    """Returns a bool array with ``True`` where m is close to badval,
    NaN or inf.

    Parameters
    ----------
    m : array or sequence of arrays
        A map (may be a sequence of maps)
    badval : float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol : float, optional
        The relative tolerance
    atol : float, optional
        The absolute tolerance
    badnan : bool, optional
        If True, also mask NaN values
    badinf : bool, optional
        If True, also mask inf values

    Returns
    -------
    mask
      a bool array with the same shape as the input map, ``True`` where input map is
      close to badval, NaN or inf, and ``False`` elsewhere.

    See Also
    --------
    mask_good

    Examples
    --------
    >>> import healpy as hp
    >>> import numpy as np
    >>> m = np.arange(12.)
    >>> m[3] = hp.UNSEEN
    >>> m[4] = np.nan
    >>> mask_bad(m)
    array([False, False, False,  True,  True, False, False, False, False,
           False, False, False], dtpe=bool)
    """
    m = np.asarray(m)
    mask = np.zeros_like(m, dtype=bool)
    if badnan:
        mask |= np.isnan(m)
    if badinf:
        mask |= np.isinf(m)
    mask[~mask] = hp.mask_bad(m[~mask], badval=badval, rtol=rtol, atol=atol)
    return mask

@core.usefulfunc
def mask_good(m, badval=hp.UNSEEN, rtol=1.e-5, atol=1.e-8,
              badnan=True, badinf=True):
    """Returns a bool array with ``False`` where m is close to badval,
    NaN or inf.

    Parameters
    ----------
    m : array or sequence of arrays
        A map (may be a sequence of maps)
    badval : float, optional
        The value of the pixel considered as bad (:const:`UNSEEN` by default)
    rtol : float, optional
        The relative tolerance
    atol : float, optional
        The absolute tolerance
    badnan : bool, optional
        If True, also mask NaN values
    badinf : bool, optional
        If True, also mask inf values

    Returns
    -------
    a bool array with the same shape as the input map, ``False`` where input map is
    close to badval, NaN or inf, and ``True`` elsewhere.

    See Also
    --------
    mask_bad, ma

    Examples
    --------
    >>> import healpy as hp
    >>> m = np.arange(12.)
    >>> m[3] = hp.UNSEEN
    >>> m[4] = np.nan
    >>> mask_good(m)
    array([ True,  True,  True, False, False,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)
    """
    m = np.asarray(m)
    mask = np.ones_like(m, dtype=bool)
    if badnan:
        mask &= ~np.isnan(m)
    if badinf:
        mask &= np.isfinite(m)
    mask[mask] = hp.mask_good(m[mask], badval=badval, rtol=rtol, atol=atol)
    return mask


# Copied from healpy commit #dbd1131
#FIXME remove this and require healpy version > 1.8.5 when released
@core.usefulfunc
def write_alm(filename, alms, out_dtype=None, lmax=-1, mmax=-1, mmax_in=-1):
    """Write alms to a fits file.

    In the fits file the alms are written
    with explicit index scheme, index = l*l + l + m +1, possibly out of order.
    By default write_alm makes a table with the same precision as the alms.
    If specified, the lmax and mmax parameters truncate the input data to
    include only alms for which l <= lmax and m <= mmax.

    Parameters
    ----------
    filename : str
      The filename of the output fits file
    alms : array, complex or list of arrays
      A complex ndarray holding the alms, index = m*(2*lmax+1-m)/2+l, see Alm.getidx
    lmax : int, optional
      The maximum l in the output file
    mmax : int, optional
      The maximum m in the output file
    out_dtype : data type, optional
      data type in the output file (must be a numpy dtype). Default: *alms*.real.dtype
    mmax_in : int, optional
      maximum m in the input array
    """

    from healpy import Alm
    from healpy.fitsfunc import getformat, pf
    from healpy import cookbook as cb

    if not cb.is_seq_of_seq(alms):
        alms = [alms]

    l2max = Alm.getlmax(len(alms[0]),mmax=mmax_in)
    if (lmax != -1 and lmax > l2max):
        raise ValueError("Too big lmax in parameter")
    elif lmax == -1:
        lmax = l2max

    if mmax_in == -1:
        mmax_in = l2max

    if mmax == -1:
        mmax = lmax
    if mmax > mmax_in:
        mmax = mmax_in

    if (out_dtype == None):
        out_dtype = alms[0].real.dtype

    l,m = Alm.getlm(lmax)
    idx = np.where((l <= lmax)*(m <= mmax))
    l = l[idx]
    m = m[idx]

    idx_in_original = Alm.getidx(l2max, l=l, m=m)

    index = l**2 + l + m + 1

    hdulist = pf.HDUList()
    for alm in alms:
        out_data = np.empty(len(index),
                             dtype=[('index','i'),
                                    ('real',out_dtype),
                                    ('imag',out_dtype)])
        out_data['index'] = index
        out_data['real'] = alm.real[idx_in_original]
        out_data['imag'] = alm.imag[idx_in_original]

        cindex = pf.Column(name="index", format=getformat(np.int32),
                           unit="l*l+l+m+1", array=out_data['index'])
        creal = pf.Column(name="real", format=getformat(out_dtype),
                          unit="unknown", array=out_data['real'])
        cimag = pf.Column(name="imag", format=getformat(out_dtype),
                          unit="unknown", array=out_data['imag'])

        tbhdu = pf.BinTableHDU.from_columns([cindex,creal,cimag])
        hdulist.append(tbhdu)
    from astropy import __version__ as astropy_version
    if astropy_version >= "1.3":
        hdulist.writeto(filename, overwrite=True)
    else:
        hdulist.writeto(filename, clobber=True)

@core.usefulfunc
def synalm(cls, lmax=None, mmax=None, new=True, verbose=True,
           seed=None, prng=None):
    """Generate a set of alm given cl.
    The cl are given as a float array. Corresponding alm are generated.
    If lmax is None, it is assumed lmax=cl.size-1
    If mmax is None, it is assumed mmax=lmax.

    Parameters
    ----------
    cls : float, array or tuple of arrays
      Either one cl (1D array) or a tuple of either 4 cl
      or of n*(n+1)/2 cl.
      Some of the cl may be None, implying no
      cross-correlation. See *new* parameter.
    lmax : int, scalar, optional
      The lmax (if None or <0, the largest size-1 of cls)
    mmax : int, scalar, optional
      The mmax (if None or <0, =lmax)
    new : bool, optional
      If True, use the new ordering of cl's, ie by diagonal
      (e.g. TT, EE, BB, TE, EB, TB or TT, EE, BB, TE if 4 cl as input).
      If False, use the old ordering, ie by row
      (e.g. TT, TE, TB, EE, EB, BB or TT, TE, EE, BB if 4 cl as input).
    seed : scalar, optional
      Random number generator seed.
    prng : numpy.random.RandomState instance, optional
      If supplied, use this RandomState instance to generate maps.
      Overrides the `seed` option.

    Returns
    -------
    alms : array or list of arrays
      the generated alm if one spectrum is given, or a list of n alms
      (with n(n+1)/2 the number of input cl, or n=3 if there are 4 input cl).

    Notes
    -----
    The order of the spectra will change in a future release. The new= parameter
    help to make the transition smoother. You can start using the new order
    by setting new=True.
    In the next version of healpy, the default will be new=True.
    This change is done for consistency between the different tools
    (alm2cl, synfast, anafast).
    In the new order, the spectra are ordered by diagonal of the correlation
    matrix. Eg, if fields are T, E, B, the spectra are TT, EE, BB, TE, EB, TB
    with new=True, and TT, TE, TB, EE, EB, BB if new=False.
    """
    import healpy._healpy_sph_transform_lib as sphtlib

    if (not new) and verbose:
        warnings.warn("The order of the input cl's will change in a future "
                      "release.\n"
                      "Use new=True keyword to start using the new order.\n"
                      "See documentation of healpy.synalm.",
                      category=hp.sphtfunc.FutureChangeWarning)
    if not hp.cookbook.is_seq(cls):
        raise TypeError('cls must be an array or a sequence of arrays')

    # store state and re-seed random number generator if requested
    if seed is not None and prng is None:
        np.random.seed(seed)
        prng = np.random
    elif prng is None:
        prng = np.random.RandomState()

    if not hp.cookbook.is_seq_of_seq(cls):
        # Only one spectrum
        if lmax is None or lmax < 0:
            lmax = cls.size-1
        if mmax is None or mmax < 0:
            mmax = lmax
        cls_list = [np.asarray(cls, dtype=np.float64)]
        szalm = hp.Alm.getsize(lmax,mmax)
        alm = np.zeros(szalm,'D')
        alm.real = prng.standard_normal(szalm)
        alm.imag = prng.standard_normal(szalm)
        alms_list=[alm]
        sphtlib._synalm(cls_list,alms_list,lmax,mmax)
        return alm

    # From here, we interpret cls as a list of spectra
    cls_list = list(cls)
    maxsize = max([len(c) for c in cls])

    if lmax is None or lmax < 0:
        lmax = maxsize-1
    if mmax is None or mmax < 0:
        mmax = lmax

    Nspec = sphtlib._getn(len(cls_list))

    if Nspec <= 0:
        if len(cls_list) == 4:
            if new: ## new input order: TT EE BB TE -> TT EE BB TE 0 0
                cls_list = [cls[0], cls[1], cls[2], cls[3], None, None]
            else: ## old input order: TT TE EE BB -> TT TE 0 EE 0 BB
                cls_list = [cls[0], cls[1], None, cls[2], None, cls[3]]
            Nspec = 3
        else:
            raise TypeError("The sequence of arrays must have either 4 elements "
                            "or n(n+1)/2 elements (some may be None)")

    szalm = hp.Alm.getsize(lmax,mmax)
    alms_list = []
    for i in range(Nspec):
        alm = np.zeros(szalm,'D')
        alm.real = prng.standard_normal(szalm)
        alm.imag = prng.standard_normal(szalm)
        alms_list.append(alm)
    if new: # new input order: input given by diagonal, should be given by row
        cls_list = hp.sphtfunc.new_to_old_spectra_order(cls_list)
    # ensure cls are float64
    cls_list = [(np.asarray(cl, dtype=np.float64) if cl is not None else None)
                for cl in cls_list]
    sphtlib._synalm(cls_list, alms_list, lmax, mmax)
    return alms_list

@core.usefulfunc
def synfast(cls, nside, lmax=None, mmax=None, alm=False, pol=True,
            pixwin=False, fwhm=0.0, sigma=None, beam=None, new=True,
            verbose=True, deriv=None, seed=None, executable=None,
            mapfile=None, almfile=None, parfile=None,
            talm=None, prng=None):
    """Create a map(s) from cl(s).  Same as healpy function, but
    interfaces with f90 synfast under-the-hood to handle derivatives,
    handles arbitrary beam window functions, and realizations with
    constrained temperature modes.

    Parameters
    ----------
    cls : array or tuple of array
      A cl or a list of cl (either 4 or 6, see :func:`synalm`)
    nside : int, scalar
      The nside of the output map(s)
    lmax : int, scalar, optional
      Maximum l for alm. Default: min of 3*nside-1 or length of the cls - 1
    mmax : int, scalar, optional
      Maximum m for alm.
    alm : bool, scalar, optional
      If True, return also alm(s). Default: False.
    pol : bool, optional
      If True, assumes input cls are TEB and correlation. Output will be TQU maps.
      (input must be 1, 4 or 6 cl's)
      If False, fields are assumed to be described by spin 0 spherical harmonics.
      (input can be any number of cl's)
      If there is only one input cl, it has no effect. Default: True.
    pixwin : bool, scalar, optional
      If True, convolve the alm by the pixel window function. Default: False.
    fwhm : float, scalar, optional
      The fwhm of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    sigma : float, scalar, optional
      The sigma of the Gaussian used to smooth the map (applied on alm)
      [in radians]
    beam : array or sequence of 3 arrays, optional
      If supplied, the beam function is applied instead of a Gaussian
      beam to each alm.
    new : bool, optional, default True
      order of input cls.  SEE NOTES!
    deriv : int, scalar, optional
      If 1, calculate first derivative maps
      If 2, calculate first and second derivative maps
    seed : int or array, optional
      If set, reseed the random number generator
    executable : string, optional
      If set, override path to the synfast f90 binary.  Default: synfast
    mapfile : string, optional
      If set, store map to this destination.  Otherwise, the map will be
      written to a temporary location by the synfast executable.
    almfile : string, optional
      If set, store alm to this destination for passing to the synfast
      executable.  Otherwise, the alms will be stored and read from a
      temporary location.
    parfile : string, optional
      If set, store the parameter file to this destination for passing to
      the synfast executable.  Otherwise, the parameter file will be stored
      and read from a temporary location.
    talm : array, optional
      Input temperature alms.  If supplied, the E mode alms
      are constrained by these temperature modes.  Size must match the requested
      lmax and mmax.
    verbose : bool, optional
      print stuff or not.
    prng : numpy.random.RandomState instance, optional
      If supplied, use this RandomState instance to generate maps.
      Overrides the `seed` option.

    Returns
    -------
    maps : array or tuple of arrays
      The output map (possibly list of maps if polarized input).
      or, if alm is True, a tuple of (map,alm)
      (alm possibly a list of alm if polarized input)

    Notes
    -----
    The order of the spectra will change in a future release. The new= parameter
    help to make the transition smoother. You can start using the new order
    by setting new=True.
    In the next version of healpy, the default will be new=True.
    This change is done for consistency between the different tools
    (alm2cl, synfast, anafast).
    In the new order, the spectra are ordered by diagonal of the correlation
    matrix. Eg, if fields are T, E, B, the spectra are TT, EE, BB, TE, EB, TB
    with new=True, and TT, TE, TB, EE, EB, BB if new=False.

    Intermediate files not explicitly stored elsewhere will be written to
    a temporary directory, created using the tempfile.mkdtemp() utility.
    The temporary directory will be cleared upon completion.
    """

    if executable is None:
        executable = synfast_exe

    # temporary directory
    if deriv:
        tmproot = tf.mkdtemp()

    try:
        if deriv and deriv not in [1,2]:
            raise ValueError("Wrong deriv value (must be 1 or 2 if set)")

        if not hp.isnsideok(nside):
            raise ValueError("Wrong nside value (must be a power of two).")

        if lmax is None or lmax < 0:
            cls_lmax = hp.cookbook.len_array_or_arrays(cls) - 1
            lmax = min(cls_lmax, 3 * nside - 1)
        else:
            lmax = lmax

        if deriv:
            # create and populate parameter file
            if parfile is None:
                parfile = os.path.join(tmproot, 'synfast.par')
            f = open(parfile, 'w')

            simul_type = deriv + (1 + (deriv > 0)) * (1 + (pol is True))
            f.write('simul_type = %d\n' % simul_type)

            if almfile is None:
                almfile = os.path.join(tmproot, 'alm_in.fits')
            f.write('almsfile = %s\n' % almfile)

            f.write('nsmax = %d\n' % nside)
            f.write('nlmax = %d\n' % lmax)

            # be explicit
            f.write("infile = ''\n")
            f.write('apply_windows = 0\n')
            f.write("windowfile = ''\n")

            if mapfile is None:
                mapfile = os.path.join(tmproot, 'map.fits')
            f.write('outfile = %s\n' % mapfile)
            mapfile = mapfile.lstrip('!')

            f.close()

        # generate alms and apply windows
        alms = synalm(cls, lmax=lmax, mmax=mmax, new=new, verbose=verbose,
                      seed=seed, prng=prng)

        # constrained realization
        # following eq. 1 of bicep2 2014 paper
        # alme = clte / cltt * almt + sqrt(clee - clte**2/cltt) * n
        # so just subtract original almt term and replace with constrained
        if pol and talm is not None:
            # Catch case when spec is 1d but pol is True (i.e. source_pol=False)
            if np.asarray(cls).ndim == 2:
                # clte / cltt
                fac = (cls[3] if new else cls[1]).copy()
                msk = cls[0] != 0
                fac[~msk] = 0
                fac[msk] /= cls[0][msk]

                # replace almt
                # NB: size of talm must match here
                alms[1] += hp.almxfl(talm - alms[0], fac, mmax=mmax, 
                                     inplace=False)
                alms[0][:] = talm
            else:
                alms = talm

        smoothalm(alms, fwhm=fwhm, sigma=sigma, beam=beam, pol=pol,
                  mmax=mmax, verbose=verbose, inplace=True)

        # apply pixwin if necessary
        if pixwin is True:
            if not pol:
                alms = [alms]
            pw = hp.pixwin(nside, True)
            for ialm, alm in enumerate(alms):
                pixelwindow = pw[1] if ialm > 0 and pol else pw[0]
                hp.almxfl(alm, pixelwindow, mmax=mmax, inplace=True)
            if not pol:
                alms = alms[0]

        # create map and return if no derivatives required
        if not deriv:
            m = np.asarray(hp.alm2map(alms, nside, pixwin=False, lmax=lmax,
                                        mmax=mmax, verbose=verbose))
            if alm is True:
                return m, alms
            return m

        # store alms to pass to fortran
        write_alm(almfile, alms, lmax=lmax, mmax=mmax)

        # run synfast
        stdout = None if verbose else open(os.devnull, 'w')
        sp.check_call([executable, parfile], stdout=stdout, stderr=sp.STDOUT)
        if not verbose:
            stdout.close()

        # read in output
        n_pols = 1 + 2 * (pol is True)
        fmax = max(n_pols, 3 * deriv * n_pols)
        maps = read_map(mapfile, range(fmax), verbose=verbose)

        # return
        if alm is True:
            return maps, alms
        return maps

    finally:
        # cleanup
        if deriv:
            shutil.rmtree(tmproot)

@core.usefulfunc
def alteralm(alm_in, nside_in, nlmax=None, fwhm_in=None, fwhm_out=None,
             beam_in=None, beam_out=None, nside_out=None,
             coord_in=None, coord_out=None, epoch_in=None, epoch_out=None,
             double_precision=False, alm_out=None, verbose=True, **kwargs):
    """
    Modify a set of Alm spherical harmonics for use with anafast/synfast.

    Parameters
    ----------
    alm_in : array, tuple of array or string
        Input alms to be modified.  Can be numpy arrays or a path
        to a fits file.
    nside_in : int
        resolution parameter of the source map.  The pixel window function
        will be divided out if `nside_out` is set.
    nlmax : int, optional
        Maximum ell value for output alms.  Default: maximum input ell.
    fwhm_in, fwhm_out : float, optional
        Beam FWHM to divide out (fwhm_in) or multiply by (fwhm_out),
        in **radians**.
    beam_in, beam_out : array or string, optional
        Beam window function to divide out (beam_in) or multiply by (beam_out).
    nside_out : int, optional
        resolution parameter whose window function to multiply by
        Default: same as `nside_in`. Set to zero to only correct for input
        map pixel window function.
    coord_in, coord_out : 'C', 'E' or 'G'; optional
        Coordinate system to rotate from/to.
    epoch_in, epoch_out : float, optional
        Epoch of input/output coordinate system
    double_precision : bool, optional
        If True, use double precision floats, otherwise single precision.
        Default: False.
    alm_out : string, optional
        Filename of output alm fits file, if storing for further use.
    executable : string, optional
        Path to `alteralm` binary.
    parfile : string, optional
        Path to parameter file.

    Returns
    -------
    alm_out : array or tuple of arrays
    """

    executable = [kwargs.pop('executable', alteralm_exe)]
    if double_precision:
        executable += ['--double']
    tmproot = tf.mkdtemp()
    params = {}

    # process parameters
    if isinstance(alm_in, str):
        params['infile_alms'] = alm_in
    else:
        alm_file = os.path.join(tmproot, 'alm_in.fits')
        write_alm(alm_file, alm_in)
        params['infile_alms'] = alm_file

    if alm_out is None:
        alm_out = os.path.join(tmproot, 'alm_out.fits')
    params['outfile_alms'] = alm_out
    alm_out = alm_out.lstrip('!')

    if fwhm_in is not None:
        params['fwhm_arcmin_in'] = np.degrees(fwhm_in) * 60.
    if fwhm_out is not None:
        params['fwhm_arcmin_out'] = np.degrees(fwhm_out) * 60.

    if beam_in is not None:
        if isinstance(beam_in, str):
            params['beam_file_in'] = beam_in
        else:
            beam_file = os.path.join(tmproot, 'beam_in.fits')
            hp.write_cl(beam_in)
            params['beam_file_in'] = beam_file
    if beam_out is not None:
        if isinstance(beam_out, str):
            params['beam_file_out'] = beam_out
        else:
            beam_file = os.path.join(tmproot, 'beam_out.fits')
            hp.write_cl(beam_out)
            params['beam_file_out'] = beam_file

    if not hp.isnsideok(nside_in):
        raise ValueError('Invalid nside_in')
    params['nsmax_in'] = nside_in
    if nside_out is not None:
        if hp.isnsideok(nside_out) or nside_out == 0:
            params['nsmax_out'] = nside_out
        else:
            raise ValueError('Invalid nside_out')

    if coord_in is not None:
        if coord_in not in 'CEG':
            raise ValueError('Invalid coord_in')
        params['coord_in'] = coord_in
    if coord_out is not None:
        if coord_out not in 'CEG':
            raise ValueError('Invalid coord_out')
        params['coord_out'] = coord_out

    if epoch_in is not None:
        params['epoch_in'] = epoch_in
    if epoch_out is not None:
        params['epoch_out'] = epoch_out

    # populate parameter file
    parfile = kwargs.pop('parfile', os.path.join(tmproot, 'alteralm.par'))
    f = open(parfile, 'w')
    for k, v in params.items():
        if v is not None:
            f.write('{} = {}\n'.format(k, v))
    f.close()

    # run alteralm
    stdout = None if verbose else open(os.devnull, 'w')
    try:
        sp.check_call(executable + [parfile], stdout=stdout, stderr=sp.STDOUT)
    except (sp.CalledProcessError, KeyboardInterrupt):
        shutil.rmtree(tmproot)
        raise
    if stdout:
        stdout.close()

    # read in output
    from healpy.fitsfunc import pf
    hdulist = pf.open(alm_out)
    npol = len(hdulist)-1
    hdulist.close()
    alm_out = tuple(hp.read_alm(alm_out, hdu=hdu) for hdu in range(1, npol+1))
    if npol == 1:
        alm_out = alm_out[0]

    # cleanup
    shutil.rmtree(tmproot)

    # return
    return alm_out

@core.usefulfunc
def smoothalm(alms, fwhm=0.0, sigma=None, beam=None, pol=True,
              mmax=None, verbose=True, inplace=True):
    """Smooth alm with a Gaussian symmetric beam or custom window function.

    Parameters
    ----------
    alms : array or sequence of 3 arrays
      Either an array representing one alm, or a sequence of arrays.
      See *pol* parameter.
    fwhm : float, optional
      The full width half max parameter of the Gaussian. Default:0.0
      [in radians]
    sigma : float, optional
      The sigma of the Gaussian. Override fwhm.
      [in radians]
    beam : array or sequence of 3 arrays, optional
      If supplied, the beam function is applied instead of a Gaussian
      beam to each alm.
    pol : bool, optional
      If True, assumes input alms are TEB. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, apply spin 0 harmonic transform to each alm.
      (input can be any number of alms)
      If there is only one input alm, it has no effect. Default: True.
    mmax : None or int, optional
      The maximum m for alm. Default: mmax=lmax
    inplace : bool, optional
      If True, the alm's are modified inplace if they are contiguous arrays
      of type complex128. Otherwise, a copy of alm is made. Default: True.
    verbose : bool, optional
      If True prints diagnostic information. Default: True

    Returns
    -------
    alms : array or sequence of 3 arrays
      The smoothed alm. If alm[i] is a contiguous array of type complex128,
      and *inplace* is True the smoothing is applied inplace.
      Otherwise, a copy is made.
    """

    # make imports identical to healpy source for easy porting
    from healpy.sphtfunc import almxfl, Alm
    from healpy import cookbook as cb

    if beam is None:
        if sigma is None:
            sigma = fwhm / (2.*np.sqrt(2.*np.log(2.)))

        if verbose:
            print("Sigma is {0:f} arcmin ({1:f} rad) ".format(sigma*60*180/np.pi,sigma))
            print("-> fwhm is {0:f} arcmin".format(sigma*60*180/np.pi*(2.*np.sqrt(2.*np.log(2.)))))

    # Check alms
    if not cb.is_seq(alms):
        raise ValueError("alm must be a sequence")

    if sigma == 0 and beam is None:
        # nothing to be done
        return alms

    lonely = False
    if not cb.is_seq_of_seq(alms):
        alms = [alms]
        lonely = True

    # check beam
    if beam is not None:
        if not cb.is_seq(beam):
            raise ValueError("beam must be a sequence")
        if not lonely:
            if not cb.is_seq_of_seq(beam):
                beam = [beam]*len(alms)
            else:
                if len(beam) != len(alms):
                    raise ValueError("alm and beam shape mismatch")
        else:
            if cb.is_seq_of_seq(beam):
                raise ValueError("alm and beam shape mismatch")
            else:
                beam = [beam]

    # we have 3 alms -> apply smoothing to each map.
    # polarization has different B_l from temperature
    # exp{-[ell(ell+1) - s**2] * sigma**2/2}
    # with s the spin of spherical harmonics
    # s = 2 for pol, s=0 for temperature
    retalm = []
    for ialm, alm in enumerate(alms):
        lmax = Alm.getlmax(len(alm), mmax)
        if lmax < 0:
            raise TypeError('Wrong alm size for the given '
                            'mmax (len(alms[%d]) = %d).'%(ialm, len(alm)))
        if beam is None:
            ell = np.arange(lmax + 1.)
            s = 2 if ialm >= 1 and pol else 0
            fact = np.exp(-0.5 * (ell * (ell + 1) - s ** 2) * sigma ** 2)
        else:
            fact = beam[ialm]
        res = almxfl(alm, fact, mmax=mmax, inplace=inplace)
        retalm.append(res)
    # Test what to return (inplace/not inplace...)
    # Case 1: 1d input, return 1d output
    if lonely:
        return retalm[0]
    # case 2: 2d input, check if in-place smoothing for all alm's
    for i in range(len(alms)):
        samearray = alms[i] is retalm[i]
        if not samearray:
            # Case 2a:
            # at least one of the alm could not be smoothed in place:
            # return the list of alm
            return retalm
    # Case 2b:
    # all smoothing have been performed in place:
    # return the input alms
    return alms

@core.usefulfunc
def smoothing(map_in, fwhm=0.0, sigma=None, beam=None, pol=True,
              iter=3, lmax=None, mmax=None, use_weights=False,
              fill=np.nan, datapath=None, verbose=True):
    """Smooth a map with a Gaussian symmetric beam or custom window function.

    No removal of monopole or dipole is performed.

    Parameters
    ----------
    map_in : array or sequence of 3 arrays
      Either an array representing one map, or a sequence of
      3 arrays representing 3 maps, accepts masked arrays
    fwhm : float, optional
      The full width half max parameter of the Gaussian [in radians].
      Default:0.0
    sigma : float, optional
      The sigma of the Gaussian [in radians]. Override fwhm.
    beam : array or sequence of 3 arrays, optional
      If supplied, the beam function is applied instead of a Gaussian
      beam to each alm.
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TQU maps.
      (input must be 1 or 3 alms)
      If False, each map is assumed to be a spin 0 map and is
      treated independently (input can be any number of alms).
      If there is only one input map, it has no effect. Default: True.
    iter : int, scalar, optional
      Number of iteration (default: 3)
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    fill : scalar, optional
      Fill the bad pixels with this value, if supplied.  Default: NaN.
    datapath : None or str, optional
      If given, the directory where to find the weights data.
    verbose : bool, optional
      If True prints diagnostic information. Default: True

    Returns
    -------
    maps : array or list of 3 arrays
      The smoothed map(s)
    """

    # make imports identical to healpy source for easy porting
    from healpy import pixelfunc
    from healpy.sphtfunc import map2alm, alm2map
    from healpy import cookbook as cb

    if not cb.is_seq(map_in):
        raise TypeError("map_in must be a sequence")

    # save the masks of inputs
    masks = mask_bad(map_in, badnan=True, badinf=True)
    if np.any(masks):
        map_in = np.array(map_in, copy=True)
        map_in[masks] = hp.UNSEEN

    if cb.is_seq_of_seq(map_in):
        nside = pixelfunc.npix2nside(len(map_in[0]))
        n_maps = len(map_in)
    else:
        nside = pixelfunc.npix2nside(len(map_in))
        n_maps = 0

    if pol or n_maps in (0, 1):
        # Treat the maps together (1 or 3 maps)
        alms = map2alm(map_in, lmax=lmax, mmax=mmax, iter=iter,
                       pol=pol, use_weights=use_weights,
                       datapath=datapath)
        smoothalm(alms, fwhm=fwhm, sigma=sigma, beam=beam,
                  inplace=True, verbose=verbose)
        output_map = alm2map(alms, nside, pixwin=False, verbose=verbose)
    else:
        # Treat each map independently (any number)
        output_map = []
        for m in map_in:
            alm = map2alm(m, lmax=lmax, mmax=mmax, iter=iter, pol=pol,
                          use_weights=use_weights, datapath=datapath)
            smoothalm(alm, fwhm=fwhm, sigma=sigma, beam=beam,
                      inplace=True, verbose=verbose)
            output_map.append(alm2map(alm, nside, pixwin=False, verbose=verbose))

    output_map = np.asarray(output_map)
    output_map[masks] = fill
    return output_map

@core.usefulfunc
def anafast(map1, map2=None, nspec=None, lmax=None, mmax=None,
            iter=3, alm=False, pol=True, use_weights=False,
            datapath=None):
    """Computes the power spectrum of an Healpix map, or the cross-spectrum
    between two maps if *map2* is given.
    No removal of monopole or dipole is performed.
    NaN and inf values are masked before computing.

    Parameters
    ----------
    map1 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps
    map2 : float, array-like shape (Npix,) or (3, Npix)
      Either an array representing a map, or a sequence of 3 arrays
      representing I, Q, U maps
    nspec : None or int, optional
      The number of spectra to return. If None, returns all, otherwise
      returns cls[:nspec]
    lmax : int, scalar, optional
      Maximum l of the power spectrum (default: 3*nside-1)
    mmax : int, scalar, optional
      Maximum m of the alm (default: lmax)
    iter : int, scalar, optional
      Number of iteration (default: 3)
    alm : bool, scalar, optional
      If True, returns both cl and alm, otherwise only cl is returned
    pol : bool, optional
      If True, assumes input maps are TQU. Output will be TEB cl's and
      correlations (input must be 1 or 3 maps).
      If False, maps are assumed to be described by spin 0 spherical harmonics.
      (input can be any number of maps)
      If there is only one input map, it has no effect. Default: True.
    use_weights: bool, scalar, optional
      If True, use the ring weighting. Default: False.
    datapath : None or str, optional
      If given, the directory where to find the weights data.

    Returns
    -------
    res : array or sequence of arrays
      If *alm* is False, returns cl or a list of cl's (TT, EE, BB, TE, EB, TB for
      polarized input map)
      Otherwise, returns a tuple (cl, alm), where cl is as above and
      alm is the spherical harmonic transform or a list of almT, almE, almB
      for polarized input
    """

    # mask bad values
    if not hp.pixelfunc.is_ma(map1):
        if np.any(np.isnan(map1) | np.isinf(map1)):
            map1 = np.array(map1, copy=True)
            map1[mask_bad(map1, badnan=True, badinf=True)] = hp.UNSEEN
    if map2 is not None and not hp.pixelfunc.is_ma(map2):
        if np.any(np.isnan(map2) | np.isinf(map2)):
            map2 = np.array(map2, copy=True)
            map2[mask_bad(map2, badnan=True, badinf=True)] = hp.UNSEEN

    # run anafast
    ret = hp.anafast(map1, map2=map2, nspec=nspec, lmax=lmax, mmax=mmax,
                      iter=iter, alm=alm, pol=pol, use_weights=use_weights,
                      datapath=datapath)

    # return
    return ret

@core.usefulfunc
def write_map(filename, m, nest=False, dtype=np.float64, fits_IDL=True,
              coord=None, partial=False, mask=None, pixels=None, nside=None,
              column_names=None, column_units=None, extra_header=(),
              append=False, return_hdu=False):
    """Writes an healpix map into an healpix file.

    Parameters
    ----------
    filename : str
      the fits file name
    m : array or sequence of 3 arrays
      the map to write. Possibly a sequence of 3 maps of same size.
      They will be considered as I, Q, U maps.
      Supports masked maps, see the `ma` function.
    nest : bool, optional
      If True, ordering scheme is assumed to be NESTED, otherwise, RING. Default: RING.
      The map ordering is not modified by this function, the input map array
      should already be in the desired ordering (run `ud_grade` beforehand).
    fits_IDL : bool, optional
      If True, reshapes columns in rows of 1024, otherwise all the data will
      go in one column. Default: True
    coord : str
      The coordinate system, typically 'E' for Ecliptic, 'G' for Galactic or 'C' for
      Celestial (equatorial)
    partial : bool, optional
      If True, fits file is written as a partial-sky file with explicit indexing.
      Otherwise, implicit indexing is used.  Default: False.
    mask : bool array, optional
      If supplied, mask (1=good, 0=bad) is applied to the input map, and the result
      is stored as a partial map.  Overrides `partial` option.
    pixels : index array, optional
      If supplied, the input map is assumed to be a partial map containing only
      these pixels.  Overrides `mask` and `partial` options.
    nside : int, optional
      If `pixels` is supplied, this argument is required to verify the map shape.
    column_names : str or list
      Column name or list of column names, if None we use:
      I_STOKES for 1 component,
      I/Q/U_STOKES for 3 components,
      II, IQ, IU, QQ, QU, UU for 6 components,
      COLUMN_0, COLUMN_1... otherwise
    column_units : str or list
      Units for each column, or same units for all columns.
    extra_header : list or dict
      Extra records to add to FITS header.
    dtype : np.dtype or list of np.dtypes, optional
      The datatype in which the columns will be stored. Will be converted
      internally from the numpy datatype to the fits convention. If a list,
      the length must correspond to the number of map arrays.
    append : bool
      Set this option to append the map to an existing file as a new HDU.
    return_hdu : bool
      Set this option to return the BinTableHDU that would be written, rather
      that writing it to disk.
    """
    from healpy.fitsfunc import getformat, standard_column_names, pf
    from healpy import pixelfunc

    standard_column_names.update({
            4: ['{}_STOKES'.format(comp) for comp in 'IQUV'],
            10: ['II', 'IQ', 'IU', 'IV', 'QQ', 'QU', 'QV', 'UU', 'UV', 'VV']
            })

    if not hasattr(m, '__len__'):
        raise TypeError('The map must be a sequence')

    m = pixelfunc.ma_to_array(m)
    if pixels is None:
        if pixelfunc.maptype(m) == 0: # a single map is converted to a list
            m = [m]
    else:
        m = np.atleast_2d(m)

    # check the dtype and convert it
    try:
        fitsformat = []
        for curr_dtype in dtype:
            fitsformat.append(getformat(curr_dtype))
    except TypeError:
        #dtype is not iterable
        fitsformat = [getformat(dtype)] * len(m)

    if column_names is None:
        column_names = standard_column_names.get(
            len(m), ["COLUMN_%d" % n for n in range(len(m))])
    else:
        assert len(column_names) == len(m), \
            "Length column_names != number of maps"

    if column_units is None or isinstance(column_units, str):
        column_units = [column_units] * len(m)

    # maps must have same length
    assert len(set(map(len, m))) == 1, "Maps must have same length"
    if pixels is None:
        nside = pixelfunc.npix2nside(len(m[0]))
    elif nside is None:
        raise ValueError('Invalid healpix map : nside required')

    if nside < 0:
        raise ValueError('Invalid healpix map : wrong number of pixel')

    cols=[]

    if mask is not None or pixels is not None:
        partial = True

    if partial:
        fits_IDL = False
        if pixels is not None:
            pix = pixels
            if any([mm.shape != pix.shape for mm in m]):
                raise ValueError('Invalid healpix map : pixel index mismatch')
        else:
            if mask is None:
                mask = mask_good(m[0])
            m = [mm[mask] for mm in m]
            pix = np.where(mask)[0]
        if len(pix) == 0:
            raise ValueError('Invalid healpix map : empty partial map')
        ff = getformat(np.min_scalar_type(-pix.max()))
        if ff is None:
            ff = 'I'
        cols.append(pf.Column(name='PIXEL',
                              format=ff,
                              array=pix,
                              unit=None))

    for cn, cu, mm, curr_fitsformat in zip(column_names, column_units, m,
                                           fitsformat):
        if len(mm) > 1024 and fits_IDL:
            # I need an ndarray, for reshape:
            mm2 = np.asarray(mm)
            cols.append(pf.Column(name=cn,
                                  format='1024%s' % curr_fitsformat,
                                  array=mm2.reshape(mm2.size/1024,1024),
                                  unit=cu))
        else:
            cols.append(pf.Column(name=cn,
                                  format='%s' % curr_fitsformat,
                                  array=mm,
                                  unit=cu))

    tbhdu = pf.BinTableHDU.from_columns(cols)
    # add needed keywords
    tbhdu.header.set('PIXTYPE', 'HEALPIX', 'HEALPIX pixelisation')
    tbhdu.header.set('ORDERING', 'NESTED' if nest else 'RING',
                     'Pixel ordering scheme, either RING or NESTED')
    if coord:
        tbhdu.header.set('COORDSYS', coord,
                         'Ecliptic, Galactic or Celestial (equatorial)')
    tbhdu.header.set('EXTNAME', 'xtension',
                     'name of this binary table extension')
    tbhdu.header.set('NSIDE', nside, 'Resolution parameter of HEALPIX')
    if not partial:
        tbhdu.header.set('FIRSTPIX', 0, 'First pixel # (0 based)')
        tbhdu.header.set('LASTPIX', pixelfunc.nside2npix(nside) - 1,
                         'Last pixel # (0 based)')
    tbhdu.header.set('INDXSCHM', 'EXPLICIT' if partial else 'IMPLICIT',
                     'Indexing: IMPLICIT or EXPLICIT')
    tbhdu.header.set('OBJECT', 'PARTIAL' if partial else 'FULLSKY',
                     'Sky coverage, either FULLSKY or PARTIAL')

    if not isinstance(extra_header, dict):
        for args in extra_header:
            if args[0] == 'COMMENT':
                tbhdu.header.add_comment(*args[1:])
            else:
                tbhdu.header.set(*args)
    else:
        tbhdu.header.update(extra_header)

    if return_hdu:
        return tbhdu

    if not append:
        from astropy import __version__ as astropy_version
        if astropy_version >= "1.3":
            tbhdu.writeto(filename, overwrite=True)
        else:
            tbhdu.writeto(filename, clobber=True)
    else:
        if isinstance(filename, str):
            if not os.path.exists(filename):
                # doesn't exist yet. write normally, with dummy Primary HDU
                tbhdu.writeto(filename)
        pf.append(filename, tbhdu.data, tbhdu.header, verify=False)

@core.usefulfunc
def read_map(filename, field=0, dtype=np.float64, nest=False, partial=False,
             fill=np.nan, hdu=1, h=False, verbose=False, memmap=False,
             return_part=False, return_pix=True, return_names=False):
    """Read an healpix map from a fits file.  Partial sky files are expanded
    to full size and filled with UNSEEN.

    Parameters
    ----------
    filename : str
      The fits file name.
      Can also be an HDUList object from astropy.io.fits.open or a
      particular HDU from the list.
    field : int or tuple of int, or None, optional
      The column to read. Default: 0.
      By convention 0 is temperature, 1 is Q, 2 is U.
      Field can be a tuple to read multiple columns (0,1,2)
      If the fits file is a partial-sky file, field=0 corresponds to the
      first column after the pixel index column.
      If None, all columns are read in.
    dtype : data type or list of data types, optional
      Force the conversion to some type. Passing a list allows different
      types for each field. In that case, the length of the list must
      correspond to the length of the field parameter. Default: np.float64
    nest : bool, optional
      If True return the map in NEST ordering, otherwise in RING ordering;
      use fits keyword ORDERING to decide whether conversion is needed or not
      If None, no conversion is performed.
    partial : bool, optional
      If True, fits file is assumed to be a partial-sky file with explicit indexing,
      if the indexing scheme cannot be determined from the header.
      If False, implicit indexing is assumed.  Default: False.
      A partial sky file is one in which OBJECT=PARTIAL and INDXSCHM=EXPLICIT,
      and the first column is then assumed to contain pixel indices.
      A full sky file is one in which OBJECT=FULLSKY and INDXSCHM=IMPLICIT.
      At least one of these keywords must be set for the indexing
      scheme to be properly identified.
    return_part : bool, optional
      If the map is a partial-sky file (see 'partial' above), don't fill
      out to full-sky. Return the map and the pixels array (if `return_pix` is True).
    return_pix : bool, optional
      If the map is a partial-sky file (see 'partial' above), and `return_part`
      is True, return the pixel array.
    return_names : bool, optional
      If True, return the names of fields that have been read.
    fill : scalar, optional
      Fill the bad pixels with this value, if supplied.  Default: NaN.
    hdu : int, optional
      the header number to look at (start at 0)
    h : bool, optional
      If True, return also the header. Default: False.
    verbose : bool, optional
      If True, print a number of diagnostic messages
    memmap : bool, optional
      Argument passed to astropy.io.fits.open, if True, the map is not read into memory,
      but only the required pixels are read when needed. Default: False.

    Returns
    -------
    m | (m0, m1, ...) [, header] : 1D or 2D array
      The map(s) read from the file
    pix : (If return_part is True and return_pix is True) 1D array
      List of pixels contained in partial-sky map. If return_part=True but
      the map is not partial-sky, this will be None
    nside : (If return_part is True) int
      healpix nside of the map. Needed for partial-sky
    header : (if h is True)
      The FITS header
    """
    import warnings
    from healpy.fitsfunc import pf, HealpixFitsWarning
    from healpy import pixelfunc, UNSEEN

    if isinstance(filename, str):
        hdulist = pf.open(filename, memmap=memmap)
        fits_hdu = hdulist[hdu]
    elif isinstance(filename, pf.HDUList):
        fits_hdu = filename[hdu]
    else:
        # assume it's an HDU directly
        fits_hdu = filename

    if not isinstance(fits_hdu, pf.BinTableHDU):
        raise TypeError("FITS error: Healpix map must be a binary table")

    # check nside
    nside = fits_hdu.header.get('NSIDE')
    if nside is None:
        warnings.warn(
            "No NSIDE in the header file : will use length of array",
            HealpixFitsWarning)
    else:
        nside = int(nside)
    if verbose:
        print('NSIDE = {0:d}'.format(nside))
    if not pixelfunc.isnsideok(nside):
        raise ValueError('Wrong nside parameter.')
    sz = pixelfunc.nside2npix(nside)

    # check ordering
    ordering = fits_hdu.header.get('ORDERING', 'UNDEF').strip()
    if ordering == 'UNDEF':
        ordering = (nest and 'NESTED' or 'RING')
        warnings.warn("No ORDERING keyword in header file : "
                      "assume {}".format(ordering))
    if verbose:
        print('ORDERING = {0:s} in fits file'.format(ordering))

    # partial sky: check OBJECT, then INDXSCHM
    obj = fits_hdu.header.get('OBJECT', 'UNDEF').strip()
    if obj != 'UNDEF':
        if obj == 'PARTIAL':
            partial = True
        elif obj == 'FULLSKY':
            partial = False

    schm = fits_hdu.header.get('INDXSCHM', 'UNDEF').strip()
    if schm != 'UNDEF':
        if schm == 'EXPLICIT':
            if obj == 'FULLSKY':
                raise ValueError('Incompatible INDXSCHM keyword')
            partial = True
        elif schm == 'IMPLICIT':
            if obj == 'PARTIAL':
                raise ValueError('Incompatible INDXSCHM keyword')
            partial = False

    if schm == 'UNDEF':
        schm = 'EXPLICIT' if partial else 'IMPLICIT'
        #warnings.warn("No INDXSCHM keyword in header file : "
                       #"assume {}".format(schm))
    if verbose:
        print('INDXSCHM = {0:s}'.format(schm))

    # check field
    if field is None:
        field = range(len(fits_hdu.data.columns) - 1*partial)
    if not (hasattr(field, '__len__') or isinstance(field, str)):
        field = (field,)
    ret = []

    if return_names:
        names = fits_hdu.data.names

    if not return_part:
        return_pix = False

    if partial:
        # increment field counters
        field = tuple(f if isinstance(f, str) else f+1
                      for f in field)
        if return_pix or not return_part:
            try:
                pix = fits_hdu.data.field(0).astype(int).ravel()
            except pf.VerifyError as e:
                print(e)
                print("Trying to fix a badly formatted header")
                fits_hdu.verify("fix")
                pix = fits_hdu.data.field(0).astype(int).ravel()
        else:
            pix = None
    else:
        pix = None

    if return_names:
        rnames = [f if isinstance(f, str) else names[f] for f in field]
    else:
        rnames = None

    try:
        assert len(dtype) == len(field), \
            "The number of dtypes are not equal to the number of fields"
    except TypeError:
        dtype = [dtype] * len(field)

    for ff, curr_dtype in zip(field, dtype):
        try:
            m = fits_hdu.data.field(ff).astype(curr_dtype).ravel()
        except pf.VerifyError as e:
            print(e)
            print("Trying to fix a badly formatted header")
            fits_hdu.verify("fix")
            m = fits_hdu.data.field(ff).astype(curr_dtype).ravel()

        if partial and not return_part:
            mnew = fill * np.ones(sz, dtype=curr_dtype)
            mnew[pix] = m
            m = mnew

        if (not pixelfunc.isnpixok(m.size) or \
            (sz>0 and sz != m.size)) and verbose:
            print('nside={0:d}, sz={1:d}, m.size={2:d}'.format(nside, sz, m.size))
            raise ValueError('Wrong nside parameter.')

        if not nest is None: # no conversion with None
            if nest and ordering == 'RING':
                idx = pixelfunc.nest2ring(
                    nside, np.arange(m.size, dtype=np.int32))
                m = m[idx]
                if verbose:
                    print('Ordering converted to NEST')
            elif (not nest) and ordering == 'NESTED':
                idx = pixelfunc.ring2nest(
                    nside, np.arange(m.size, dtype=np.int32))
                m = m[idx]
                if verbose:
                    print('Ordering converted to RING')

        try:
            m[mask_bad(m)] = fill
        except OverflowError:
            pass
        ret.append(m)

    # convert list of map arrays to 1D or 2D
    if len(ret) == 1:
        ret = ret[0]
    else:
        ret = np.asarray(ret)
    # append pixel array or FITS header as requested
    ret = ((ret,) + ((pix,) * return_pix + (nside,)) * return_part +
           (fits_hdu.header.items(),) * h + (rnames,) * return_names)
    if len(ret) == 1:
        return ret[0]
    return ret

@core.usefulfunc
def hotspot(m, return_values=False, extrema_outfile=None,
            maxima_outfile=None, minima_outfile=None, verbose=False,
            **kwargs):
    """
    Python interface to HEALPIX `hotspot` Fortran utility.

    Arguments
    ---------
    m : 1D array or map file
        Map in which to search for peaks
    return_values : bool, optional
        If True, return arrays of map values at each extremum along with
        the position arrays.
    extrema_outfile : string, optional
        If supplied, store map of extrema to this fits file.
        Map is nonzero for pixels at minima or maxima in the input map.
    maxima_outfile, minima_outfile : string, optional
        If supplied, store positions and values of minima/maxima in
        these data files.
    verbose : bool, optional
        If True, print output of the fortran utility to STDOUT.

    Returns
    -------
    minpos, maxpos : 1D arrays
        Pixel numbers of minima/maxima in the input map
    minval, maxval : 1D arrays, optional
        Map values at minima/maxima in the input map
    """

    executable = [kwargs.pop('executable', hotspot_exe)]
    tmproot = tf.mkdtemp()
    params = {}

    if isinstance(m, str):
        params['infile'] = m
    else:
        if m.ndim > 1:
            raise ValueError('Input map must be 1D')
        infile = os.path.join(tmproot, 'map.fits')
        write_map(infile, m)
        params['infile'] = infile

    if extrema_outfile is None:
        extrema_outfile = os.path.join(tmproot, 'pixlminmax.fits')
    params['extrema_outfile'] = extrema_outfile

    if maxima_outfile is None:
        maxima_outfile = os.path.join(tmproot, 'maxima.dat')
    params['maxima_outfile'] = maxima_outfile

    if minima_outfile is None:
        minima_outfile = os.path.join(tmproot, 'minima.dat')
    params['minima_outfile'] = minima_outfile

    parfile = kwargs.pop('parfile', os.path.join(tmproot, 'hotspot.par'))
    with open(parfile, 'w') as f:
        for k, v in params.items():
            if v is not None:
                f.write('{} = {}\n'.format(k, v))

    stdout = None if verbose else open(os.devnull, 'w')
    try:
        sp.check_call(executable + [parfile], stdout=stdout, stderr=sp.STDOUT)
    except (sp.CalledProcessError, KeyboardInterrupt):
        shutil.rmtree(tmproot)
        raise
    if stdout:
        stdout.close()

    maxpos, maxval = np.loadtxt(maxima_outfile, unpack=True)
    minpos, minval = np.loadtxt(minima_outfile, unpack=True)
    maxpos = maxpos.astype(int)
    minpos = minpos.astype(int)

    if return_values:
        return minpos, maxpos, minval, maxval
    return minpos, maxpos

@core.usefulfunc
def stack_map(m, pix, size=5, radial=False, return_norm=False):
    """
    Stack the data at the given map pixels into a map of the given radius about
    the center.  This is not particularly fast.

    Arguments
    ---------
    m : array_like
        Input map, T or TQU.
    pix : array_like
        Map positions to stack.
    size : float, optional
        Stacked map radius, in degrees.
    radial : bool, optional
        If True, rotate the Q/U maps into radial Q/U coordinates.
    return_norm : bool, optional
        If True, return the stack normalization along with the stack data

    Returns
    -------
    stack : array_like
        Stacked map of the same shape as `m`.  The stack is centered at
        `ra, dec = (0, 0)`, and may be plotted with:

        >>> sa.map.cartview(stack, lonra=[-5, 5], latra=[-5, 5])
    norm : scalar, optional
        Stack normalization.  Use for coaddition.
    """

    m = np.atleast_2d(m)
    nside = hp.get_nside(m)
    stack = np.full_like(m, np.nan)

    # pixel vectors
    theta, phi = hp.pix2ang(nside, pix)
    dec = 90. - np.degrees(theta)
    ra = np.degrees(phi)

    # pixels within radius about ra=dec=0
    disc = hp.query_disc(nside, [1., 0., 0.], np.radians(size))
    stack[..., disc] = 0
    norm = 0

    # stack
    for ii, (r, d) in enumerate(zip(ra, dec)):
        # rotate to common pixel at ra=dec=0
        mr = rotate_map(m, rot=[r, d, 0], pixels=disc)[..., disc]
        # ignore any points that are not entirely within the mask
        if np.any(mask_bad(mr[0])):
            continue
        # remove local monopole
        mr -= np.mean(mr, axis=-1)[:, None]
        # add to stack
        stack[..., disc] += mr
        norm += 1

    # normalize and mask to radius
    stack[..., disc] /= norm

    # rotate QU to radial coordinates (flat-sky approximation)
    if radial and stack.shape[0] == 3:
        rotate_qu_stack(stack, inplace=True)

    # return
    if return_norm:
        return stack.squeeze(), norm
    return stack.squeeze()

@core.usefulfunc
def rotate_map(m, coord=['C', 'G'], rot=None, mask=None, pixels=None,
               pol_axis=[0.,0.,1.]):
    """
    Rotate an input map from one coordinate system to another or to place a
    particular point at centre in rotated map. This does the proper Q and U
    Stokes rotation. Sign of Q U rotation should be correct for inverse
    rotation back to original coords (psi -> -psi)

    e.g. m = rotate_map(m, rot=[phi,90.-theta,0.])

    takes point at original theta, phi to new coord ra=dec=0

    Arguments
    ---------
    m : array_like
        A single map or two (Q,U) or three (I,Q,U) maps
    coord : list of two coordinates, optional.
        Coordinates to rotate between.  Default: ['C', 'G']
    rot : scalar or sequence, optional
        Describe the rotation to apply.
        In the form (lon, lat, psi) (unit: degrees) : the point at
        longitude lon and latitude lat will be at the center of the rotated
        map. An additional rotation of angle psi around this direction is applied
    mask : 1D array
        If supplied, only pixels in the *rotated map* that fall within the mask
        are handled.
    pixels : 1D array
        If supplied, only pixels in the *rotated map* that are also in this list
        are handled. Overrides `mask`.
    pol_axis : 3-vector, optional
        Axis normal to the plane in which the Q/U coordinates are defined.

    Returns
    -------
    mr : array_like
        The rotated map.
    """

    m = np.atleast_2d(m)

    if m.ndim > 2:
        raise ValueError(
            'Input map array must have no more than two dimensions')
    if m.shape[0] not in [1,2,3]:
        raise ValueError(
            'Input map must have 1 (T only), 2 (Q/U only) or 3 (T/Q/U) columns')

    res = hp.get_nside(m)
    pol = m.shape[0] in [2,3]
    pol_only = m.shape[0] == 2

    if rot is None:
        #use default coord transform C->G
        R = hp.Rotator(coord=coord, inv=True)
    else:
        R = hp.Rotator(rot=rot, inv=True)

    try:
        # qpoint is a 50% speedup
        import qpoint as qp
        Q = qp.QPoint()
        use_qpoint = True
    except ImportError:
        use_qpoint = False

    if pixels is None:
        if mask is not None:
            pixels, = np.where(mask)
        else:
            pixels = np.arange(len(m[0]))

    # rotate new coordinate system to original coordinates
    theta, phi = hp.pix2ang(res, pixels)
    mtheta, mphi = R(theta, phi)
    del theta
    del phi

    mr = np.full_like(m, np.nan)

    if use_qpoint:
        ra = np.degrees(mphi)
        dec = 90. - np.degrees(mtheta)
        del mtheta
        del mphi
        if not pol_only:
            mr[0, pixels] = Q.get_interp_val(m[0], ra, dec)
    else:
        if not pol_only:
            mr[0, pixels] = hp.get_interp_val(m[0], mtheta, mphi)

    if not pol:
        return mr.squeeze()

    #interpolate Q and U (better before or after rot?)
    if use_qpoint:
        mr[-2, pixels] = Q.get_interp_val(m[-2], ra, dec)
        mr[-1, pixels] = Q.get_interp_val(m[-1], ra, dec)
    else:
        mr[-2, pixels] = hp.get_interp_val(m[-2], mtheta, mphi)
        mr[-1, pixels] = hp.get_interp_val(m[-1], mtheta, mphi)

    rotate_qu_radial(
        mr, coord=coord, rot=rot, pol_axis=pol_axis, inplace=True)

    return mr

@core.usefulfunc
def rotate_proj(p, coord=['C', 'G'], rot=None, mask=None, pixels=None,
               pol_axis=[0.,0.,1.]):
    """
    Rotate an input proj map from one coordinate system to another or to place a
    particular point at centre in rotated map. This does the proper Q and U
    Stokes rotation. Sign of Q U rotation should be correct for inverse
    rotation back to original coords (psi -> -psi)

    e.g. p = rotate_proj(p, rot=[phi,90.-theta,0.])

    takes point at original theta, phi to new coord ra=dec=0

    Arguments
    ---------
    m : array_like
        A single hits map or six (II,IQ,IU,QQ,QU,UU) maps
    coord : list of two coordinates, optional.
        Coordinates to rotate between.  Default: ['C', 'G']
    rot : scalar or sequence, optional
        Describe the rotation to apply.
        In the form (lon, lat, psi) (unit: degrees) : the point at
        longitude lon and latitude lat will be at the center of the rotated
        map. An additional rotation of angle psi around this direction is applied
    mask : 1D array
        If supplied, only pixels in the *rotated proj* that fall within the mask
        are handled.
    pixels : 1D array
        If supplied, only pixels in the *rotated proj* that are also in this list
        are handled. Overrides `mask`.
    pol_axis : 3-vector, optional
        Axis normal to the plane in which the Q/U coordinates are defined.

    Returns
    -------
    pr : array_like
        The rotated proj.
    """

    p = np.atleast_2d(p)

    if p.ndim > 2:
        raise ValueError(
            'Input map array must have no more than two dimensions')
    if p.shape[0] not in [1,6]:
        raise ValueError(
            'Input proj must have 1 (T only), 6 (T/Q/U) columns')

    if p.shape[0] == 1:
        return rotate_map(p, coord=coord, rot=rot, mask=mask, pixels=pixels,
                          pol_axis=pol_axis)

    res = hp.get_nside(p)

    if rot is None:
        # use default coord transform C->G
        R = hp.Rotator(coord=coord, inv=True)
    else:
        R = hp.Rotator(rot=rot, inv=True)

    try:
        # qpoint is a 50% speedup
        import qpoint as qp
        Q = qp.QPoint()
        use_qpoint = True
    except ImportError:
        use_qpoint = False

    if pixels is None:
        if mask is not None:
            pixels, = np.where(mask)
        else:
            pixels = np.arange(len(p[0]))

    # rotate new coordinate system to original coordinates
    theta, phi = hp.pix2ang(res, pixels)
    mtheta, mphi = R(theta, phi)
    del theta
    del phi

    pr = np.full_like(p, np.nan)

    if use_qpoint:
        ra = np.degrees(mphi)
        dec = 90. - np.degrees(mtheta)
        del mtheta
        del mphi

    # interpolate Q and U (better before or after rot?)
    for idx in range(6):
        if use_qpoint:
            pr[idx, pixels] = Q.get_interp_val(p[idx], ra, dec)
        else:
            pr[idx, pixels] = hp.get_interp_val(p[idx], mtheta, mphi)

    rotate_qu_radial(
        pr, coord=coord, rot=rot, pol_axis=pol_axis, inplace=True)

    return pr

@core.usefulfunc
def rotate_map_alm(m, coord=['C', 'G'], verbose=True):
    """
    Rotate a map from one coordinate system to another using the `alteralm`
    fortran utility.

    Arguments
    ---------
    m : array_like
        A single map or three (I,Q,U) maps
    coord : list of two coordinates, optional.
        Coordinates to rotate between.  Default: ['C', 'G']

    Returns
    -------
    mr : array_like
        The input map rotated to the new coordinates.
    """

    nside = hp.get_nside(m)
    lmax = 3 * nside - 1
    alm = hp.map2alm(m, lmax=lmax, pol=True)
    alm_rot = alteralm(alm, nside, coord_in=coord[0], coord_out=coord[1],
                       double_precision=True, verbose=verbose)
    return np.asarray(
        hp.alm2map(alm_rot, nside, lmax=lmax, pol=True, fwhm=0., verbose=False))

@core.usefulfunc
def rotate_qu(m, ang, fill=np.nan):
    """
    Rotate the Q and U maps as if by a fixed instrument angle offset.
    Note: operation is in-place if m is a numpy array.

    Arguments:
    ----------
    m : array_like
        A map with shape (3, npix) for (T, Q, U) or (2, npix) for (Q, U)
    ang : float
        The angle by which to rotate (degrees).
    fill : scalar, optional
      Fill the bad pixels with this value, if supplied.  Default: NaN.

    Returns
    -------
    mr : array_like
        The map with Q and U rotated.
    """
    mr = np.atleast_2d(m)
    # don't need to handle NaN, since it actually handles this normally
    mask = np.logical_or(mask_bad(mr[-2]), mask_bad(mr[-1]))
    c = np.cos(2 * np.deg2rad(instang))
    s = np.sin(2 * np.deg2rad(instang))

    mr[-2], mr[-1] =  (mr[-2] * c + mr[-1] * s), (-mr[-2] * s + mr[-1] * c)
    mr[-2:, mask] = fill

    return mr

@core.usefulfunc
def rotate_qu_radial(m, m_Q = None, m_U = None, coord=['C', 'G'], rot=None,
                     pol_axis=[0.,0.,1.], inplace=False):
    """
    Rotate a map's Q/U components to radial Qr/Ur coordinates.

    Arguments
    ---------
    m : array_like
        map (Q/U or T/Q/U) or proj (TT/TQ/TU/QQ/QU/UU) to which the rotation is applied.
    m_Q: array_like, optional 
        same but stacked on Q.  Cannot be a proj map.
    m_U: array_like, optional 
        same but stacked on U.  Cannot be a proj map.
    coord : list of two coordinates, optional.
        Coordinates to rotate between.  Default: ['C', 'G']
    rot : scalar or sequence, optional
        Describe the rotation to apply.
        In the form (lon, lat, psi) (unit: degrees) : the point at
        longitude lon and latitude lat will be at the center of the rotated
        map. An additional rotation of angle psi around this direction is applied
    pol_axis : 3-vector, optional
        Axis normal to the plane in which the Q/U coordinates are defined.
    inplace : bool, optional
        If True, the rotation is applied in-place in memory.

    Returns
    -------
    m : array_like
        The input map, with the Q/U components rotated to radial
        Qr/Ur coordinates.
    """

    m = np.atleast_2d(m)
    if m.shape[0] not in [2, 3, 6]:
        raise ValueError("Input m must contain 2, 3, or six columns")

    proj = m.shape[0] == 6
    if proj:
        if m_Q is not None or m_U is not None:
            raise ValueError("Cannot use m_Q or m_U if input m is a proj map")
        iq, iu, iqq, iqu, iuu = range(-5, 0, 1)
    else:
        iq, iu = -2, -1

    if not inplace:
        m = m.copy()
    nside = hp.get_nside(m)
    pixels, = np.where(mask_good(m[0]))

    mq, mu = m[iq, pixels], m[iu, pixels]
    if proj:
        mqq, mqu, muu = m[iqq, pixels], m[iqu, pixels], m[iuu, pixels]
    elif m_Q is not None and m_U is not None:
        mq_q, mu_q = m_Q[iq, pixels], m_Q[iu, pixels]
        mq_u, mu_u = m_U[iq, pixels], m_U[iu, pixels]

    if rot is None:
        R = hp.Rotator(coord=coord, inv=True)
    else:
        R = hp.Rotator(rot=rot, inv=True)

    vec = hp.pix2vec(nside, pixels)
    mvec = R(vec)
    del vec

    vec0 = np.asarray(pol_axis)
    mvec0 = R(vec0)

    # calculate orientation of local meridian
    # based on Healpix rotate_coord
    mvec = np.asarray(mvec).T
    x = np.cross(vec0, mvec)
    cos_psi = np.dot(np.cross(x, mvec), mvec0)
    del mvec
    sin_psi = np.dot(x, mvec0)
    del x

    norm = sin_psi * sin_psi + cos_psi * cos_psi
    s2psi = 2. * sin_psi * cos_psi / norm
    del cos_psi
    c2psi = 1. - 2. * sin_psi * sin_psi / norm
    del norm
    del sin_psi

    # Rotate to new Q and U wrt to local meridian, in place
    m[iq, pixels], m[iu, pixels] = mq * c2psi + mu * s2psi, mu * c2psi - mq * s2psi

    s4psi = 2. * s2psi * c2psi
    del c2psi
    c4psi = 1. - 2. * s2psi * s2psi
    del s2psi

    if proj:
        ms = (mqq + muu) / 2.
        md = (mqq - muu) / 2.
        delta = md * c4psi + mqu * s4psi
        m[iqq, pixels], m[iqu, pixels], m[iuu, pixels] = (
            ms + delta, mqu * c4psi - md * s4psi, ms - delta )

    if m_Q is not None and m_U is not None:
        mqu = mq_u + mu_q
        md = mq_q - mu_u
        m_Q[iq, pixels], m_Q[iu, pixels] = mq_q + mu_u, mqu * s4psi + md * c4psi
        m_U[iq, pixels], m_U[iu, pixels] = mqu * c4psi - md * s4psi, md
        return m, m_Q, m_U
    else:
        return m

@core.usefulfunc
def rotate_qu_stack(m, m_Q = None, m_U = None, inplace=False):
    """
    Rotate a stack map's Q/U components to radial Qr/Ur coordinates.
    Assumes that the data are stacked along (theta, phi) = (pi/2, 0),
    and that the Q/U coordinate axes are defined normal to
    (theta, phi) = (0, 0).

    Arguments
    ---------
    m : array_like
        1D (T-only) or 2D (TQU) map to which the rotation is applied.
    m_Q: array_like, optional 
        same but stacked on Q
    m_U: array_like, optional 
        same but stacked on U
    inplace : bool, optional
        If True, the rotation is applied in-place in memory.

    Returns
    -------
    m : array_like
        The input map, with the Q/U components rotated to radial
        Qr/Ur coordinates.
    """
    return rotate_qu_radial(m, m_Q=m_Q, m_U=m_U, rot=[0., 90., 0.],
                            pol_axis=[1., 0., 0.], inplace=inplace)

@core.usefulfunc
def outline_map(m, low=None, high=None):
    """
    Create an outline of a map given a low and/or high amplitude threshold.
    The outline is defined as any pixels which have at least one neighbor
    with values within the threshold and not UNSEEN.

    Arguments
    ---------
    m : array_like
        A single map
    low : float, optional
        Value above which to include pixels inside the outline.
        Default: None (ignore)
    high : float, optional
        Value below which to include pixels inside the outline.
        Default: None (ignore)

    Returns
    -------
    mo : array_like
        Map the same size as `m` set to zero everywhere except the outline
        of the region in which the map values are above the threshold.
    """
    mo = np.zeros_like(m)
    nside = hp.get_nside(m)
    mask = mask_good(m)
    mm = m[mask]
    if low is None:
        low = mm.min() - 1
    if high is None:
        high = mm.max() + 1
    mb = np.zeros(m.shape, dtype=bool)
    mb[mask] = ((mm > low) & (mm < high)).astype(bool)
    pix = np.where(mb)[0]
    for p in pix:
        n = hp.get_all_neighbours(nside, p)
        n = n[n>=0]
        if np.any(mb[n]==0):
            mo[p] = 1
    return mo

@core.usefulfunc
def fill_outline(m, tol=None):
    """
    Fill an outline map.

    TODO: handle edge cases (outline wrapped around poles or meridians)

    Arguments
    ---------
    m : array_like
        Outline map, (0 where empty, 1 on the outline), such as one produced by
        `outline_map`.
    tol : float, optional
        Tolerance factor.  If `None`, determined from the map `nside`.

    Returns
    -------
    m : array_like
        Map with outlined region filled in.
    """

    nside = hp.get_nside(m)

    # tolerance pixel separation
    if tol is None:
        tol = np.sqrt(4. * np.pi / hp.nside2npix(nside) / 2.)

    m = np.asarray(m, dtype=bool)
    mf = m.copy()

    # outline pixels
    opix = np.where(m)[0]
    otheta, ophi =hp.pix2ang(nside, opix)

    # other pixels
    fpix = np.where(~m)[0]
    ftheta, fphi = hp.pix2ang(nside, fpix)

    # coarse crop
    inner = (fphi > ophi.min()) & (fphi < ophi.max()) & \
            (ftheta > otheta.min()) & (ftheta < otheta.max())
    fpix = fpix[inner]
    ftheta = ftheta[inner]
    fphi = fphi[inner]

    # better crop
    for pix, theta, phi in zip(fpix, ftheta, fphi):
        pidx = np.abs(ophi - phi) * np.sin(otheta) < tol
        if not np.any(pidx):
            continue
        theta_range = otheta[pidx]
        if theta < theta_range.min():
            continue
        if theta > theta_range.max():
            continue
        tidx = np.abs(otheta - theta) < tol
        if not np.any(tidx):
            continue
        phi_range = ophi[tidx]
        if phi < phi_range.min():
            continue
        if phi > phi_range.max():
            continue
        mf[pix] = True
    return mf

@core.usefulfunc
def latlon_mask(nside=512, latrange=(-70, -42), lonrange=(-50, 50), coord=None):
    """
    Create a mask that is True where lat and lon lie in a given range.  Default
    result is for the nominal SPT 1500d field in equatorial coordinates.

    Arguments
    ---------
    nside : int
        The nside at which to make the mask
    latrange : (min, max) 2-tuple of floats
        Min and max latitude in degrees. In the range -90 to +90
    lonrange : (min, max) 2-tuple of floats
        Min and max longitude in degrees. In the range -180 to +180
    coord : 2-tuple of strings (or list, etc)
        Pair of (input_coord, output_coord) as per rotate_map.
        For example, define a mask in equatorial coordinates and rotate
        into Galactic coordinates by setting this argument to ['C', 'G']

    Returns
    -------
    Boolean array/healpix map that is True in selected region
    """
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    lat = 90. - np.degrees(theta)
    lon = np.degrees(phi)
    lon[lon > 180.] -= 360.
    ret = np.logical_and(
            np.logical_and(lon > lonrange[0], lon < lonrange[1]),
            np.logical_and(lat > latrange[0], lat < latrange[1]))
    if coord is not None and len(coord) == 2:
        ret = rotate_map(ret, coord=coord).astype(bool)
    return ret

@core.usefulfunc
def circle_mask(nside=512, center=(-90.0, 0.0), radius=15.0, coord=None):
    """
    Create a mask that is True where lat and lon lie in a given range.
    Default result is a circle of radius 15 degress about the southern
    equatorial pole.

    Arguments
    ---------
    nside : int
        The nside at which to make the mask
    center : (lat, lon) 2-tuple of floats
        Coordinate to center the circular mask, in degrees.
    radius : float
        radius of the circle in degrees.
    coord : 2-tuple of strings (or list, etc)
        Pair of (input_coord, output_coord) as per rotate_map.
        For example, define a mask in equatorial coordinates and rotate
        into Galactic coordinates by setting this argument to ['C', 'G']

    Returns
    -------
    Boolean array/healpix map that is True in selected region
    """
    ret = np.zeros(hp.nside2npix(nside), dtype=bool)
    vec = hp.ang2vec(center[1], center[0], lonlat=True)
    ret[hp.query_disc(nside, vec, np.radians(radius))] = True
    if coord is not None and len(coord) == 2:
        ret = rotate_map(ret, coord=coord).astype(bool)
    return ret

@core.usefulfunc
def oval_mask(nside=512, center=(-90.0, 0.0), radius_lat=10.0, radius_lon=20.0,
              coord=None):
    """
    Create a mask that is True where lat and lon lie in a given range.
    Default result is an oval with lat radius 10 degrees, lon radius 20 degrees,
    centered on the southern equatorial pole.

    Arguments
    ---------
    nside : int
        The nside at which to make the mask
    center : (lat, lon) 2-tuple of floats
        Coordinate to center the ellipse, in degrees.
    radius_lat : float
        radius of the ellipse along the lat axis, in degrees.
    radius_lon : float
        radius of the ellipse along the lon axis, in degrees.
    coord : 2-tuple of strings (or list, etc)
        Pair of (input_coord, output_coord) as per rotate_map.
        For example, define a mask in equatorial coordinates and rotate
        into Galactic coordinates by setting this argument to ['C', 'G']

    Returns
    -------
    Boolean array/healpix map that is True in selected region
    """
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    lat = 90. - np.degrees(theta)
    lon = np.degrees(phi)
    lon[lon > 180.] -= 360.
    d1 = (center[0] - lat) / radius_lat
    d2 = (center[1] - lon) / radius_lon
    ret = 1 > d1**2. + d2**2.
    if coord is not None and len(coord) == 2:
        ret = rotate_map(ret, coord=coord).astype(bool)
    return ret

@core.usefulfunc
def ud_grade(map_in, nside_out, pess=False, order_in='RING', order_out=None,
             power=None, dtype=None, fill=np.nan):
    """Upgrade or degrade resolution of a map (or list of maps).

    in degrading the resolution, ud_grade sets the value of the superpixel
    as the mean of the children pixels.

    Parameters
    ----------
    map_in : array-like or sequence of array-like
      the input map(s) (if a sequence of maps, all must have same size)
    nside_out : int
      the desired nside of the output map(s)
    pess : bool
      if ``True``, in degrading, reject pixels which contains
      a bad sub_pixel. Otherwise, estimate average with good pixels
    order_in, order_out : str
      pixel ordering of input and output ('RING' or 'NESTED')
    power : float
      if non-zero, divide the result by (nside_in/nside_out)**power
      Examples:
      power=-2 keeps the sum of the map invariant (useful for hitmaps),
      power=2 divides the mean by another factor of (nside_in/nside_out)**2
      (useful for variance maps)
    dtype : type
      the type of the output map
    fill : scalar
      fill value for bad pixels

    Returns
    -------
    map_out : array-like or sequence of array-like
      the upgraded or degraded map(s)

    Examples
    --------
    >>> import healpy as hp
    >>> hp.ud_grade(np.arange(48.), 1)
    array([  5.5 ,   7.25,   9.  ,  10.75,  21.75,  21.75,  23.75,  25.75,
            36.5 ,  38.25,  40.  ,  41.75])
    """

    # mask bad values correctly
    map_in = np.asarray(map_in)
    bad = mask_bad(map_in)
    map_bad = map_in[bad]
    map_in[bad] = hp.UNSEEN

    map_out = np.asarray(hp.ud_grade(
        map_in, nside_out, pess=pess, order_in=order_in, order_out=order_out,
        power=power, dtype=dtype))
    map_out[mask_bad(map_out)] = fill

    # undo masking
    map_in[bad] = map_bad

    # return array
    return map_out


def spice(map_in, nlmax=None, mask=None, weight=None, fwhm=None, beam=None,
          apodizetype=None, apodizesigma=None, thetamax=None,
          decouple=False, pixwin=False, subav=False, subdipole=False,
          polarization=None, parfile=None, clfile=None, verbose=True,
          outroot=None, use_temp=False,
          kernelsfile=None, return_kernel=False, **kwargs):
    """Estimate cl(s) from input map(s).  Python interface to PolSpice.
    Input maps as arrays or filenames.  Intermediate files are stored
    in a temporary root that is discarded upon completion.

    Parameters
    ----------
    map_in : array or tuple of array
      A map or list of maps (3 maps for polarization)
    nlmax : int, scalar, optional
      Maximum l for cl.  Default: 3*nside - 1
    mask : array or string, optional
      Mask array, same shape as map, or string path to fits file
      Should be True for pixels to include, False for pixels to exclude.
    maskp : array or string, optional
      Polarization mask array (if different than mask), same shape as map,
      or string path to fits file.
    weight : array or string, optional
      Weight array, same shape as map, or string path to fits file
    weightp : array or string, optional
      Weight array for polarization maps (if different than weight),
      same shape as map, or string path to fits file.
    apodizetype : int, optional
      Apodization type, Gaussian (0) or cosine (1)
    apodizesigma : float, optional
      Apodization width (degrees). Default: disabled
      Values close to thetamax are recommended.
    thetamax : float, optional
      Maximum integration angle (degrees). Default: 180
    fwhm : float, optional
      beam width to deconvolve, in **radians**: Default: disabled
    beam : array or string, optional
      beam window function B_ell for correcting the output spectra,
      starting with the ell=0 mode.
    polarization : bool, optional
      If True, treat the input map as polarized.  If not supplied,
      assume True if `map_in` is a list of 3 maps.
    decouple : bool, optional
      If True, return T/E/B spectra, otherwise return T/Q/U spectra
      Default: False
    pixwin : bool or string, optional
      If True, apply default pixel window function for the map nside
      If string, apply supplied window function
      Default: False
    subav : bool, optional
      subtract best-fit monopole from the map
    subdipole : bool, optional
      subtract best-fit monopole and dipole from the map
    windowfile : string, optional
      If supplied and exists, sets the windowfilein parameter to avoid
      recomputing the window.  If supplied and missing, sets the
      windowfileout parameter to store the computed window for future use.
    parfile, mapfile, maskfile[p], weightfile[p], beamfile, clfile, pixelfile :
        string, optional
      Explicit file names for input and output files, stored
      in a location other than the temporary root.
    map[file]2, mask[file][p]2, weight[file][p]2, beam[file]2 :
        array_like or string, optional
      Second map(s) and weight functions for calculating cross-spectra
    symmetric_cl : bool, optional
    verbose : bool, optional
      print stuff or not.
    outroot : string, optional
      output prefix.  if not set, files are stored in a temporary directory
      created using tempfile.mkdtemp()
    use_temp : bool, optional
      allows one to specify an outroot, but still retain the temp directories
    kernelsfile : string or bool, optional
      If True, store the spice kernel in the output root.  Otherwise, should
      be an explicit filename to store in a location other than the output root.
      If False or None, the kernel is not stored.
    return_kernel : bool, optional
      If True, return the spice kernel along with the spectrum.

    Returns
    -------
    cls : array or list of arrays
      power spectrum computed by spice.
      Polarization order (with decouple=True): [TT, EE, BB, TE, TB, EB]
    kernel : array, optional
        Estimator kernel, returned if return_kernel is True.
    """

    executable = kwargs.pop('executable', spice_exe)

    # temporary directory
    if outroot is None or use_temp:
        tmproot = tf.mkdtemp(dir=outroot)
        filetag = ''
    else:
        tmproot, filetag = os.path.split(outroot)
        if not os.path.exists(tmproot):
            os.mkdir(tmproot)
    if filetag and not filetag.endswith('_'):
        filetag += '_'

    def opt(val):
        return 'YES' if val is True else 'NO' if val is False else val

    params = {}

    def get_file(filename, arg, filearg, altarg=None,
                 write_func=write_map):
        if isinstance(arg, basestring) and arg in kwargs:
            arg = kwargs.pop(arg)
        argname = filearg
        filearg = kwargs.pop(filearg, kwargs.pop(altarg, None))
        filename = os.path.join(tmproot, filetag + filename)
        if isinstance(arg, basestring):
            arg = None
            filearg = arg
        elif arg is True or arg is False:
            arg = None
            filearg = opt(arg)
        elif filearg is None and arg is not None:
            filearg = filename
            if write_func is write_map:
                arg = np.array(arg)
                arg[mask_bad(arg)] = hp.UNSEEN
            write_func(filearg, arg)
        if filearg is not None and filearg not in ['YES', 'NO'] and \
                not os.path.exists(filearg):
            raise OSError('{} not found'.format(filename))
        params[argname] = opt(filearg)

    # deal with input files
    get_file('map.fits', map_in, 'mapfile')
    if not isinstance(map_in, basestring) and polarization is None:
        polarization = len(map_in) == 3

    if clfile is None:
        clfile = os.path.join(tmproot, filetag + 'cl.fits')
    params['clfile'] = clfile

    get_file('mask.fits', mask, 'maskfile')
    get_file('maskp.fits', 'maskp', 'maskfilep')
    get_file('weight.fits', weight, 'weightfile')
    get_file('weightp.fits', 'weightp', 'weightfilep')
    get_file('beam.fits', beam, 'beam_file', 'beamfile',
             write_func=hp.write_cl)
    get_file('pixwin.fits', pixwin, 'pixelfile',
             write_func=hp.write_cl)
    if np.isscalar(fwhm):
        params['beam'] = np.degrees(fwhm) * 60

    get_file('map2.fits', 'map2', 'mapfile2')
    get_file('mask2.fits', 'mask2', 'maskfile2')
    get_file('maskp2.fits', 'maskp2', 'maskfilep2')
    get_file('weight2.fits', 'weight2', 'weightfile2')
    get_file('weightp2.fits', 'weightp2', 'weightfilep2')
    get_file('beam2.fits', 'beam2', 'beam_file2', 'beamfile2',
             write_func=hp.write_cl)
    fwhm2 = kwargs.pop('fwhm2', None)
    if np.isscalar(fwhm2):
        params['beam2'] = np.degrees(fwhm2) * 60
    params['symmetric_cl'] = opt(kwargs.pop('symmetric_cl', None))
    tolerance = kwargs.pop('tolerance', None)
    if np.isscalar(tolerance):
        params['tolerance'] = tolerance

    # create and populate parameter file
    if parfile is None:
        parfile = os.path.join(tmproot, filetag + 'spice.par')
    f = open(parfile, 'w')

    for k, v in params.items():
        if v is not None:
            f.write('{} = {}\n'.format(k, v))

    if nlmax is not None:
        f.write('nlmax = {:d}\n'.format(nlmax))
    f.write('polarization = {}\n'.format(opt(polarization)))
    f.write('decouple = {}\n'.format(opt(decouple)))
    f.write('subav = {}\n'.format(opt(subav)))
    f.write('subdipole = {}\n'.format(opt(subdipole)))
    if apodizetype in [0,1]:
        f.write('apodizetype = {:d}\n'.format(apodizetype))
    if np.isscalar(apodizesigma):
        f.write('apodizesigma = {:f}\n'.format(apodizesigma))
    if np.isscalar(thetamax):
        f.write('thetamax = {:f}\n'.format(thetamax))
    f.write('fits_out = YES\n')

    if verbose in [True, False]:
        f.write('verbosity = {}\n'.format(opt(verbose)))
    elif verbose in range(2):
        f.write('verbosity = {:d}\n'.format(verbose))

    windowfile = kwargs.pop('windowfile', None)
    if windowfile is not None:
        if windowfile is True:
            windowfile = os.path.join(tmproot, filetag + 'window.fits')
        if os.path.exists(windowfile):
            f.write('windowfilein = {}\n'.format(windowfile))
        else:
            f.write('windowfileout = {}\n'.format(windowfile))

    if return_kernel:
        if kernelsfile is None:
            kernelsfile = True
    if kernelsfile is not None and kernelsfile is not False:
        if kernelsfile is True:
            kernelsfile = os.path.join(tmproot, filetag + 'kernels.fits')
        f.write('kernelsfileout = {}\n'.format(kernelsfile))

    f.close()

    if len(kwargs.keys()):
        if outroot is None:
            shutil.rmtree(tmproot)
        raise TypeError("spice got unexpected keyword argument(s): {}".format(
                ", ".join(kwargs.keys())))

    # run spice
    stdout = None if verbose else open(os.devnull, 'w')
    try:
        sp.check_call([executable, '-optinfile', parfile], stdout=stdout,
                       stderr=sp.STDOUT)
    except (sp.CalledProcessError, KeyboardInterrupt):
        if outroot is None:
            shutil.rmtree(tmproot)
        raise
    if stdout:
        stdout.close()

    # read in output
    cls = np.asarray(hp.read_cl(clfile))
    if return_kernel:
        # do this by hand because kernel isn't stored in the standard hdu...
        import astropy.io.fits as pf
        hdus = pf.open(kernelsfile)
        kernel = hdus[0].data.astype(float)
        hdus.close()

    # cleanup
    if outroot is None or use_temp:
        shutil.rmtree(tmproot)

    # return
    if return_kernel:
        return cls, kernel
    return cls
