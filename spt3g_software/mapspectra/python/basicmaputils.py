"""
Basic map analysis utilities and functions.
All functions should work on plain numpy ndarrays 
if resolution is specified.

For more advanced functions, see map_analysis.py
"""
import numpy as np
import scipy.fftpack
from spt3g import core, coordinateutils
from spt3g.util.genericutils import is_float, is_int
from spt3g.util import maths

rad = core.G3Units.rad
deg = core.G3Units.deg

fft2_func = scipy.fftpack.fft2
ifft2_func = scipy.fftpack.ifft2


def check_arrs_size(*args):
    """
    Asserts that all of *args are either numeric scalars or
    arrays with the same shape.
    """
    assert len(args) > 0
    sh = np.shape(args[0])
    for a in args:
        assert is_float(a) or is_int(a) or sh == np.shape(a)


def _get_queb_trans_angle(m, res=None):
    """
    Helper function that calculates np.arctan2(kx, -ky) for all of the
    k-modes in the 2d fft of q and u for a given map size.  This is used
    with the flat sky definition of E/B.
    """
    if res is None:
        reso_rad = m.res / rad
    else:
        reso_rad = res / rad
    ny, nx = np.shape(m)
    npix = ny * nx
    lx, ly = np.meshgrid(
        np.fft.fftfreq(nx, reso_rad) * 2.0 * np.pi,
        np.fft.fftfreq(ny, reso_rad) * 2.0 * np.pi,
    )
    angles = np.arctan2(-lx, ly)
    return angles


@core.usefulfunc
def make_ellgrid(res, shape):
    """
    For a given resolution and shape computes the ell values corresponding to
    the elements in the 2d FFT.
    In particular this gives the ell bins computed by map_to_ft.
    """
    reso_rad = res / rad
    ell_x, ell_y = np.meshgrid(
        2 * np.pi * np.fft.fftfreq(shape[1], reso_rad),
        2 * np.pi * np.fft.fftfreq(shape[0], reso_rad),
    )
    ellgrid = np.hypot(ell_x, ell_y)
    return ellgrid


@core.usefulfunc
def map_to_ft(
    t, apod_mask=None, res=None, padded_map_shape=None
):
    """
    Computes the normalized 2d FFT of a map.

    Computes fft2(t*apod_mask)/scale_fac,
    where scale_fac corrects for the sky area and apodization mask to have the
    output values properly normalized in sqrt(C_l) units.

    Parameters
    ----------
    t : FlatSkyMap or ndarray
        The map for which to compute the 2d fourier transform

    apod_mask : FlatSkyMap or ndarray
        Multiplies the input map before the fourier transform is taken.

    res : float
        Resolution of the input map in G3 Units.
        Only required if input map is not a FlatSkyMap

    padded_map_shape : 2-tuple or 'square'
        The dimension to which q and u will be padded.
        Can be tuple like (npix_y, npix_x), or the string 'square'.

    Returns
    -------
    The 2d FFT of the input map, properly normalized in sqrt(C_l) units.

    See Also
    --------
    get_fft_scal_fac
    """
    if isinstance(t, coordinateutils.FlatSkyMap):
        res = t.res
        units = t.units
        proj = t.proj
    else:
        assert(res is not None), "Resolution not given."
        units = getattr(core.G3TimestreamUnits, "None")
        proj = getattr(coordinateutils.MapProjection, "ProjNone")
    if padded_map_shape is None:
        padded_map_shape = np.shape(t)
    elif padded_map_shape == "square":
        padded_map_shape = (max(np.shape(t)), max(np.shape(t)))
    t = pad_map(t, padded_map_shape)
    if apod_mask is None:
        apod_mask = np.ones(padded_map_shape)
    else:
        apod_mask = pad_map(apod_mask, padded_map_shape)
    n1, n2 = padded_map_shape
    scale_fac = get_fft_scale_fac(res, n1, n2, apod_mask)
    t_ft = fft2_func(t * apod_mask) / scale_fac
    map_nx = t.shape[1]
    map_ny = t.shape[0]
    from .map_spectrum_classes import MapSpectrum2D
    return MapSpectrum2D(map_nx, res, spec=t_ft, map_ny=map_ny, dy=res,
                         units=units, proj=proj)


def qu_to_eb_ft(q, u, apod_mask, res=None):
    """
    Computes the flat sky E/B 2d fft maps.  Scales them to be normalized to
    sqrt(C_l) based off of the apodization mask and covered sky area.

    Similar to map_to_ft but applies the flat sky definition of E/B from Q/U
    maps to compute the 2d k-space E/B maps.

    Uses the "basic" or naive QU combination for both E and B.
    Known to leak B->E (not a problem) and E->B (a problem).
    """
    if res is None:
        res = q.res
    check_arrs_size(q, u, apod_mask)
    n1, n2 = np.shape(q)
    scale_fac = get_fft_scale_fac(res, n1, n2, apod_mask)
    qu_to_eb_angle = _get_queb_trans_angle(q, res)
    q_ft = fft2_func(q * apod_mask) / scale_fac
    u_ft = fft2_func(u * apod_mask) / scale_fac
    e_ft = q_ft * np.cos(2 * qu_to_eb_angle) + u_ft * np.sin(
        2 * qu_to_eb_angle
    )
    b_ft = u_ft * np.cos(2 * qu_to_eb_angle) - q_ft * np.sin(
        2 * qu_to_eb_angle
    )
    return e_ft, b_ft


def eb_ft_to_qu(e_ft, b_ft, res, apod_mask=None):
    """
    Computes the inverse of qu_to_eb_ft.  Converts E/B 2d fft k space maps to
    real space q/u maps
    """
    if apod_mask is not None:
        check_arrs_size(e_ft, b_ft, apod_mask)
    else:
        check_arrs_size(e_ft, b_ft)
    n1, n2 = np.shape(e_ft)
    scale_fac = get_fft_scale_fac(res, n1, n2, apod_mask)
    qu_to_eb_angle = _get_queb_trans_angle(e_ft, res)
    qmap = e_ft * np.cos(2 * qu_to_eb_angle) - b_ft * np.sin(
        2 * qu_to_eb_angle
    )
    umap = e_ft * np.sin(2 * qu_to_eb_angle) + b_ft * np.cos(
        2 * qu_to_eb_angle
    )
    qmap = np.real(ifft2_func(qmap) * scale_fac)
    umap = np.real(ifft2_func(umap) * scale_fac)
    return qmap, umap


def get_fft_scale_fac(res, n1, n2, apod_mask=None):
    """
    Returns a scale factor to convert a 2d fft of a map into
    sqrt(C_l) normalized units by scaling it by the sky area.
    In particular it returns:
      np.mean(apod_mask**2)**0.5 * (n1 * n2)**0.5 / reso_rad

    For forward transforms: np.fft.fft(mp)  / scale_fac
    For inverse transforms: np.fft.ifft(mp) * scale_fac

    res: resolution of the map
    n1, n2: the map dimensions
    apod_mask: an apodization mask as an optional argument
    """
    reso_rad = res / rad
    if apod_mask is not None:
        maskfac = np.mean(apod_mask ** 2)
    else:
        maskfac = 1
    return maskfac ** 0.5 * (n1 * n2) ** 0.5 / reso_rad


def get_reg_spaced_ell_bins(cen_ells):
    """
    Returns list of (lower, upper) edges of bins centered on values
    in cen_ells. Bin width is fixed and set by the separation between
    the first and second values in cen_ells.
    """
    delta_ell = cen_ells[1] - cen_ells[0]
    return list(
        map(lambda x: (x - delta_ell / 2, x + delta_ell / 2), cen_ells)
    )


def bin_cls(ells, cls, ell_bins, ell_weights_1d=None):
    """
    Given Cls and the corresponding ells, average the Cls within the
    specified ell bins.
    """
    out_sums = np.zeros(len(ell_bins))
    wgt_sums = np.zeros(len(ell_bins))
    if ell_weights_1d is None:
        ell_weights_1d = cls * 0.0 + 1.0
    cls = cls * ell_weights_1d
    for i, eb in enumerate(ell_bins):
        inds = np.where(np.logical_and(ells >= eb[0], ells < eb[1]))
        out_sums[i] = np.sum(cls[inds])
        wgt_sums[i] = np.sum(ell_weights_1d[inds])
    out_sums /= wgt_sums
    return out_sums


def fill_center(m_to_fill, m):
    """
    Put m in the middle of m_to_fill.

    Assumes that the len of the the arrays along an axis are both even or
    both odd. That is to say if m has an even/odd axis len, m_to_fill should too,
    otherwise, with the simple rebinning the notion of center is a bit off.
    """
    s = np.shape(m)
    sf = np.shape(m_to_fill)
    np.asarray(m_to_fill)[
        sf[0] // 2 - s[0] // 2 : s[0] // 2 + sf[0] // 2,
        sf[1] // 2 - s[1] // 2 : s[1] // 2 + sf[1] // 2,
    ] = np.asarray(m)


def pad_map(m, new_size, fill=0):
    """
    Pads border of m with 'fill'.
    """
    if new_size == np.shape(m):
        return m
    new_map = np.zeros(new_size) + fill
    fill_center(new_map, m)
    return new_map


def degrade_map(m, df):
    """
    Reduce the resolution of m by integer factor df.
    """
    ms = np.shape(m)
    assert ms[0] % df == 0 or ms[1] % df == 0
    ns = list(np.shape(m))
    ns[0] = int(ns[0] / df)
    ns[1] = int(ns[1] / df)
    n = np.zeros(ns)
    for i in range(0, df):
        for j in range(0, df):
            n[:, :] += m[i::df, j::df]
    return n


def av_ps(
    mp_ft_sq, res, ell_bins, ell_weights_2d=None, tf2d=None, ellgrid=None
):
    """
    Averages the input 2d fourier map in the specified ell_bins
    """
    if ellgrid is None:
        ellgrid = make_ellgrid(res, np.shape(mp_ft_sq))
    if (ell_weights_2d) is None:
        if tf2d is None:
            spectrum, bins = np.histogram(
                ellgrid.ravel(), bins=ell_bins, weights=(mp_ft_sq).ravel()
            )
        else:
            spectrum, bins = np.histogram(
                ellgrid.ravel(),
                bins=ell_bins,
                weights=(mp_ft_sq / tf2d).ravel(),
            )
        norm, temp = np.histogram(ellgrid.ravel(), bins=ell_bins)

    else:
        if tf2d is None:
            spectrum, bins = np.histogram(
                ellgrid.ravel(),
                bins=ell_bins,
                weights=(mp_ft_sq * ell_weights_2d).ravel(),
            )
        else:
            spectrum, temp = np.histogram(
                ellgrid.ravel(),
                bins=ell_bins,
                weights=(mp_ft_sq * ell_weights_2d / tf2d).ravel(),
            )
        norm, temp = np.histogram(
            ellgrid.ravel(), bins=ell_bins, weights=ell_weights_2d.ravel()
        )
    spectrum[np.nonzero(norm)] /= norm[np.nonzero(norm)]
    return np.nan_to_num(spectrum)


@core.usefulfunc
def simple_cls(
    map1,
    map2=None,
    apod_mask=None,
    smooth=True,
    ell_bins=None,
    delta_ell=50,
    ell_min=None,
    ell_max=None,
    res=None,
    return_2d=False,
):
    """
    Computes the 2d power spectrum of map1, or if map2 is provided,
    the cross-power spectrum of map1 and map2.
    If return_2d = False, the 2d power spectrum is azimuthally averaged
    within the ell_bins, and the resulting 1d power spectrum is returned.

    If you're working with G3 objects instead of arrays or want
    a fancier function, use map_analysis.calculate_powerspectra

    Parameters:
    -----------
    map1 : array-like
        The 2-d array of which you want the power spectrum.
        If map2 = None, the autospectrum of map1 is returned.
    map2 [None]: array-like
        An optional second 2-d array used to calculate the
        cross-spectrum of map1 and map2. Must be same size as map1.
    apod_mask [None]: array-like
        An optional array, same size as map1, that smoothly
        goes to zero at the edges.
    smooth [True]: bool
        If apod_mask = None and smooth = True, smooth the
        input maps with a Hanning window before taking the FFT.
    ell_bins [None]: array-like
        A list of bin edges, like bins in np.histogram
    delta_ell [50]: float
        The spacing between ell bins.
    ell_min [None]: float
        If specified, create ell_bins down to ell_min.
    ell_max [None]: float
        If specified, create ell_bins up to ell_max.
    res [None]: float
        The resolution of the map, in G3 units. Default 1 arcminute
    return_2d [False]: bool
        If True, return the 2-d power spectrum.

    Returns:
    --------
    ells, cls : the ell values and either the 2-d or 1-d power spectrum.
    """
    if res is None:
        if hasattr(map1, "res"):
            res = map1.res
        else:
            res = 1 * core.G3Units.arcmin

    npix_y, npix_x = np.shape(map1)
    if ell_bins is None:
        if ell_min is None:
            ell_min = np.ceil(
                2 * np.pi / (np.max((npix_x, npix_y)) * res / core.G3Units.rad)
            )
        if ell_max is None:
            ell_max = 2 * np.pi / (res / core.G3Units.rad)
        ell_bins = np.arange(ell_min, ell_max + 1, delta_ell)
    ells = (np.array(ell_bins[1:]) + np.array(ell_bins[:-1]))/2

    if apod_mask is None:
        apod_mask = np.ones((npix_y, npix_x))
        if smooth:
            apod_mask *= maths.hanning(npix_y, npix_x)
    if np.shape(apod_mask) != np.shape(map1):
        raise ValueError("Map and apodization mask are different sizes.")

    ft = map_to_ft(map1, apod_mask, res=res)

    if map2 is None:
        ft = np.abs(ft) ** 2
    else:
        if np.shape(map2) != np.shape(map1):
            raise ValueError("Maps are different sizes.")
        xft = map_to_ft(map2, apod_mask, res=res)
        ft = (ft * np.conjugate(xft)).real

    if return_2d:
        ells = make_ellgrid(res, np.shape(map1))
        return ells, ft
    cls = av_ps(ft, res, ell_bins)
    return ells, cls

@core.usefulfunc
def smooth_flatsky_map(m, b_l, res=None, ell=None):
    """
    Smooth a flat sky map by an arbitrary function in ell-space

    Parameters:
    -----------
    m : G3 map object or array-like
        Map to be smoothed.
    b_l : array
        Ell-space vector to smooth m by.
    res : float
        Resolution in radians of m. If not provided, m must have res
        as an attribute.
    ell : array of shape b_l
        Ell vector corresponding to b_l. If None, assume b_l is for every
        ell starting from zero.
    """
    if res is None:
        res = m.res / core.G3Units.radians

    nx, ny = m.shape
    lx, ly = np.meshgrid(np.fft.fftfreq(nx, res), 
                         np.fft.fftfreq(ny, res))
    L = 2 * np.pi * np.sqrt(lx**2. + ly**2.)
    
    if ell is None:
        ell = np.arange(len(b_l))

    b_l_2d = np.interp(L.flatten(), ell, b_l, right=0.).reshape(L.shape).T
    map_filtered = np.fft.ifft2(np.fft.fft2(m) * b_l_2d).real 

    return map_filtered
