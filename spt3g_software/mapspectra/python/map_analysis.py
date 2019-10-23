"""
map_analysis.py

Functions for analyzing maps
"""
import numpy as np
import numexpr as ne
import scipy.fftpack
from spt3g import core, coordinateutils, mapmaker
from . import basicmaputils, fancyestimators, apodmask
from spt3g.util import files
import warnings
from .map_spectrum_classes import (
    MapSpectrum2D,
    MapSpectrum1D,
    MapSpectrum1DDict,
)


@core.usefulfunc
def flatten_pol(frame, invert=False):
    """
    For maps defined on the sphere the direction of the polarization angle is
    is defined relative to the direction of North.  When making maps we follow
    this definition.

    For any flat sky estimators, the polarization angle is defined relative to
    the vertical axis.  For some map projections the direction of north is not
    the same as the vertical axis.  This function applies a rotation to the Q
    and U values to switch the curved sky Q/U definition to the flat sky Q/U
    definition.

    If for whatever reason you want to reverse the process set the invert
    argument to True.

    """
    rad = core.G3Units.rad
    deg = core.G3Units.deg
    assert(
        "Q" in frame.keys() and "U" in frame.keys()
    ), "Flatten pol is requested, but the map doesn't have Q and U"
    q, u = frame['Q'], frame['U']
    # do not do anything if flatten/inverse flatten is done
    if hasattr(q, 'is_flattened') and q.is_flattened == True:
        if not invert:
            return frame
    elif hasattr(q, 'is_flattened') and q.is_flattened == False:
        if invert:
            return frame
    else:
    # this is the only function doing flatten_pol and set the attribute
    # if is_flattened is not set, we assume is_flattened is False
        if invert:
            return frame
    del frame['Q'], frame['U']
    ra0 = q.alpha_center
    dec0 = q.delta_center
    ra, dec = coordinateutils.maputils.get_ra_dec_map(q)

    # avoid branch cut at ra=0
    ra += 180.0 * deg - ra0
    ra %= 360.0 * deg
    ra /= rad
    # sanity test that no pixels are close to branch cut
    assert np.max(ra) < 350.0 * deg
    assert np.min(ra) > 10.0 * deg

    # convert to co-latitude
    dec *= -1.0 / rad
    dec += np.pi / 2.0

    # compute gradients
    dec_grd_x, dec_grd_y = np.gradient(dec)
    del dec
    arg = -1.0j * np.arctan2(dec_grd_y, dec_grd_x)
    del dec_grd_x, dec_grd_y
    ra_grd_x, ra_grd_y = np.gradient(ra)
    del ra
    arg += -1.0j * np.arctan2(-ra_grd_x, ra_grd_y)
    del ra_grd_x, ra_grd_y

    # compute rotation factor
    if invert:
        arg *= -1
    arg[:] = np.exp(arg)

    # apply rotation
    arg *= np.asarray(q) + 1.0j * np.asarray(u)

    qc = coordinateutils.FlatSkyMap(q)
    uc = coordinateutils.FlatSkyMap(u)
    np.asarray(qc)[:] = np.real(arg)
    np.asarray(uc)[:] = np.imag(arg)
    del arg
    if invert:
        qc.is_flattened = False
        uc.is_flattened = False
    else:
        qc.is_flattened = True
        uc.is_flattened = True
    frame['Q'], frame['U'] = qc, uc
    return frame


def deproject_tp(mp, band=None, q_scale=0.0, u_scale=0.0):
    """
    Correct for T-to-P monopole leakage

    Unweights the maps before the deprojection is done.

    Parameters
    ----------
    mp : Map frame
        Map frame containing T, Q, U, and weights.

    band : {None, 90, 150, 220}, optional
        If specified, overrides `q_scale` and `u_scale` with
        precomputed numbers.

    q_scale : float
        Multiply T map by this factor and subtract from Q map

    u_scale : float
        Multiply T map by this factor and subtract from U map

    Returns
    -------
    Input map frame with the unweighted, T-to-P-deprojected maps
    """
    mapmaker.mapmakerutils.RemoveWeightModule(mp)
    mapmaker.mapmakerutils.ZeroMapNans(mp)

    if band is not None:
        if int(band) == 90:
            q_scale = 0.009
            u_scale = 0.007
        elif int(band) == 150:
            q_scale = 0.011
            u_scale = 0.011
        elif int(band) == 220:
            q_scale = 0.015
            u_scale = 0.012
        else:
            raise ValueError(band)

    qmap = mp.pop('Q')
    mp['Q'] = qmap + -q_scale * mp['T']
    del qmap

    umap = mp.pop('U')
    mp['U'] = umap + -u_scale * mp['T']
    del umap

    return mp


def qhatUhatFromTwoMaps(
    map_fr1,
    map_fr2,
    apod_mask=None,
    uv_mask=None,
    laxis=350,
    lmin=500,
    lmax=2500,
    verbose=False,
):
    """
    Calculate the dimensionless T->P leakage coefficients, TQ/TT and TU/TT.

    Parameters
    ----------
    map_fr1: (Map frame) The map coadded from half of the data.
    map_fr2: (Map frame) The map coadded from the other half of the data.
    apod_mask: the apodization mask for these maps.

    uv_mask [None]: (2D array) The uv-space mask to apply when crossing maps.
        By default this program constructs a mask that is one everywhere
        except within LAXIS of the LX or LY axis. If you provide a uv_mask
        input and also set laxis>0, then the mask used will be the product of
        the uv_mask and the laxis mask.
        NB: I (RK) recommend using *only* the default uv_mask rather than
        supplying your own, because any guarantees about not biasing the TE
        spectrum are based on that default uv_mask (specifically, laxis=350).
    laxis [350]: used in constructing the default UV_MASK (see above).
    lmin [500]: the minimum ell used in estimating the TP coefficients.
    lmax [2500]: the maximum ell used in estimating the TP coefficients.

    Returns
    -------
    The output is a list [scaleQhat, scaleUhat], where:
        scaleQhat - the weighted average of TQ/TT.
        scaleUhat - the weighted average of TU/TT.
    AUTHOR
        Ryan Keisler, 23 October 2013
    CHANGES
        2018: Migrated to 3G software. DPD
    """
    # Figure out what input we were handed
    if isinstance(map_fr1, str):
        for frame in core.G3File(map_fr1):
            if frame.type == core.G3FrameType.Map:
                map_fr1 = frame
                break
    elif not isinstance(map_fr1, core.G3Frame):
        raise TypeError(type(map_fr1))

    if isinstance(map_fr2, str):
        for frame in core.G3File(map_fr2):
            if frame.type == core.G3FrameType.Map:
                map_fr2 = frame
                break
    elif not isinstance(map_fr2, core.G3Frame):
        raise TypeError(type(map_fr2))

    # Make the uv_mask.  By default this will mask out areas of
    # uv space within LAXIS of either the LX or LY axis.
    tshape = uv_mask.shape if uv_mask is not None else map_fr1["T"].shape
    uv_axis_mask = np.ones(tshape)
    if laxis > 0:
        if verbose:
            print(
                "Masking out all modes with |ell_X| < "
                + "%.2f or |ell_Y| < %.2f." % (laxis, laxis)
            )
        delta_l0 = (
            2.0 * np.pi / (tshape[0] * map_fr1["T"].res / core.G3Units.rad)
        )
        delta_l1 = (
            2.0 * np.pi / (tshape[1] * map_fr1["T"].res / core.G3Units.rad)
        )
        laxis_pix0 = int(laxis / delta_l0)
        laxis_pix1 = int(laxis / delta_l1)
        uv_axis_mask[0:laxis_pix0, :] = 0.0
        uv_axis_mask[-laxis_pix0:, :] = 0.0
        uv_axis_mask[:, 0:laxis_pix1] = 0.0
        uv_axis_mask[:, -laxis_pix1:] = 0.0
    if uv_mask is None:
        uv_mask = uv_axis_mask
    else:
        uv_mask *= uv_axis_mask

    # Take all cross-spectra.
    cl = calculate_powerspectra(
        map_fr1,
        input2=map_fr2,
        apod_mask=apod_mask,
        return_2d=True,
        flatten=False,
        qu_eb="qu",
    )
    l = cl["ell"]

    # We can measure TP/TT at every ell, but we'd like to return a single value.
    weight_th = uv_mask
    wh_zero = np.where((l < lmin) | (l > lmax))
    weight_th[wh_zero[0], wh_zero[1]] = 0.0
    weight_th /= np.sum(weight_th)

    # Calculate and return scalePhat==weighted sum of TP/TT.
    scaleQhat = np.sum(cl["TQ"] * weight_th) / np.sum(
        cl["TT"] * weight_th
    )
    scaleUhat = np.sum(cl["TU"] * weight_th) / np.sum(
        cl["TT"] * weight_th
    )

    if verbose:
        print(" ")
        print("TQ/TT = %0.4f" % scaleQhat)
        print("TU/TT = %0.4f" % scaleUhat)
        print(" ")

    return scaleQhat, scaleUhat


@core.usefulfunc
def constructEB(
    frame,
    res=None,
    apod_mask=None,
    e_mode_method="basic",
    b_mode_method="chi",
    padded_map_shape=None,
):
    """
    Takes (flattened) Q, U maps and construcs the 2-d Fourier transforms
    of the E and B maps via the specified methods.

    Parameters
    ----------
    frame: dict-like
        a map frame or a dictionary that is similar to a map frame.
    apod_mask [None]: array
        a 2-d array corresponding to an apodization mask.
    padded_map_shape [None]: string to tuble
        The dimension to which q and u will be padded.
        Can be tuple like (npix_y, npix_x), or the string 'square'.

    b_mode_method ['chi']: (string) How to calculate B?
        The preferred mode is "chi". We want to avoid E->B leakage, so we
        need to use a pure B-mode estimator. The "chi" and "smith" modes
        are mathematically identical, but the "chi" estimator is numerically
        nicer. Remember to use good apdodization windows!

        'chi' : Smith & Zaldarriaga chiB-method.
            This mode uses 'mapspectra.fancyestimators.calculateChi'.
        'smith': Original Kendrick Smith pure B-mode estimator.
            From 2006 Phys. Rev D paper.
            This mode uses 'mapspectra.fancyestimators.smith_pure_ft'.
        'basic' or anything else: Naive Q/U combination. Note that this will
            leak E->B.

    e_mode_method ['basic']: (string) How to calculate E?
        The preferred mode is "basic". We don't care about B->E leakage, and
        the pure estimators tend to have increased mode mixing.

        'chi':  Smith & Zaldarriaga chiB-method.
            This mode uses 'mapspectra.fancyestimators.calculateChi'.
        'smith': Original Kendrick Smith pure B-mode estimator(switched for E).
            From 2006 Phys. Rev D paper.
            This mode uses 'mapspectra.fancyestimators.smith_pure_ft'.
        'basic' or anything else: Naive Q/U combination. Note that this will
            leak B->E (but you probably don't care).

    Returns
    -------
    E, B 2d Fourier transforms
    """
    assert(
        'Q' in frame.keys() and 'U' in frame.keys()
    ), "Q and U not in the input frame."
    q, u = frame['Q'], frame['U']
    if isinstance(q, coordinateutils.FlatSkyMap):
        res = q.res
        units = q.units
        proj = q.proj
    else:
        assert(res is not None), "Resolution not given."
        units = getattr(core.G3TimestreamUnits, "None")
        proj = getattr(coordinateutils.MapProjection, "ProjNone")
    if padded_map_shape is None:
        padded_map_shape = np.shape(q)
    elif str(padded_map_shape) == "square":
        padded_map_shape = (max(np.shape(q)), max(np.shape(q)))
    # pad_map does nothing if the shape doesn't change
    q = basicmaputils.pad_map(q, padded_map_shape)
    u = basicmaputils.pad_map(u, padded_map_shape)
    if apod_mask is None:
        apod_mask = np.ones(padded_map_shape)
    else:
        apod_mask = basicmaputils.pad_map(apod_mask, padded_map_shape)
    ellgrid = basicmaputils.make_ellgrid(res, padded_map_shape)
    ellgrid[0, 0] = 1e6

    if e_mode_method == "basic" or b_mode_method == "basic":
        e_ft, b_ft = basicmaputils.qu_to_eb_ft(
            q, u, apod_mask, res
        )
    if e_mode_method == "chi":
        chiE = fancyestimators.calculateChi(
            q, u, res, which_chi="E", apodization_window=apod_mask
        )
        e_ft = basicmaputils.map_to_ft(chiE, apod_mask, res=res)
        e_ft *= fancyestimators.chi_scaling(ellgrid)
    elif e_mode_method == "smith":
        e_ft = fancyestimators.smith_pure_ft(
            q, u, apod_mask, res=res, do_e=True
        )
    elif e_mode_method == "basic":
        pass
    else:
        core.log_warn("Input E-mode method not recognized. Using basic method")

    if b_mode_method == "chi":
        chiB = fancyestimators.calculateChi(
            q, u, res, which_chi="B", apodization_window=np.asarray(apod_mask)
        )
        b_ft = basicmaputils.map_to_ft(chiB, apod_mask, res=res)
        b_ft *= fancyestimators.chi_scaling(ellgrid)
    elif b_mode_method == "smith":
        b_ft = fancyestimators.smith_pure_ft(
            q, u, apod_mask, res=res, do_e=False
        )
    elif b_mode_method == "basic":
        pass
    else:
        core.log_warn("Input B-mode method not recognized. Using basic method")
    map_nx = q.shape[1]
    map_ny = q.shape[0]
    dx = res
    dy = dx
    return {
        "E": MapSpectrum2D(
            map_nx, dx, spec=e_ft, map_ny=map_ny, dy=dy, units=units, proj=proj
        ),
        "B": MapSpectrum2D(
            map_nx, dx, spec=b_ft, map_ny=map_ny, dy=dy, units=units, proj=proj
        ),
    }


@core.usefulfunc
def mapffts_to_powerspectra(
    input1,
    lbins=None,
    input2=None,
    return_2d=False,
    realimag="real",
    tf_2d=None,
    ell_weights_2d=None,
    symmetric_cross_spectra=True,
    lmin=0,
    lmax=10000,
    delta_l=1,
    calculate_dls=False,
):
    """ Calculates the 2D or 1D power spectrum

    Calculate the annulus-averaged auto- or cross-spectrum of
    input maps or map ffts, for bins described by lbins.

    Parameters
    ----------
    input1: MapSpectrum2D object, or a dictionary.
        The map FFTs can be a MapSpectrum2D object, or a dictionary
        of MapSpectrum2D objects. The dictionary can have keys
        t/q/u, t/e/b, or t/q/u/e/b.
    lbins: 1D numpy array, or 'henning'
        A list of bin edges. If 'henning' is given, will do
        henning et al 2018 binning.
    input2 [None]: Same as input1
        Optional second input of the same type to cross-correlate
        with input1 (otherwise will return an auto-spectrum).
    return_2d [False]: boolean.
        Whether to return the 2d spectra
    realimag ['real']: string.
        Option to return the real or imaginary part (or both).
        Options are 'real', 'imag', 'both'. By default the code returns
        the complex power spectrum if return_2d is True.
    tf_2d [None]: a 2d array or a dictionary
        Used when dividing out the 2-d transfer function to unbias
        the power spectrum. Either a 2-d array to be used for all power
        spectra, or a dictionary of 2-d arrays keyed by map polarization
        ('T', 'E', 'B', etc.)
    ell_weights_2d [None] : a 2d array or a dictionary
        Used when summing the 2d power spectrum. Either a 2-d array
        to be used for all power spectra, or a dictionary
        of 2-d arrays keyed by map polarization ('T', 'E', 'Q', etc.)
    symmetric_cross_spectra [True]: boolean.
        Average equivalent cross-spectra, e.g. te and et, and
        return only one of them.
    lmin [0]: int or float
        The minimum ell to use, if lbins not provided
    lmax [10000]: int or float
        The maximum ell to use, if lbins not provided
    delta_l [1]: int or float
        The width of ell bins to use, if lbins not provided
    calculate_dls: boolean
        whether to calculate dls

    Returns
    -------
        1D or 2D auto-spectra or cross-spectra
            The 1D spectra are in format of Mapspectrum1D or
            a dictionary filled with Mapspectrum1Ds
            The 2D spectra are in format of
            {'ell': 2D array, 'TT': MapSpectrum2D, ...}
    """

    equivalent_ps = {
        "TE": "ET", "EB": "BE", "TB": "BT",
        "TQ": "QT", "QU": "UQ", "TU": "UT",
    }
    useful_ps = {
        "TE", "ET", "EB", "BE", "TB", "BT",
        "TQ", "QT", "QU", "UQ", "TU", "UT",
        "TT", "EE", "BB", "QQ", "UU",
    }
    # get transfer function
    if not isinstance(tf_2d, dict):
        tmp = dict()
        for k in ["T", "Q", "U", "E", "B"]:
            tmp[k] = tf_2d
        tf_2d = tmp
    # get the ell weights for 2d
    if not isinstance(ell_weights_2d, dict):
        tmp = dict()
        for k in ["T", "Q", "U", "E", "B"]:
            tmp[k] = ell_weights_2d
        ell_weights_2d = tmp

    if realimag in ['real', 'imag']:
        realimag = '.'+realimag
    elif realimag == 'complex':
        realimag = ''
    else:
        raise ValueError("realimag must be 'real', 'imag', or 'complex'.")

    # Define lbins
    if lbins is None:
        # By default the binning is from 0 to 10000 with spacing one
        lbins = np.arange(lmin, lmax + 1, delta_l)
    elif str(lbins) == "henning":
        # Henning et al. (2018) bins
        lbins = []
        div1, div2, div3, div4 = 50, 100, 500, 1000
        np.concatenate((lbins, np.arange(51, 2001, div1)))
        np.concatenate((lbins, np.arange(2001, 3001, div2)))
        np.concatenate((lbins, np.arange(3001, 5001, div3)))
        np.concatenate((lbins, np.arange(5001, 12001, div4)))

    if input2 is None:
        input2 = input1
    # If one single MapSpectrum2D instead of a dict is given as input
    # assume it's a T map by default.
    # Even if you have one single input map, it's recommended to give it a key
    # and store it in a dictionary to indicate what it is.
    if isinstance(input1, MapSpectrum2D) and isinstance(input2, MapSpectrum2D):
        input1 = {"T": input1}
        input2 = {"T": input2}
    # if the inputs are dictionaries of MapSpectrum2D objects
    if isinstance(input1, dict) and isinstance(input2, dict):
        k = next(iter(input1))
        assert (
            input1[k].compatible(input2[k])
        ), "Input1 and input2 are not compatible."
        res = input1[k].dx
        units = input1[k].units
        proj = input1[k].proj
        # calculate 2D the power spectrum
        ps = dict()
        ell = input1[k].get_ell()
        map_nx, map_ny = input1[k].map_nx, input1[k].map_ny
        dx, dy = input1[k].dx, input1[k].dy
        pi = np.pi
        for pol, ft in input1.items():
            for xpol, xft in input2.items():
                # do not do things such as QB spectrum, only useful ones
                if pol + xpol in useful_ps:
                    xft = np.conjugate(xft)
                    if calculate_dls:
                        spec2d = ne.evaluate(
                            "(ft*xft*ell*(ell+1)/2/pi)%s" % realimag)
                        spec_type = "dl"
                    else:
                        spec2d = ne.evaluate("(ft*xft)%s" % realimag)
                        spec_type = "cl"
                    del xft
                    if return_2d:
                        ps["ell"] = ell
                        ps[pol + xpol] = MapSpectrum2D(
                            map_nx, dx, spec=spec2d, map_ny=map_ny, dy=dy,
                            units=units, proj=proj
                        )
                    else:
                        # set ell weight of the cross spectra to be the
                        # product of the ell weights of the two maps.
                        if ell_weights_2d[pol] is not None:
                            ell_weight = (
                                ell_weights_2d[pol] * ell_weights_2d[xpol]
                            )
                        else:
                            ell_weight = None
                        # set transfer function of the cross spectra to be the
                        # product of the transfer functions of the two maps.
                        if tf_2d[pol] is not None:
                            tf = tf_2d[pol] * tf_2d[xpol]
                        else:
                            tf = None

                        spectrum = basicmaputils.av_ps(
                            spec2d, dx, lbins, ell_weight, tf, ell
                        )

                        ps[pol + xpol] = MapSpectrum1D(
                            lbins,
                            spectrum,
                            spec_type=spec_type,
                            dx=dx,
                            map_nx=map_nx,
                            map_ny=map_ny,
                            units=units,
                            proj=proj,
                            is_asd=False
                        )
            del ft

        if symmetric_cross_spectra:
            for spec_a, spec_b in equivalent_ps.items():
                if spec_a in ps and spec_b in ps:
                    ps[spec_a] = (ps[spec_a] + ps[spec_b]) / 2.0
                    del ps[spec_b]
        return ps


@core.usefulfunc
def calculate_powerspectra(
    input1,
    input2=None,
    lbins=None,
    delta_l=1,
    lmin=0,
    lmax=10000,
    return_2d=False,
    realimag="real",
    tf_2d=None,
    ell_weights_2d=None,
    symmetric_cross_spectra=True,
    apod_mask=None,
    padded_map_shape=None,
    flatten=True,
    e_mode_method="basic",
    b_mode_method="chi",
    calculate_dls=False,
    qu_eb="eb",
):
    """ Convert a map frame into power spectra

    Parameters
    ----------
    input1: G3Frame.
        The input raw maps.
    input2 [None]: Same as input1
        Optional second input of the same type to cross-correlate
        with input1 (otherwise will return an auto-spectrum).
    lbins: 1D numpy array
        A list of bin edges.
    delta_l [1]: int or float
        The width of ell bins to use, if lbins not provided
    lmin [0]: int or float
        The minimum ell to use, if lbins not provided
    lmax [10000]: int or float
        The maximum ell to use, if lbins not provided
    return_2d [False]: boolean.
        Whether to return the 2d spectra
    realimag ['real']: string.
        Option to return the real or imaginary part (or both).
        Options are 'real', 'imag', 'both'. By default the code returns
        the complex power spectrum if return_2d is True.
    tf_2d [None]: a 2d array or a dictionary
        Used when dividing out the 2-d transfer function to unbias
        the power spectrum. Either a 2-d array to be used for all power
        spectra, or a dictionary of 2-d arrays keyed by map polarization
        ('T', 'E', 'B', etc.)
    ell_weights_2d [None] : a 2d array or a dictionary
        Used when summing the 2d power spectrum. Either a 2-d array
        to be used for all power spectra, or a dictionary
        of 2-d arrays keyed by map polarization ('T', 'E', 'Q', etc.)
    symmetric_cross_spectra [True]: boolean.
        Average equivalent cross-spectra, e.g. te and et, and
        return only one of them.
    apod_mask [None]: FlatSkyMap, 2d array, or 'from_weight'.
        Apodization mask. Could be 'from_weight', in which
        case a nominal apodization mask is calculated.
        If left as None, no masking is done.
    padded_map_shape [None]: 2d array [npix_y, npix_x] or 'square'
        Specify the shape of the fft map. Pad maps to this size.
        If set to square, the larger dimension of the map is used
        for both npix_x, npix_y.
    flatten [True]: boolean
        Flatten the input polarized maps before calculating
        power spectra. Set to False for calculating T->P leakage.
    e_mode_method ['basic']: string
        How to calculate the E mode. Options are "chi", "smith",
        and "basic". See constructEB docstring
    b_mode_method ['chi']: string
        How to calculate the B mode. Options are "chi", "smith",
        and "basic". See constructEB docstring
    calculate_dls: boolean
        whether to calculate dls
    qu_eb ['eb']: string
        Options are 'qu', 'eb', or 'both'. Return the MapSpectrum2D
        objects for qu, eb, or queb.

    Returns
    -------
        1D or 2D auto-spectra or cross-spectra
            The 1D spectra are in format of Mapspectrum1D or
            Mapspectrum1D
            The 2D spectra are in format of
            {'ell': 2D array, 'TT': MapSpectrum2D, ...}
    """
    # Set up apod mask
    if str(apod_mask) == "from_weight":
        # radius_arcmin=90 is used instead of default 20
        # in case these are single-observation maps
        for k in input1:
            if isinstance(input1[k], coordinateutils.G3SkyMapWeights):
                wgt_k = k
        apod_mask = apodmask.make_border_apodization(
            input1[wgt_k], radius_arcmin=90.0
        )
    elif apod_mask is not None and not hasattr(apod_mask, 'shape'):
        raise TypeError(
            "apod_mask must be array_like, None, or 'from_weight'.")
    # RemoveWeightModule will validate the input frame.
    # It won't do anything more if map is not weighted
    mapmaker.mapmakerutils.RemoveWeightModule(input1)
    mapmaker.mapmakerutils.ZeroMapNans(input1)
    if flatten and 'Q' in input1.keys() and 'U' in input1.keys():
        flatten_pol(input1)
    mapffts1 = {}
    if qu_eb == "qu" or qu_eb == "both":
        pol_types = {"T", "Q", "U"}
    if qu_eb == "eb":
        pol_types = {"T"}
    for pol_type in pol_types:
        if pol_type in input1.keys():
            mapffts1[pol_type] = basicmaputils.map_to_ft(
                input1[pol_type], apod_mask,
                padded_map_shape=padded_map_shape
            )
    if qu_eb == 'eb' or qu_eb == 'both':
        if "Q" in input1.keys() and "U" in input1.keys():
            mapffts1.update(constructEB(
                input1,
                apod_mask=apod_mask,
                e_mode_method=e_mode_method,
                b_mode_method=b_mode_method,
                padded_map_shape=padded_map_shape))
    del input1

    if input2 is not None:
        mapmaker.mapmakerutils.RemoveWeightModule(input2)
        mapmaker.mapmakerutils.ZeroMapNans(input2)
        if flatten and 'Q' in input2.keys() and 'U' in input2.keys():
            flatten_pol(input2)
        mapffts2 = {}
        if qu_eb == "qu" or qu_eb == "both":
            pol_types = {"T", "Q", "U"}
        if qu_eb == "eb":
            pol_types = {"T"}
        for pol_type in pol_types:
            if pol_type in input2.keys():
                mapffts2[pol_type] = basicmaputils.map_to_ft(
                    input2[pol_type], apod_mask,
                    padded_map_shape=padded_map_shape
                )
        if qu_eb == 'eb' or qu_eb == 'both':
            if "Q" in input2.keys() and "U" in input2.keys():
                mapffts2.update(constructEB(
                    input2,
                    apod_mask=apod_mask,
                    e_mode_method=e_mode_method,
                    b_mode_method=b_mode_method,
                    padded_map_shape=padded_map_shape))
        del input2
    else:
        mapffts2 = None

    ps = mapffts_to_powerspectra(
        mapffts1,
        lbins=lbins,
        input2=mapffts2,
        return_2d=return_2d,
        realimag=realimag,
        tf_2d=tf_2d,
        ell_weights_2d=ell_weights_2d,
        symmetric_cross_spectra=symmetric_cross_spectra,
        lmin=lmin,
        lmax=lmax,
        delta_l=delta_l,
        calculate_dls=calculate_dls,
    )
    return ps


@core.usefulfunc
def filterMap(map_fr, uv_filter, apod_mask=None):
    """
    Takes an input Map frame, applies the specified filter in k-space,
    then returns a Map frame containing the result.

    Parameters
    ----------
    map_fr: dict-like
        The input map. It will not be modified.
    uv_filter: 2D array
        Must be the same size as the map. We'll multiply
        the result of FFTing the map by this array, then FFT back.
        Note that this uv_filter array should be in the conventional FFT
        ordering.
    apod_mask [None]: 2D array
        If not None, then apodize the input map by
        this mask before FFTing. You should use a good apodization!

    Returns
    -------
    A Map frame in which the T, Q, and U maps have been filtered in k-space
    by the input uv_filter.
    """

    if (
        isinstance(map_fr, core.G3Frame)
        and map_fr.type == core.G3FrameType.Map
    ):
        pass
    else:
        raise ValueError("I only know how to filter a Map frame.")

    filtered_map = core.G3Frame(core.G3FrameType.Map)
    if "Id" in map_fr:
        filtered_map["Id"] = "filtered_" + map_fr["Id"]

    # Create a default "mask" of all ones if we didn't get an input mask
    if apod_mask is None:
        apod_mask = np.ones(map_fr["T"].shape)

    if "Wunpol" in map_fr:
        if map_fr["T"].is_weighted:
            t = mapmaker.mapmakerutils.remove_weight_t(
                map_fr["T"], map_fr["Wunpol"]
            )
        else:
            t = fr["T"]
        maps = {"T": t}
        filtered_map["Wunpol"] = map_fr["Wunpol"]

    elif "Wpol" in map_fr:
        if (
            map_fr["T"].is_weighted
            and map_fr["Q"].is_weighted
            and map_fr["U"].is_weighted
        ):
            t, q, u = mapmaker.remove_weight(
                map_fr["T"], map_fr["Q"], map_fr["U"], map_fr["Wpol"]
            )
        elif not (
            map_fr["T"].is_weighted
            or map_fr["Q"].is_weighted
            or map_fr["U"].is_weighted
        ):
            t, q, u = map_fr["T"], map_fr["Q"], map_fr["U"]
        else:
            raise ValueError(
                "Map has combination of weighted and unweighted maps"
            )
        maps = {"T": t, "Q": q, "U": u}
        filtered_map["Wpol"] = fr["Wpol"]
    else:
        raise KeyError("No weights in map frame (Wpol or Wunpol)")

    for k, pol in maps.items():
        ft = basicmaputils.map_to_ft(pol, apod_mask=apod_mask)
        ft *= uv_filter
        marr = np.asarray(pol)
        marr[:] = np.fft.ifft2(ft).real
        filtered_map[k] = pol

    return filtered_map


def subtract_two_maps(
    map1, map2, divide_by_two=False, in_place=False, mapid=None
):
    """
    Compute the difference of two maps.

    Parameters
    ----------
    map1: str or G3Frame
        The first map. Can be .g3 file or Map frame
    map2: str or G3Frame
        The second map. Can be .g3 file or Map frame
    divide_by_two: bool
        Does what it says. Used for e.g. calculating noise.
    in_place: boool
        If the input frames contain weighted maps, this will
        modify them in-place with the unweighted maps.
        Reduces memory usage by 30%.
    mapid : str
        if mapid is not none and map file given, read in
        the frame in the file with the specified mapid.

    Returns
    -------
    Map frame containing the (unweighted) results of map1 - map2
    """
    # This function may look long, but it's just doing a lot of checks
    # so that the user can't do anything too stupid.

    # Load in first map
    if isinstance(map1, str):
        obs = core.G3File(map1)
        for frame in obs:
            if frame.type == core.G3FrameType.Map and (
                mapid == None or frame["Id"] == mapid
            ):
                map1 = frame
    elif not isinstance(map1, core.G3Frame):
        raise TypeError()
    # Remove weights if needed
    pol = False
    if map1["T"].is_weighted:
        if in_place:
            mapmaker.mapmakerutils.RemoveWeightModule(map1)
            mapmaker.mapmakerutils.ZeroMapNans(map1)
            t1 = map1["T"]
            if "Q" in map1:
                pol = True
                q1, u1 = map1["Q"], map1["U"]
        elif "Wunpol" in map1:
            t1 = mapmaker.mapmakerutils.remove_weight_t(
                map1["T"], map1["Wunpol"]
            )
        elif "Wpol" in map1:
            pol = True
            if not map1["Q"].is_weighted or not map1["U"].is_weighted:
                raise ValueError(
                    "Map frame has combination of weighted, unweighted maps"
                )
            t1, q1, u1 = mapmaker.remove_weight(
                map1["T"], map1["Q"], map1["U"], map1["Wpol"]
            )
        else:
            raise KeyError("Weights not found in frame of weighted map")
    else:
        t1 = map1["T"]
        if "Wpol" in map1:
            pol = True
            q1, u1 = map1["Q"], map1["U"]
    # Load in second map
    if isinstance(map2, str):
        obs = core.G3File(map2)
        for frame in obs:
            if frame.type == core.G3FrameType.Map and (
                mapid == None or frame["Id"] == mapid
            ):
                map2 = frame
    elif not isinstance(map2, core.G3Frame):
        raise TypeError()

    # Remove weights if needed
    if map2["T"].is_weighted:
        if in_place:
            mapmaker.mapmakerutils.RemoveWeightModule(map2)
            mapmaker.mapmakerutils.ZeroMapNans(map2)
            t2 = map2["T"]
            if "Q" in map2:
                q2, u2 = map2["Q"], map2["U"]
        elif "Wunpol" in map2:
            t2 = mapmaker.mapmakerutils.remove_weight_t(
                map2["T"], map2["Wunpol"]
            )
        elif "Wpol" in map2:
            if not map2["Q"].is_weighted or not map2["U"].is_weighted:
                raise ValueError(
                    "Map frame has combination of weighted, unweighted maps"
                )
            t2, q2, u2 = mapmaker.remove_weight(
                map2["T"], map2["Q"], map2["U"], map2["Wpol"]
            )
        else:
            raise KeyError("Weights not found in frame of weighted map")
    else:
        t2 = map2["T"]
        if "Wpol" in map2:
            q2, u2 = map2["Q"], map2["U"]
    norm = 1
    if divide_by_two:
        norm = 2
    # Compute and return difference maps
    difference_map = core.G3Frame(core.G3FrameType.Map)
    t1 -= t2
    difference_map["T"] = t1 / norm
    if not pol:
        try:
            difference_map["Wunpol"] = map1["Wunpol"] + map2["Wunpol"]
        except KeyError:
            pass
    if pol:
        try:
            q1 -= q2
            u1 -= u2
            difference_map["Q"] = q1 / norm
            difference_map["U"] = u1 / norm
            difference_map["Wpol"] = map1["Wpol"] + map2["Wpol"]
        except NameError:
            print("Are you differencing polarized and unpolarized maps?")
            print("Only including T in the output.")

    return difference_map


def apply_weight(frame):
    """
    Function to weight a map.
    If the map is already weighted, return the original map.

    Parameters
    ----------
    frame: weighted or unweighted map frame

    Returns
    -------
    weighted map frame
    """
    if frame["T"].is_weighted:
        return frame
    if "Wpol" in frame:
        t = (
            frame["T"] * frame["Wpol"].TT
            + frame["Q"] * frame["Wpol"].TQ
            + frame["U"] * frame["Wpol"].TU
        )
        q = (
            frame["T"] * frame["Wpol"].TQ
            + frame["Q"] * frame["Wpol"].QQ
            + frame["U"] * frame["Wpol"].QU
        )
        u = (
            frame["T"] * frame["Wpol"].TU
            + frame["Q"] * frame["Wpol"].QU
            + frame["U"] * frame["Wpol"].UU
        )
        del frame["T"]
        del frame["Q"]
        del frame["U"]
        frame["T"] = t
        frame["Q"] = q
        frame["U"] = u
        frame["T"].is_weighted = True
        frame["Q"].is_weighted = True
        frame["U"].is_weighted = True
    elif "Wunpol" in frame:
        t = frame["T"] * frame["Wunpol"].TT
        del frame["T"]
        frame["T"] = t
        frame["T"].is_weighted = True
    return frame


def add_two_maps(map1, map2, mapid=None, recover_input=False):
    """
    Compute the sum of two maps.

    Parameters
    ----------
    map1: str or G3Frame
        The first map. Can be .g3 file or Map frame
    map2: str or G3Frame
        The second map. Can be .g3 file or Map frame
    mapid : str
        if mapid is not none and map file given, read in
        the frame in the file with the specified mapid.
    recover_input : bool
        the maps needs to be weighted before summing.
        if set to true, unweight the input maps if they
        were unweighted.
        if you do not care, set to false to save time

    Returns
    -------
    Map frame containing the (weighted) results of map1 - map2
    """
    # This function may look long, but it's just doing a lot of checks
    # so that the user can't do anything too stupid.

    # Load in first map
    if isinstance(map1, str):
        obs = core.G3File(map1)
        for frame in obs:
            if frame.type == core.G3FrameType.Map and (
                mapid == None or frame["Id"] == mapid
            ):
                map1 = frame
    elif not isinstance(map1, core.G3Frame):
        raise TypeError()

    # Weight first map if needed
    weighted_1 = True
    if not map1["T"].is_weighted:
        map1 = apply_weight(map1)
        weighted_1 = False

    # Load in second map
    if isinstance(map2, str):
        obs = core.G3File(map2)
        for frame in obs:
            if frame.type == core.G3FrameType.Map and (
                mapid == None or frame["Id"] == mapid
            ):
                map2 = frame
    elif not isinstance(map2, core.G3Frame):
        raise TypeError()

    # Weight second map if needed
    weighted_2 = True
    if not map2["T"].is_weighted:
        map2 = apply_weight(map2)
        weighted_2 = False

    # Calculate the summed map
    summed = core.G3Frame(core.G3FrameType.Map)
    if "Wpol" in map1 and "Wpol" in map2:
        summed["T"] = map1["T"] + map2["T"]
        summed["Q"] = map1["Q"] + map2["Q"]
        summed["U"] = map1["U"] + map2["U"]
        summed["Wpol"] = map1["Wpol"] + map2["Wpol"]
    elif "Wunpol" in map1 and "Wunpol" in map2:
        summed["T"] = map1["T"] + map2["T"]
        summed["Wunpol"] = map1["Wunpol"] + map2["Wunpol"]
    else:
        raise ValueError("Adding polarized and unpolarized maps?")

    # Recover the input maps into the original state before weighting
    if recover_input:
        if not weighted_1:
            mapmaker.mapmakerutils.RemoveWeightModule(map1)
            mapmaker.mapmakerutils.ZeroMapNans(map1)
        if not weighted_2:
            mapmaker.mapmakerutils.RemoveWeightModule(map2)
            mapmaker.mapmakerutils.ZeroMapNans(map2)
    # If one of the input map is unweighted, unweight the sum
    if (not weighted_1) or (not weighted_2):
        mapmaker.mapmakerutils.RemoveWeightModule(summed)
        mapmaker.mapmakerutils.ZeroMapNans(summed)
        return summed
    else:
        return summed


def psd_2d_to_1d(
    psd_2d, lmin=0, lmax=10000, delta_l=1, lbins=None, spec_type="cl"
):
    """
    Bin MapSpectrum2D (2d power spectrum)
    to MapSpectrum1D (1d power spectrum)
    """
    if lbins is None:
        # By default the binning is from 0 to 10000 with spacing one
        lbins = np.arange(lmin, lmax + 1, delta_l)
    if spec_type == "dl":
        ell = psd_2d.get_ell()
        psd_2d = psd_2d * ell * (ell + 1) / 2 / np.pi
    spectrum = basicmaputils.av_ps(psd_2d, psd_2d.dx, lbins)
    return MapSpectrum1D(
        lbins,
        spectrum,
        spec_type="cl",
        dx=psd_2d.dx,
        map_nx=psd_2d.map_nx,
        map_ny=psd_2d.map_ny,
        units=psd_2d.units,
        proj=psd_2d.proj,
        is_asd=False,
    )


@core.usefulfunc
def calculateNoise(
    map1,
    map2=None,
    apod_mask="from_weight",
    ell_range=[5500.0, 6500.0],
    return_all=False,
    verbose=True,
):
    """
    Calculate the noise in a map, measured in uK-arcmin.
    The noise is calculated over a high enough ell range such that the
    cls should be noise-dominated.
    If two maps are provided, a difference map will be computed and that
    map will be used to calculate noise.
    Maps needs to be G3Frames.
    """

    if map2 is not None:
        map1 = subtract_two_maps(map1, map2, divide_by_two=True, in_place=True)
    # note here by default delta_l = 1

    cls = calculate_powerspectra(
        map1,
        lmin=300,
        lmax=8000,
        apod_mask=apod_mask,
        qu_eb="qu",
    )
    idx = np.where((cls["TT"].bin_centers > ell_range[0])
                   & (cls["TT"].bin_centers < ell_range[1]))[0]
    tt_noise = np.sqrt(np.mean((cls["TT"])[idx])) / (
        core.G3Units.arcmin * core.G3Units.uK
    )
    qq_noise = np.sqrt(np.mean((cls["QQ"])[idx])) / (
        core.G3Units.arcmin * core.G3Units.uK
    )
    uu_noise = np.sqrt(np.mean((cls["UU"])[idx])) / (
        core.G3Units.arcmin * core.G3Units.uK
    )
    if verbose:
        print(
            "The average noise between ells %0d and %0d is:"
            % (ell_range[0], ell_range[1])
        )
        print("TT: %.1f uK-arcmin" % tt_noise)
        print("QQ: %.1f uK-arcmin" % qq_noise)
        print("UU: %.1f uK-arcmin" % uu_noise)

    if return_all:
        return tt_noise, qq_noise, uu_noise, cls
    else:
        return tt_noise, qq_noise, uu_noise
