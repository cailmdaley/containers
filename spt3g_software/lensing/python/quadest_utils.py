# spt3g_software/lensing/python/quadest_utils.py
# --
# this module contains classes and routines for applying quadratic
# anisotropy estimators to CMB maps.
#
# these estimators are motivated as follows: first, we start with
# a source of statistical anisotropy 's', which is a 2D field having
# a Harmonic transform s(L). 's' induces couplings between modes of the
# CMB temperature and polarization observables X, Y \in {T, E, B} as
#
# < X(l_X) Y(l_Y) > = \int{d^2 L} W^{s, XY}(l_X, l_Y, L) s(L),
#
# where <> represents an ensemble average over realizations of the CMB
# with s(L) fixed and W^{s, XY} is a weight function associated with the
# A-B coupling induced by 's'. we can try to measure s(L) using a
# quadratic estimator q, defined by
#
# q^{XY}(L) = 1/2 \int{d^2 l_X} \int{d^2 l_Y}
#                    W^{s, XY}(l_X, l_Y, L) \bar{X}(l_X) \bar{Y}(l_Y)
#
# where \bar{X} and \bar{Y} are filtered observations of the X and Y
# fields. nominally, the cost of evaluating q(L) is O(lmax^{6}) (because
# l_X, l_Y, and L are each 2D fields with lmax^2 elements), however if
# the weight function can be written as a sum of separable terms in
# l_X, l_Y, and L then this can be done much faster. for a flat-sky map
# in Fourier space, for example, the weight function could be written as
#
# W^{s, XY} = \sum_{i=0}^{N} \int{d^2 z}
#                   (e^{+i*2\pi*s^{i,X}+i*(l_X.z)} W^{i,X}(l_X)) *
#                    (e^{+i*2\pi*s^{i,Y}+i*(l_Y.z)} W^{i,Y}(l_Y)) *
#                     (e^{-i*2\pi*s^{i,L}+i*( -L.z)} W^{i,L}( L ))
#
# where s^{i,X/Y/L} are integers representing a spin parameter for the
# components of the weight function. with such a weight function q(L)
# may be evaluated in O(N_i * lmax^2 * log(lmax)) using fast Fourier
# transforms (FFTs). a similar result holds for full-sky maps, using
# fast Spherical Harmonic Transforms (SHTs) to obtain an algorithm
# which is slightly slower, at O(N_i * lmax^3) algorithm. separable
# weight functions and spins are encapsulated by the 'qest' class.
#
# in addition to q^{XY}(L), we often want to evaluate
#
#   * the response of the estimator q^{XY}(L) to s(L). if the
#     filtering which relates \bar{X} to X is diagonal in
#     Fourier space with \bar{X}(L) = F(L)X(L) then this
#     can also be calculated quickly with FFTs on the flat-sky,
#     or Wigner D-matrices on the full-sky. these calculations
#     are performed by the 'quadest_utils.fill_response' method.
#
#   * the ensemble-averaged cross-power <q1^{UV}(L) q2^{*XY}(L)>,
#     given estimates of the spectral cross-powers
#     <\bar{XU(l)\bar{X}^*(l)>, <\bar{V}(l),\bar{Y}^*(l)>, etc.
#     again, this cross-spectrum can be calculated quickly using
#     FFTs or D-matrices. this calculation is performed by the
#     'quadest_utils.fill_clqq' method.
#
#############################
# Classes for calculating quadratic estimates.
# In order to write the QE convolutions as FFT's
# we separate the equations into terms.
#   where each term can be expressed as components of the form
#        F(l) * G(L-l) * H(L)
#   Here, l is the integration variable.
# components self.wl[:][0] correspond to F(l)
# components self.wl[:][1] correspond to G(L-l)
# components self.wl[:][2] are independend of l, corresponding to H(L)
#
# self.npad_conv is the factor for padding the convolution-by-FFT calculations
#############################

import numpy as np
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum2D
import pdb
from spt3g import core

def pad_ft(a, npad=2):
    """ Pad an array in Fourier-space for FFT convolutions. """
    if npad == 1:
        return a
    nx, ny = a.shape
    p = np.zeros([nx * npad, ny * npad], dtype=a.dtype)
    p[0:nx, 0:ny] = np.fft.fftshift(a)
    p = np.roll(np.roll(p, -nx // 2, axis=0), -ny // 2, axis=1)
    return p


def unpad_ft(a, npad=2):
    """ Un-pad an array in Fourier-space. """
    if npad == 1:
        return a
    nx_pad, ny_pad = a.shape
    nx = int(nx_pad / npad)
    ny = int(ny_pad / npad)
    return np.roll(
        np.roll(
            (
                np.roll(np.roll(a, nx // 2, axis=0), ny // 2, axis=1)[
                    0:nx, 0:ny
                ]
            ),
            nx // 2,
            axis=0,
        ),
        ny // 2,
        axis=1,
    )


def iconvolve_padded(f, g, npad=2):
    """ Perform convolution of two harmonic functions.
    Calculate the convolution between harmonic functions f and g, defined as
       ret(L) = iint{d^2\vec{l} f(l) \times g(L-l)},
    by taking the product of the functions in real space.
    """
    return (
        unpad_ft(
            np.fft.fft2(
                np.fft.ifft2(pad_ft(f, npad=npad))
                * np.fft.ifft2(pad_ft(g, npad=npad))
            ),
            npad=npad,
        )
        * npad ** 2
    )


class QuadEst:
    """ Generic quadratic CMB anistropy estimator class.
    This class contains routines for calculating quadratic estimator
    coefficients, evaluating quadratic estimators on real-space maps/fields, and
    calculating the covariance between Fourier-transformed maps/fields.

    Attributes:
    ---------
    cl : dictionary
        Power spectrum dictionary holding e.g., TT, EE, TE spectra.
    npad : int
        Integer specifying padding for FFT operations.
    """

    def __init__(self, cl, npad=2):
        self.cl = cl
        self.npad = npad

        self.estimator_dict = {
            "Phi_TT": {"cls": (cl["TT"],), "fls": "TT"},
            "Phi_TT_s0": {"cls": (cl["TT"],), "fls": "TT"},
            "Phi_TT_curl": {"cls": (cl["TT"],), "fls": "TT"},
            "T_TPhi": {"cls": (cl["TT"], cl.get("PP", None)), "fls": "TT"},
            "Tau_TT": {"cls": (cl["TT"],), "fls": "TT"},
            "Phi_TE": {"cls": (cl["TE"],), "fls": "TE"},
            "Phi_ET": {"cls": (cl["TE"],), "fls": "ET"},
            "Phi_TE_curl": {"cls": (cl["TE"],), "fls": "TE"},
            "Phi_ET_curl": {"cls": (cl["TE"],), "fls": "ET"},
            "Phi_TB": {"cls": (cl["TE"],), "fls": "TB"},
            "Phi_BT": {"cls": (cl["TE"],), "fls": "BT"},
            "Phi_TB_curl": {"cls": (cl["TE"],), "fls": "TB"},
            "Phi_BT_curl": {"cls": (cl["TE"],), "fls": "BT"},
            "Tau_TB": {"cls": (cl["TE"],), "fls": "TB"},
            "Tau_TE": {"cls": (cl["TE"],), "fls": "TE"},
            "Phi_EE": {"cls": (cl["EE"],), "fls": "EE"},
            "Phi_EE_curl": {"cls": (cl["EE"],), "fls": "EE"},
            "Phi_EB": {"cls": (cl["EE"],), "fls": "EB"},
            "Phi_BE": {"cls": (cl["EE"],), "fls": "BE"},
            "Phi_EB_curl": {"cls": (cl["EE"],), "fls": "EB"},
            "Phi_BE_curl": {"cls": (cl["EE"],), "fls": "BE"},
            # not sure about fields for these
            "B_EPhi": {"cls": (cl["EE"], cl.get("PP", None))},
            "B_EPhi_neg": {"cls": (cl["EE"], cl.get("PP", None))},
            "B_EPhi_curl": {"cls": (cl["EE"], cl.get("PP", None))},
            "Noise_NoisePhi": {"cls": (cl["EE"], cl.get("PP", None))},
            "Tau_EB": {"cls": (cl["EE"],)},
            "Tau_EE": {"cls": (cl["EE"],)},
        }

    def calc_weights(self, estimator_id, lx, ly):
        """Calculates quadratic estimator coefficients.

        Given an estimator_id and l coordinates, return the correct wl and sl
        coefficients. Called by QuadEst.eval() and QuadEst.covariance().

        Attributes:
        --------
        estimator_id : str
            Unique ID that specifies both the quantity to be estimated and
            the quantities the estimate is based on. The quantity to be
            estimated comes first, and is separated from the quantities the
            estimate is based on by an underscore. For example, to estimate
            the potential Phi from T and B, one would use the estimator_id
            Phi_TB. In some cases, additional descriptors are added to the
            end separated by an underscore, e.g. Phi_TB_curl.
        lx, ly : ndarray
            1D array containing x- and y-axis l coordinates.
        """
        assert (
            estimator_id in self.estimator_dict
        ), "{} is not a valid estimator ID.".format(estimator_id)

        l = np.sqrt(lx ** 2 + ly ** 2)
        cls = self.estimator_dict[estimator_id]["cls"]
        interp_cls = [np.interp(l, np.arange(len(cl)), cl, right=0) for cl in cls]

        if estimator_id == "Phi_TT":
            # phi estimated from TT
            wl = np.array(
                [
                    [l * interp_cls[0], -0.5, l],
                    [l * interp_cls[0], -0.5, l],
                    [-0.5, l * interp_cls[0], l],
                    [-0.5, l * interp_cls[0], l],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+1, +0, +1], [-1, +0, -1], [+0, +1, +1], [+0, -1, -1]]
            )

        elif estimator_id == "Phi_TT_curl":
            # phi curl-mode estimated from TT

            wl, sl = self.calc_weights("Phi_TT", lx, ly)
            wl[:, 2] *= np.array([1.0j, -1.0j, 1.0j, -1.0j])

        elif estimator_id == "Phi_TT_s0":
            wl = np.array(
                [
                    [lx * interp_cls[0], +1, lx],
                    [ly * interp_cls[0], +1, ly],
                    [+1.0, lx * interp_cls[0], lx],
                    [+1.0, ly * interp_cls[0], ly],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+0, +0, +0], [+0, +0, +0], [+0, +0, +0], [+0, +0, +0]]
            )

        elif estimator_id == "Phi_TE":
            # phi estimated from TE

            wl = np.array(
                [
                    [l * interp_cls[0], -0.25, l],
                    [l * interp_cls[0], -0.25, l],
                    [l * interp_cls[0], -0.25, l],
                    [l * interp_cls[0], -0.25, l],
                    [-0.5, l * interp_cls[0], l],
                    [-0.5, l * interp_cls[0], l],
                ],
                dtype=object,
            )
            sl = np.array(
                [
                    [+3, -2, +1],
                    [-3, +2, -1],
                    [-1, +2, +1],
                    [+1, -2, -1],
                    [+0, +1, +1],
                    [+0, -1, -1],
                ]
            )

        elif estimator_id == "Phi_ET":
            # phi estimated from ET
            wl, sl = self.calc_weights("Phi_TE", lx, ly)
            wl[:, [0, 1]] = wl[:, [1, 0]]  # switch 0th and first columns
            sl[:, [0, 1]] = sl[:, [1, 0]]  # switch 0th and first columns

        elif estimator_id == "Phi_TE_curl":
            # phi curl-mode estimated from TE
            wl, sl = self.calc_weights("Phi_TE", lx, ly)
            wl[:, 2] *= np.array([1.0j, -1.0j, 1.0j, -1.0j, 1.0j, -1.0j])

        elif estimator_id == "Phi_ET_curl":
            # phi curl-mode estimated from ET
            wl, sl = self.calc_weights("Phi_TE_curl", lx, ly)
            wl[:, [0, 1]] = wl[:, [1, 0]]
            # switch 0th and first columns

        elif estimator_id == "Phi_TB":
            # phi estimated from TB
            wl = np.array(
                [
                    [l * interp_cls[0], +0.25j, l],
                    [l * interp_cls[0], -0.25j, l],
                    [l * interp_cls[0], -0.25j, l],
                    [l * interp_cls[0], +0.25j, l],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+3, -2, +1], [-3, +2, -1], [-1, +2, +1], [+1, -2, -1]]
            )

        elif estimator_id == "Phi_BT":
            # phi estimated from BT

            wl, sl = self.calc_weights("Phi_TB", lx, ly)
            # switch 0th and first columns
            wl[:, [0, 1]] = wl[:, [1, 0]]
            sl[:, [0, 1]] = sl[:, [1, 0]]

        elif estimator_id == "Phi_TB_curl":
            # Phi curl-mode estimated from TB
            wl, sl = self.calc_weights("Phi_TB", lx, ly)
            wl[:, 2] *= np.array([1.0j, -1.0j, 1.0j, -1.0j])

        elif estimator_id == "Phi_BT_curl":
            # phi curl-mode estimated from BT
            wl, sl = self.calc_weights("Phi_TB_curl", lx, ly)

        elif estimator_id == "Phi_EE":
            # phi estimated from EE
            wl = np.array(
                [
                    [-0.25, l * interp_cls[0], l],
                    [-0.25, l * interp_cls[0], l],
                    [l * interp_cls[0], -0.25, l],
                    [l * interp_cls[0], -0.25, l],
                    [-0.25, l * interp_cls[0], l],
                    [-0.25, l * interp_cls[0], l],
                    [l * interp_cls[0], -0.25, l],
                    [l * interp_cls[0], -0.25, l],
                ],
                dtype=object,
            )
            sl = np.array(
                [
                    [-2, +3, +1],
                    [+2, -3, -1],
                    [+3, -2, +1],
                    [-3, +2, -1],
                    [+2, -1, +1],
                    [-2, +1, -1],
                    [-1, +2, +1],
                    [+1, -2, -1],
                ]
            )

        elif estimator_id == "Phi_EE_curl":
            # phi curl-mode estimated from EE
            wl, sl = self.calc_weights("Phi_EE", lx, ly)
            wl[:, 2] *= np.array(
                [1.0j, -1.0j, 1.0j, -1.0j, 1.0j, -1.0j, 1.0j, -1.0j]
            )

        elif estimator_id == "Phi_EB":
            # phi estimated from EB
            wl = np.array(
                [
                    [l * interp_cls[0], +0.25j, l],
                    [l * interp_cls[0], -0.25j, l],
                    [l * interp_cls[0], -0.25j, l],
                    [l * interp_cls[0], +0.25j, l],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+3, -2, +1], [-3, +2, -1], [-1, +2, +1], [+1, -2, -1]]
            )

        elif estimator_id == "Phi_BE":
            # phi estimated from BE
            wl, sl = self.calc_weights("Phi_EB", lx, ly)
            wl[:, [0, 1]] = wl[:, [1, 0]]  # switch 0th and first columns
            sl[:, [0, 1]] = sl[:, [1, 0]]  # switch 0th and first columns

        elif estimator_id == "Phi_EB_curl":
            # phi curl-mode estimated from EB
            wl, sl = self.calc_weights("Phi_EB", lx, ly)
            wl[:, 2] *= np.array([1.0j, -1.0j, 1.0j, -1.0j])

        elif estimator_id == "Phi_BE_curl":
            # phi curl-mode estimated from BE
            wl, sl = self.calc_weights("Phi_EB_curl", lx, ly)
            wl[:, [0, 1]] = wl[:, [1, 0]]  # switch 0th and first columns

        elif estimator_id == "B_EPhi":
            wl = np.array(
                [
                    [l * interp_cls[0], l * interp_cls[1], -0.25j],
                    [l * interp_cls[0], l * interp_cls[1], +0.25j],
                    [l * interp_cls[0], l * interp_cls[1], +0.25j],
                    [l * interp_cls[0], l * interp_cls[1], -0.25j],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+3, -1, +2], [-3, +1, -2], [-1, -1, -2], [+1, +1, +2]]
            )

        elif estimator_id == "B_EPhi_neg":
            # the H(ell) wl coefficients in this case are -ve of those in
            # QuadEst_blm_EP for phi fields that are -ve of outputs from
            # cache_sim_ecp_Phi_map() e.g. those from healpy.alm2map()
            wl, sl = self.calc_weights("B_EPhi", lx, ly)
            wl[:, 2] *= -1.0

        elif estimator_id == "B_EPhi_curl":
            wl = np.array(
                [
                    [l * interp_cls[0], l * interp_cls[1], -0.25],
                    [l * interp_cls[0], l * interp_cls[1], -0.25],
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+3, -1, +2], [-3, +1, -2], [-1, -1, -2], [+1, +1, +2]]
            )

        elif estimator_id == "Noise_NoisePhi":
            # test class that calculates the effect on b-template if noise is
            # not filtered. blm_EP has the W(ell,ell;) term; this is the Z(ell,ell') term.
            wl = np.array(
                [
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                    [l * interp_cls[0], l * interp_cls[1], +0.25],
                ],
                dtype=object,
            )
            sl = np.array(
                [[+3, -1, +2], [-3, +1, -2], [-1, -1, -2], [+1, +1, +2]]
            )

        elif estimator_id == "T_TPhi":
            wl = np.array(
                [
                    [l * interp_cls[0], l * interp_cls[1], +1],
                    [l * interp_cls[0], l * interp_cls[1], +1],
                ],
                dtype=object,
            )
            sl = np.array([[+1, -1, +0], [-1, +1, +0]])

        elif estimator_id == "Tau_TT":
            wl = np.array(
                [[interp_cls[0], +1, -1], [+1, interp_cls[0], -1]],
                dtype=object,
            )
            sl = np.array([[+0, +0, +0], [+0, +0, +0]])

        elif estimator_id == "Tau_EB":
            wl = np.array(
                [[interp_cls[0], +1, -0.5j], [interp_cls[0], +1, -0.5j]],
                dtype=object,
            )
            sl = np.array([[-2, +2, +0], [+2, -2, +0]])

        elif estimator_id == "Tau_TB":
            wl = np.array(
                [[interp_cls[0], +1, -0.5j], [interp_cls[0], -1, -0.5j]],
                dtype=object,
            )
            sl = np.array([[-2, +2, +0], [+2, -2, +0]])

        elif estimator_id == "Tau_EE":
            wl = np.array(
                [
                    [interp_cls[0], +1, -0.5],
                    [interp_cls[0], +1, -0.5],
                    [+1, interp_cls[0], -0.5],
                    [+1, interp_cls[0], -0.5],
                ],
                dtype=object,
            )
            sl = np.array(
                [[-2, +2, +0], [+2, -2, +0], [-2, +2, +0], [+2, -2, +0]]
            )

        elif estimator_id == "Tau_TE":
            wl = np.array(
                [
                    [interp_cls[0], +1, -0.5],
                    [interp_cls[0], +1, -0.5],
                    [+1, interp_cls[0], -1],
                ],
                dtype=object,
            )
            sl = np.array([[-2, +2, +0], [+2, -2, +0], [+0, +0, +0]])

        return wl, sl

    def eval(self, ft1, ft2, estimator_id, npad=None):
        """
        Evaluate the chosen quadradic estimator from real-space maps ft1 and ft2.
        Note:  This function is written to return phi for estimators in which ft1=ft2, for example, TT.
            If Cl^{ft1, ft2} = 0 (e.g. the EB or EP estimators), then the correct filter must be multiplied by 2, for example:
              my_qest = 2. * quadest_utils.eval(ebar, bbar)
            See Equations 14 v.s. 15 in Hu & Okamoto 2002 (arXiv:0111606).

        Attributes:
        ---------
        estimator_id : str
            Unique ID that specifies both the quantity to be estimated and
            the quantities the estimate is based on. The quantity to be
            estimated comes first, and is separated from the quantities the
            estimate is based on by an underscore. For example, to estimate
            the potential Phi from T and B, one would use the estimator_id
            Phi_TB. In some cases, additional descriptors are added to the
            end separated by an underscore, e.g. Phi_TB_curl.
        ft1, ft2 : map objects
            Real-space maps on which the quadratic estimator is evaluated.
        npad: int
            Integer specifying padding for FFT operations.
        """

        npad = self.npad if npad is None else npad

        # assert r1.compatible(r2)
        # TODO: make eval and covariance consistent
        # in the way they preprocess that maps

        cfft = MapSpectrum2D(ft1.map_nx, ft1.dx, map_ny=ft1.map_ny, dy=ft1.dy)
        lx, ly = cfft.get_lxly()
        psi = np.arctan2(lx, -ly)

        wl, sl = self.calc_weights(estimator_id, lx, ly)
        for i in range(wl.shape[0]):  # loop over terms
            term1 = wl[i, 0] * ft1 * np.exp(+1.0j * sl[i, 0] * psi)
            term2 = wl[i, 1] * ft2 * np.exp(+1.0j * sl[i, 1] * psi)
            cfft += (
                iconvolve_padded(term1, term2, npad=npad)
                * wl[i, 2]
                * np.exp(-1.0j * sl[i, 2] * psi)
                * 0.5
                / np.sqrt(cfft.dx * cfft.dy)
                * np.sqrt(cfft.map_nx * cfft.map_ny)
            )
        return cfft

    def covariance(
        self,
        cfft,
        f1fft,
        f2fft,
        estimator_id1,
        estimator_id2,
        switch_12=False,
        switch_34=False,
        conj_12=False,
        conj_34=False,
        npad=None,
    ):
        """
        Calculate the covariance between two estimators with IDs `estimator_id1`
        and `estimator_id2` given Fourier-transformed fields/maps `f1fft` and `f2fft`.
        The resulting covariance is added in-place to `cfft.fft`, i.e. the
        complex-valued array containing the FFT associated with a `cfft` object.
        When used to calculate the response, the return value is half of the
        full response.
        Note: in the lensing quadratic estimator, it will always be the case
        that qe1 == qe2.

        Attributes:
        ------------
        estimator_id1, estimator_id2: str
            Unique IDs that specify both the quantity to be estimated and
            the quantities the estimate is based on. The quantity to be
            estimated comes first, and is separated from the quantities the
            estimate is based on by an underscore. For example, to estimate
            the potential Phi from T and B, one would use the estimator_id
            Phi_TB. In some cases, additional descriptors are added to the
            end separated by an underscore, e.g. Phi_TB_curl.
        cfft : complex FFT object
            The FFT object to which the covariance is added.
        f1fft, f2fft : ndarray
            The Fourier transformed map arrays used in the covariance
            calculation.
        switch_12, switch_34 : bool
            Booleans that control wether the first/second and third/fourth
            rows of the quadratic estimator coefficient matrix are switched.
        conj_12, conj_34 : bool
            Booleans that control wether to use the complex conjugates of
            the quadratic estimator coefficients.
        npad: int
            Integer specifying padding for FFT operations.
        """

        lx, ly = cfft.get_lxly()
        l = np.sqrt(lx ** 2 + ly ** 2)
        psi = np.arctan2(lx, -ly)
        nx, ny = l.shape

        if f1fft.shape != l.shape:
            assert f1fft.ndim == 1
            f1fft = np.interp(l.flatten(),
                np.arange(len(f1fft)), f1fft, right=0).reshape(l.shape)
        if f2fft.shape != l.shape:
            assert f2fft.ndim == 1
            f2fft = np.interp(l.flatten(), np.arange(len(f2fft)),
                f2fft, right=0).reshape(l.shape)
        # TODO: make eval and covariance consistent in the way they preprocess that maps

        npad = self.npad if npad is None else npad
        if npad != 2:
            core.log_warn(
                "lensing.quadest_utils.covariance():\
                npad is not equal to 2!  I hope you know what you are doing...",
                unit="Quadest_utils"
            )

        i0_12 = 0 if switch_12 is False else 1
        i0_34 = 0 if switch_34 is False else 1
        i1_12 = 1 if switch_12 is False else 0
        i1_34 = 1 if switch_34 is False else 0

        def conjfunc_12(x):
            if conj_12 is False:
                return x
            else:
                return np.conj(x)

        def conjfunc_34(x):
            if conj_34 is False:
                return x
            else:
                return np.conj(x)

        wl1, sl1 = self.calc_weights(estimator_id1, lx, ly)
        wl2, sl2 = self.calc_weights(estimator_id2, lx, ly)
        for i in range(0, wl1.shape[0]):  # loop over qe1 terms
            for j in range(0, wl2.shape[0]):  # loop over qe2 terms
                term1 = (
                    conjfunc_12(wl1[i, i0_12])
                    * conjfunc_34(wl2[j, i0_34])
                    * f1fft
                    * np.exp(
                        +1.0j
                        * (
                            (-1) ** conj_12 * sl1[i, i0_12]
                            + (-1) ** conj_34 * sl2[j, i0_34]
                        )
                        * psi
                    )
                )
                term2 = (
                    conjfunc_12(wl1[i, i1_12])
                    * conjfunc_34(wl2[j, i1_34])
                    * f2fft
                    * np.exp(
                        +1.0j
                        * (
                            (-1) ** conj_12 * sl1[i, i1_12]
                            + (-1) ** conj_34 * sl2[j, i1_34]
                        )
                        * psi
                    )
                )
                cfft += (
                    iconvolve_padded(term1, term2, npad=npad)
                    * conjfunc_12(wl1[i, 2])
                    * conjfunc_34(wl2[j, 2])
                    * np.exp(
                        -1.0j
                        * (
                            (-1) ** conj_12 * sl1[i, 2]
                            + (-1) ** conj_34 * sl2[j, 2]
                        )
                        * psi
                    )
                    * 0.25
                    / (cfft.dx * cfft.dy)
                )
        return cfft

    def fill_response(
        self, cfft, f1fft, f2fft, estimator_id1, estimator_id2, npad=2
    ):
        """ Wrapper around QuadEst.covariance() for calculating response.

        Calculate the response between two estimators with IDs `estimator_id1`
        and `estimator_id2` given Fourier-transformed fields/maps `f1fft`
        and `f2fft`.
        The resulting response overwrites `cfft`, i.e. the complex-valued array
        containing the FFT associated with a `cfft` object.
        """

        core.log_notice("fill_response", estimator_id1, estimator_id2, unit="Quadest_utils")
        cfft[:, :] = 0.0
        cfft = self.covariance(
            cfft, f1fft, f2fft, estimator_id1, estimator_id2, npad=npad
        )
        cfft *= 2.0  # because covariance returns 1/2 the response
        return cfft

    def fill_clqq(
        self, cfft, estimator_id1, estimator_id2, f11, f12, f22, npad=2
    ):
        cfft[:, :] = 0.0
        self.covariance(
            cfft,
            estimator_id1,
            estimator_id2,
            f11,
            f22,
            switch_34=False,
            conj_34=True,
            npad=npad,
        )
        self.covariance(
            cfft,
            estimator_id1,
            estimator_id2,
            f12,
            f12,
            switch_34=True,
            conj_34=True,
            npad=npad,
        )
        return cfft
