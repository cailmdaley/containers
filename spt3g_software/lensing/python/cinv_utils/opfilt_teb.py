import os
import numpy as np

from .. import utils
from .. import map_spec_utils
from ..map_spec_utils import calculate_teb, mult_weight, mult_map
from spt3g import core, coordinateutils

rad = core.G3Units.rad
arcmin = core.G3Units.arcmin
uK = core.G3Units.uK
K = core.G3Units.K

def cl2sinv(cl, clnk, tf2d, nft=0.0, nfp=0.0, ebeq=True, lmax=None):
    """
    Generate sinv filter from cl and sky noise

    Parameters:
    -----------
    cl: dictionary
        signal spectrum
    clnk: MapSpectraTEB object
        noise spectrum
    tf2d: MapSpectraTEB object
        2d transfer function
    nft: float
        noise-floor for Temperature  [uK-arcmin]
    nfp: float
        noise-floor for Polarization [uK-arcmin]

    Returns:
    -----------
    ClMatrixTEB: object
        Best estimate of sky signal, including:
        a) signal, b) "sky noise" diagonal in k, and
        c) subtracted pixel noise (white in k)
    """
    if (nft != 0.0) or (nfp != 0.0):
        clnk = clnk.copy()
        zs = np.zeros(clnk["T"].shape)

        clnk["T"] -= (nft * uK * arcmin) ** 2
        clnk["E"] -= (nfp * uK * arcmin) ** 2
        clnk["B"] -= (nfp * uK * arcmin) ** 2

        clnk["T"] = np.maximum(clnk["T"], zs)
        clnk["E"] = np.maximum(clnk["E"], zs)
        clnk["B"] = np.maximum(clnk["B"], zs)

    if ebeq == True:
        assert np.all(np.array(tf2d["E"]) == np.array(tf2d["B"]))
        clnk = clnk.copy()
        clnk["E"] = 0.5 * (clnk["E"] + clnk["B"])
        clnk["B"] = clnk["E"]

    ret = tf2d * tf2d * (tf2d * tf2d * cl + clnk).inverse()
    return ret.get_l_masked(lmax=lmax)


class SkyInverseFilter(map_spec_utils.ClMatrixTEB):
    """
    The simple sinv filter that is just the inverse of cl.
    Inherited from ClMatrixTEB class.
    """
    def __init__(self, cl):
        super(sinv_filt, self).__init__(cl)
        self.clmat = self.inverse().clmat
        self.hash = self.hashdict()

    def hashdict(self):
        return {"lmax": self.lmax, "clmat": self.clmat}


class NoiseInverseFilter(object):
    """
    The noise inverse filter object.

    Attributes
    -------
    transf: numpy array
        2D transfer function
    ninv: G3SkyMapWeights
        The inverse of the noise (weight)
    """
    def __init__(self, transf, ninv):
        self.transf = transf
        self.ninv = ninv

        self.hash = self.hashdict()

    def hashdict(self):
        ret = {}

        if type(self.transf) == np.ndarray:
            ret["transf"] = self.transf
        else:
            ret["transf"] = self.transf.hashdict()

        #ret["ninv"] = hashlib.sha1(self.ninv.TT.view(np.uint8)).hexdigest()

        return ret

    def mult_tqu(self, tqu):
        """ Multiply a tqumap by ninv, calculate the teb,
        and multiply by the transfer function
        """
        ret = calculate_teb(mult_weight( self.ninv, tqu)) * self.transf
        dx = ret["T"].dx
        dy = ret["T"].dy
        ret *= 1.0 / (dx * dy)
        return ret

    def mult_teb(self, teb):
        """Convert teb*transfer function into tqu maps,
        multiply the tqu maps by ninv,
        convert the result into teb, and multiply by the transfer function
        """
        #assert maps.pix.compatible(self.ninv, teb)
        ret = calculate_teb(mult_weight(self.ninv, (teb * self.transf).get_tqu())) * self.transf
        dx = ret["T"].dx
        dy = ret["T"].dy
        ret *= 1.0 / (dx * dy)
        return ret

    def get_lmax(self):
        if False:
            pass
        elif hasattr(self.transf, "get_ell"):
            return np.ceil(np.max(self.transf.get_ell().flatten()))
        elif hasattr(self.transf, "lmax"):
            return self.transf.lmax
        elif (getattr(self.transf, "size", 0) > 1) and (
            len(getattr(self.transf, "shape", ())) == 1
        ):
            return len(self.transf) - 1
        else:
            assert 0

    def get_fmask(self):
        if False:
            pass
        elif isinstance(self.ninv, coordinateutils.G3SkyMapWeights):
            mskt = np.flatnonzero(np.array(self.ninv.TT))
            mskq = np.flatnonzero(np.array(self.ninv.QQ))
            msku = np.flatnonzero(np.array(self.ninv.UU))
            ny, nx = np.shape(self.ninv.TT)
            if (not (np.all(mskt == mskq))) or (not (np.all(mskt == msku))):
                core.log_warn("WARNING: fmask t/q/u differ. defaulting to T.", unit="C-inverse")

            ret = np.zeros(nx * ny)
            ret[mskt] = 1.0
            return ret.reshape((ny, nx))
        elif (isinstance(self.ninv, dict) and "T" in self.ninv.keys() and
            "Q" in self.ninv.keys() and "U" in self.ninv.keys()):
            mskt = np.flatnonzero(np.array(self.ninv["T"]))
            mskq = np.flatnonzero(np.array(self.ninv["Q"]))
            msku = np.flatnonzero(np.array(self.ninv["U"]))

            assert np.all(np.array(mskt) == np.array(mskq))
            assert np.all(np.array(mskt) == np.array(msku))
            ny, nx = np.shape(ninv["T"])
            ret = np.zeros(nx * ny)
            ret[mskt] = 1.0
            return ret.reshape((ny, nx))
        else:
            return NotImplemented

    def get_fcut(self):
        fmask = self.get_fmask()
        return np.sum(fmask.flatten()) / fmask.size

    def get_fl(self):
        """The diagonal approximation of the ninv filter
        """
        if False:
            pass
        elif isinstance(self.ninv, coordinateutils.G3SkyMapWeights):
            lmax = int(self.get_lmax())
            fcut = self.get_fcut()
            dx = self.ninv.TT.x_res
            dy = self.ninv.TT.y_res
            ntt = np.average(self.ninv.TT) / fcut
            nee = (
                0.5
                * (
                    np.average(self.ninv.QQ)
                    + np.average(self.ninv.UU)
                )
                / fcut
            )

            ninv = map_spec_utils.ClMatrixTEB(
                {
                    "L": np.arange(lmax + 1),
                    "TT": ntt * np.ones(lmax + 1),
                    "EE": nee * np.ones(lmax + 1),
                    "BB": nee * np.ones(lmax + 1),
                }
            )

            return (
                ninv
                * self.transf
                * self.transf
                * (1.0 / (dx * dy))
            )
        elif (isinstance(self.ninv, dict) and "T" in self.ninv.keys() and
            "Q" in self.ninv.keys() and "U" in self.ninv.keys()):
            lmax = int(self.get_lmax())
            fcut = self.get_fcut()
            dx = ninv.TT.x_res
            dy = ninv.TT.y_res
            ntt = np.average(self.ninv["T"]) / fcut
            nee = (
                0.5
                * (np.average(self.ninv["Q"]) + np.average(self.ninv["U"]))
                / fcut
            )

            ninv = map_spec_utils.ClMatrixTEB(
                {
                    "L": np.arange(lmax + 1),
                    "TT": ntt * np.ones(lmax + 1),
                    "EE": nee * np.ones(lmax + 1),
                    "BB": nee * np.ones(lmax + 1),
                }
            )

            return (
                ninv
                * self.transf
                * self.transf
                * (1.0 / (dx * dy))
            )
        else:
            return NotImplemented


def calc_prep(tqu, sinv_filt, ninv_filt):
    return ninv_filt.mult_tqu(tqu)


def calc_fini(teb, sinv_filt, ninv_filt):
    return sinv_filt * teb


class DotOperator(object):
    """
    Multiply two dictionaries of TEB
    """
    def __init__(self, lmax=None):
        self.lmax = lmax

    def __call__(self, teb1, teb2):
        #assert teb1.compatible(teb2)

        if self.lmax != None:
            lmax = self.lmax

            return np.sum(
                (
                    teb1["T"] * np.conj(teb2["T"])
                    + teb1["E"] * np.conj(teb2["E"])
                    + teb1["B"] * np.conj(teb2["B"])
                )
                .flatten()[np.where(teb1.get_ell().flatten() <= lmax)]
                .real
            )

        else:
            return np.sum(
                (
                    teb1["T"] * np.conj(teb2["T"])
                    + teb1["E"] * np.conj(teb2["E"])
                    + teb1["B"] * np.conj(teb2["B"])
                )
                .flatten()
                .real
            )


class ForwardOperator(object):
    """ Returns [C^-1 + P^t N^-1 P] * teb,
    which I think is a 2D Fourier-space object
    """

    def __init__(self, sinv_filt, ninv_filt):
        self.sinv_filt = sinv_filt
        self.ninv_filt = ninv_filt

    def __call__(self, teb):
        return self.calc(teb)

    def calc(self, teb):
        return self.sinv_filt * teb + self.ninv_filt.mult_teb(teb)


class PreOperatorDiag(object):
    """
    (sinvfilt + ninv_filt's diagonal approximation) * teb
    """

    def __init__(self, sinv_filt, ninv_filt):
        self.filt = (sinv_filt + ninv_filt.get_fl()).inverse()

    def __call__(self, talm):
        return self.calc(talm)

    def calc(self, teb):
        return self.filt * teb
