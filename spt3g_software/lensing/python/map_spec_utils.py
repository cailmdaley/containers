import os, copy, glob, hashlib, warnings
import numpy as np
from . import utils
from spt3g import core, mapmaker, coordinateutils
from spt3g.mapspectra import basicmaputils, map_analysis, apodmask
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum2D, MapSpectrum1D
from spt3g.simulations import cmb

rad = core.G3Units.rad
arcmin = core.G3Units.arcmin
uK = core.G3Units.uK
K = core.G3Units.K


class MapSpectraTEB(dict):
    """
    This is a dictionary wrapper of three TEB MapSpectrum2D objects
    """
    def __init__(
        self,
        ffts=None,
        map_nx=None,
        map_ny=None,
        dx=None,
        dy=None,
        units=getattr(core.G3TimestreamUnits, "None"),
        proj=getattr(coordinateutils.MapProjection, "ProjNone"),
    ):
        """ class which contains the FFTs of TEB.
        temperature (T), E and B-mode polarization.
        Initialize the class with TEB MapSpectrum2Ds or numpy arrays
        """
        if ffts is None:
            spec = np.zeros((map_ny, map_nx // 2 + 1), dtype=np.complex)
            ffts = {"T": spec, "E": spec, "B": spec}
        # can initialize with a list of TEB ffts
        if isinstance(ffts, list):
            ffts = {"T": ffts[0], "E": ffts[1], "B": ffts[2]}
        # can initialize with a dictionary of TEB ffts
        for k, fft in ffts.items():
            if isinstance(fft, MapSpectrum2D):
                pass
            elif isinstance(fft, np.ndarray):
                ffts[k] = MapSpectrum2D(
                    map_nx,
                    dx,
                    spec=np.array(fft, dtype=np.complex),
                    map_ny=map_ny,
                    dy=dy,
                    units=units,
                    proj=proj,
                )
            else:
                raise TypeError(
                    "Argument {} is not a MapSpectrum2D or dict".format(k)
                )
        super().__init__(ffts)

    def hashdict(self):
        return {
            "pix": {
                "map_nx": self["T"].map_nx,
                "dx": self["T"].dx,
                "map_ny": self["T"].map_ny,
                "dy": self["T"].dy,
            },
            "T": hashlib.sha1(self["T"].view(np.uint8)).hexdigest(),
            "E": hashlib.sha1(self["E"].view(np.uint8)).hexdigest(),
            "B": hashlib.sha1(self["B"].view(np.uint8)).hexdigest(),
        }

    def copy(self):
        return MapSpectraTEB({k: self[k].copy() for k in self.keys()})

    def __add__(self, other):
        ret = self.copy()
        ret += other
        return ret

    def __iadd__(self, other):
        if isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] += other[spec]
            return self
        elif isinstance(other, MapSpectrum2D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] += other
            return self
        else:
            return NotImplemented

    def __sub__(self, other):
        ret = self.copy()
        ret -= other
        return ret

    def __isub__(self, other):
        if isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] -= other[spec]
            return self
        elif isinstance(other, MapSpectrum2D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] -= other
            return self
        else:
            return NotImplemented

    def __mul__(self, other):
        ret = self.copy()
        ret *= other
        return ret

    def __imul__(self, other):
        if np.isscalar(other) or (
            (type(other) == np.ndarray)
            and (getattr(other, "shape", None) == self["T"].shape)
        ):
            for spec in self.keys():
                self[spec] *= other
            return self
        elif (type(other) == np.ndarray) and (
            len(getattr(other, "shape", [])) == 1
        ):
            map_ny, map_nx = self["T"].map_ny, self["T"].map_nx
            tfac = np.interp(
                self.get_ell().flatten(),
                np.arange(0, len(other)),
                other,
                right=0,
            ).reshape((map_ny, map_nx // 2 + 1))
            for spec in self.keys():
                self[spec] *= tfac
            return self

        elif isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] *= other[spec]
            return self
        elif isinstance(other, MapSpectrum2D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] *= other
            return self
        # to do: finish the next two options
        elif isinstance(other, ClMatrixTEB):
            return other * self
        elif is_cambcl(other):
            return ClMatrixTEB(other) * self
        else:
            return NotImplemented

    def __truediv__(self, other):
        ret = self.copy()
        ret /= other
        return ret

    def __itruediv__(self, other):
        if isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] /= other[spec]
            return self
        elif isinstance(other, MapSpectrum2D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] /= other
                self[spec] = np.nan_to_num(self[spec])
            return self
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if np.isscalar(other):
            ret = self.copy()
            np.asarray(ret["T"])[:] = 0
            np.asarray(ret["E"])[:] = 0
            np.asarray(ret["B"])[:] = 0
            ret["T"][self["T"] != 0] = other / self["T"][self["T"] != 0]
            ret["E"][self["E"] != 0] = other / self["E"][self["E"] != 0]
            ret["B"][self["B"] != 0] = other / self["B"][self["B"] != 0]
            return ret
        else:
            return NotImplemented

    def inverse(self):
        """
        return a new object for which all elements have been set to
        their inverses, with exception of zeros which are untouched.
        """
        return 1.0 / self

    def degrade(self, fac):
        """
        reduce the resolution of this map by a factor fac.
        """
        assert np.mod(self["T"].map_nx, fac) == 0
        assert np.mod(self["T"].map_ny, fac) == 0
        assert np.mod(self["T"].map_nx // fac, 2) == 0
        map_nx = self["T"].map_nx
        map_ny = self["T"].map_ny
        return self.__class__(
            ffts=[
                np.array(
                    self["T"][0 : map_ny // fac, 0 : map_nx // fac // 2 + 1]
                ),
                np.array(
                    self["E"][0 : map_ny // fac, 0 : map_nx // fac // 2 + 1]
                ),
                np.array(
                    self["B"][0 : map_ny // fac, 0 : map_nx // fac // 2 + 1]
                ),
            ],
            map_nx=map_nx // fac,
            map_ny=map_ny // fac,
            dx=self["T"].dx * fac,
            dy=self["T"].dy * fac,
            units=self["T"].units,
            proj=self["T"].proj,
        )

    def get_pixel_window(self):
        """
        return the spectrum describing the map-level transfer function
        for the pixelization of this object.
        """
        return self["T"].get_pixel_window()

    def get_cl(self, lbins):
        """
        returns a MapSpectra object containing the power-spectra of
        T,E,B in this map.
        """
        teb_complex = self.get_complex()
        return map_analysis.calculate_powerspectra(teb_complex, lbins)

    def get_tqu(self):
        """
        returns the tqumap given by the inverse Fourier transform
        of this object.
        """
        rffts = self.get_real()
        lx, ly = rffts.get_lxly()
        tpi = 2.0 * np.arctan2(lx, -ly)
        tfac = np.sqrt(
            (rffts["T"].map_nx * rffts["T"].map_ny)
            / (rffts["T"].dx * rffts["T"].dy)
        )
        tmap = np.fft.irfft2(rffts["T"]) * tfac
        qmap = (
            np.fft.irfft2(np.cos(tpi) * rffts["E"] - np.sin(tpi) * rffts["B"])
            * tfac
        )
        umap = (
            np.fft.irfft2(np.sin(tpi) * rffts["E"] + np.cos(tpi) * rffts["B"])
            * tfac
        )
        tqumap = {}
        t_map = coordinateutils.FlatSkyMap(
            tmap,
            rffts["T"].dx,
            is_weighted=False,
            units=rffts["T"].units,
            pol_type=coordinateutils.MapPolType.T,
        )
        q_map = coordinateutils.FlatSkyMap(
            qmap,
            rffts["T"].dx,
            is_weighted=False,
            units=rffts["T"].units,
            pol_type=coordinateutils.MapPolType.Q,
        )
        u_map = coordinateutils.FlatSkyMap(
            umap,
            rffts["T"].dx,
            is_weighted=False,
            units=rffts["T"].units,
            pol_type=coordinateutils.MapPolType.U,
        )
        tqumap["T"] = t_map
        tqumap["Q"] = q_map
        tqumap["U"] = u_map
        return tqumap

    def get_real(self):
        """
        returns a (real) MapSpectraTEB object.
        """
        if self["T"].is_real and self["T"].is_real and self["T"].is_real:
            return self
        else:
            ffts = {}
            for key in self.keys():
                if not self[key].is_real:
                    ffts[key] = self[key][:, 0 : self["T"].map_nx // 2 + 1]
                    ffts[key].is_real = True
                else:
                    ffts[key] = self[key]
            return self.__class__(ffts=ffts)

    def get_complex(self):
        """
        returns a (complex) MapSpectraTEB object.
        """
        if (
            (not self["T"].is_real)
            and (not self["T"].is_real)
            and (not self["T"].is_real)
        ):
            return self
        else:
            return self.__class__(
                ffts={
                    "T": self["T"].get_complex(),
                    "E": self["E"].get_complex(),
                    "B": self["B"].get_complex(),
                }
            )

    def get_l_masked(
        self,
        lmin=None,
        lmax=None,
        lxmin=None,
        lxmax=None,
        lymin=None,
        lymax=None,
    ):
        """ returns a copy of this object which has been masked
            zero in a customizable range of Fourier space.
        """

        return self.__class__(
            ffts={
                "T": self["T"].get_l_masked(
                    lmin=lmin,
                    lmax=lmax,
                    lxmin=lxmin,
                    lxmax=lxmax,
                    lymin=lymin,
                    lymax=lymax,
                ),
                "E": self["E"].get_l_masked(
                    lmin=lmin,
                    lmax=lmax,
                    lxmin=lxmin,
                    lxmax=lxmax,
                    lymin=lymin,
                    lymax=lymax,
                ),
                "B": self["B"].get_l_masked(
                    lmin=lmin,
                    lmax=lmax,
                    lxmin=lxmin,
                    lxmax=lxmax,
                    lymin=lymin,
                    lymax=lymax,
                ),
            }
        )

    def get_l_mask(
        self,
        lmin=None,
        lmax=None,
        lxmin=None,
        lxmax=None,
        lymin=None,
        lymax=None,
    ):
        """ return a Fourier mask for the pixelization associated with
            this object which is zero over customizable ranges of L.
        """
        mask = self["T"].get_l_mask(
            lmin=lmin,
            lmax=lmax,
            lxmin=lxmin,
            lxmax=lxmax,
            lymin=lymin,
            lymax=lymax,
        )
        return self.__class__(
            ffts=[mask, mask, mask],
            map_nx=self["T"].map_nx,
            map_ny=self["T"].map_ny,
            dx=self["T"].dx,
            dy=self["T"].dy,
            units=self["T"].units,
            proj=self["T"].proj,
        )

    def get_lxly(self):
        """ returns the (lx, ly) pair associated with
            each Fourier mode in T, E, B.
        """
        return self["T"].get_lxly()

    def get_ell(self):
        """ returns the wavenumber l = sqrt(lx**2 + ly**2)
            for each Fourier mode in T, E, B.
        """
        return self["T"].get_ell()


def mult_weight(weight, map):
    """ Multiply weights to a map.
    Weight here can be generated later and be an inverse-noise covariance matrix
    """
    t = (
        map["T"] * weight.TT
        + map["Q"] * weight.TQ
        + map["U"] * weight.TU
    )
    q = (
        map["T"] * weight.TQ
        + map["Q"] * weight.QQ
        + map["U"] * weight.QU
    )
    u = (
        map["T"] * weight.TU
        + map["Q"] * weight.QU
        + map["U"] * weight.UU
    )
    return {"T":t, "Q":q, "U":u}


def mult_map(map1, map2):
    if isinstance(map1, np.ndarray):
        map1 = {"T": map1, "Q": map1, "U": map1}
    if isinstance(map2, np.ndarray):
        map2 = {"T": map2, "Q": map2, "U": map2}
    t = map1["T"] * map2["T"]
    q = map1["Q"] * map2["Q"]
    u = map1["U"] * map2["U"]
    return {"T": t, "Q": q, "U": u}


def calculate_teb(
    map_fr,
    apod_mask=None,
    padded_map_shape=None,
    e_mode_method="basic",
    b_mode_method="basic",
    qu_eb="eb",
):
    """
    Calculates the ffts from an input map frame.

    Parameters
    ----------
    map_fr: Map frame or a dictionary
       A Map frame or map dic with weight removed and maps flattened.
    apod_mask [None]: FlatSkyMap, 2d array, file path.
        Apodization mask. Could a filepath to a .g3 or .pkl of such
        If left as None, no masking is done.
    padded_map_shape [None]: 2d array [npix_y, npix_x] or 'square'
        Specify the shape of the fft map. Pad maps to this size.
        If set to square, the larger dimension of the map is used
        for both npix_x, npix_y.
    e_mode_method ['basic']: string
        How to calculate the E mode. Options are "chi", "smith",
        and "basic". See constructEB docstring
    b_mode_method ['basic']: string
        How to calculate the B mode. Options are "chi", "smith",
        and "basic". See constructEB docstring
    Returns
    -------
    MapSpectraTEB object
    """

    # Set up pixel mask
    if isinstance(apod_mask, str):
        if apod_mask.endswith(".pkl") or apod_mask.endswith(".pk"):
            apod_mask = files.load_pickle(apod_mask)
        elif apod_mask.endswith(".g3"):
            g3file = core.G3File(apod_mask)
            for frame in g3file:
                if frame.type == core.G3FrameType.Map:
                    apod_mask = frame["T"]
        else:
            raise ValueError("Cannot read apod mask file %s" % apod_mask)

    # for polarized map frame
    if (
        ("T" in map_fr.keys())
        and ("Q" in map_fr.keys())
        and ("U" in map_fr.keys())
    ):

        tfft = basicmaputils.map_to_ft(
            map_fr["T"],
            apod_mask,
            res=map_fr["T"].res,
            padded_map_shape=padded_map_shape,
        )

        ffts = map_analysis.constructEB(
            map_fr,
            res=map_fr["T"].res,
            apod_mask=apod_mask,
            e_mode_method=e_mode_method,
            b_mode_method=b_mode_method,
            padded_map_shape=padded_map_shape,
        )

        ffts.update({"T": tfft})
        return MapSpectraTEB(ffts=ffts).get_real()

    else:
        raise ValueError("The input does not have TQU.")


def make_tqumap_wt(
    map_nx,
    map_ny,
    dx,
    dy=None,
    ninv=None,
    mask=None,
    ninv_dcut=None,
    nlev_tp=None,
    maskt=None,
    maskq=None,
    masku=None,
    verbose=True,
):
    """ function to generate a tqumap_wt which describes
    an inverse-noise covariance matrix.

    Parameters
    -----------
    map_nx: int
        number of pixels in x direction
    map_ny: int
        number of pixels in y direction
    dx: float
        resolution in x direction
    dy: float
        resolution in y direction
    ninv: G3SkyMapWeights.
        pixels for which this matrix weight function has determinant
        < ninv_dcut will be masked.
    mask: np.ndarray
        global mask map to apply
        (taking noise level to infinity for pixels where mask is zero).
    ninv_dcut: float
        Used only in conjunction with ninv.
        pixels with det smaller than ninv_dcut will be cut.
    nlev_tp   = a tuple (nT, nP), in uK arcmin
        giving pixel temperature/polarization white noise levels
        to use for the noise covariance in uK.arcmin.
        If !=None, TQ,TU,QU weights are set to zero.
    maskt, maskq, masku: np.ndarray
        Optional, individual T, Q and U masks to apply.
    """
    if dy == None:
        dy = dx
    attrs = ["TT", "TQ", "TU", "QU", "QQ", "UU"]

    # if ninv is given, only mask elements with determinant out of range
    if ninv is not None:
        dets = mapmaker.mapmakerutils.make_determinant_map(ninv)
        # fixthis
        #import pickle as pk

        #dets = pk.load(open("/home/panz/test.pkl", "rb"))
        ninv = ninv.Clone(False)
    # else, generate the cov matrix from the white noise level
    else:
        ninv = coordinateutils.G3SkyMapWeights()

    if nlev_tp is not None:
        for attr in attrs:
            weight = np.zeros((map_ny, map_nx))
            if attr == "TT":
                weight[:, :] = dx * dy / (nlev_tp[0] * uK * arcmin) ** 2
            elif attr == "QQ" or attr == "UU":
                weight[:, :] = dx * dy / (nlev_tp[1] * uK * arcmin) ** 2
            weight = coordinateutils.FlatSkyMap(weight, dx)
            setattr(ninv, attr, weight)

    if ninv_dcut is not None:
        for attr in attrs:
            if verbose:
                core.log_notice(
                    "cutting %d pixels for det"
                    % len(np.where(dets < ninv_dcut)[0])
                , unit="Map_spec_utils")
            weight_cut = np.array(getattr(ninv, attr))
            weight_cut[np.where(dets < ninv_dcut)] = 0.0
            weight_cut = coordinateutils.FlatSkyMap(weight_cut, dx)
            setattr(ninv, attr, weight_cut)

    if mask is not None:
        for attr in attrs:
            temp = getattr(ninv, attr)
            np.asarray(temp)[:] = np.asarray(temp)[:] * mask
            setattr(ninv, attr, temp)
    if maskt is not None:
        # multiply TT, TQ, TU weights by maskt
        for attr in attrs:
            if "T" in attr:
                setattr(ninv, attr, getattr(ninv, attr) * maskt)
    if maskq is not None:
        for attr in attrs:
            if "Q" in attr:
                setattr(ninv, attr, getattr(ninv, attr) * maskq)
    if masku is not None:
        for attr in attrs:
            if "U" in attr:
                setattr(ninv, attr, getattr(ninv, attr) * masku)

    return ninv


    def __mul__(self, other):
        if np.isscalar(other):
            return ClMatrixT(self.clmat * other)
        elif (getattr(other, "size", 0) > 1) and (
            len(getattr(other, "shape", ())) == 1
        ):
            assert (self.lmax + 1) <= len(other)
            return ClMatrixT(self.clmat * other[0 : lmax + 1])

        elif isinstance(other, MapSpectrum2D):
            ret = other.copy()
            ell = other.get_ell()

            def fftxcl(fft, cl):
                return fft * np.interp(
                    ell.flatten(), np.arange(0, len(cl)), cl, right=0
                ).reshape(fft.shape)

            ret.fft[:, :] = fftxcl(other.fft, self.clmat[:])
            return ret
        else:
            assert 0

    def inverse(self):
        ret = ClMatrixT(np.zeros(self.lmax + 1))
        ret.clmat[np.nonzero(self.clmat)] = (
            1.0 / self.clmat[np.nonzero(self.clmat)]
        )
        return ret

    def cholesky(self):
        return ClMatrixT(np.sqrt(self.clmat))


class ClMatrixTEB(object):
    """Class to hold the 3x3 covariance matrix at each multipole
    for a set of T, E and B auto- and cross-spectra.
    """

    def __init__(self, cl):
        """ Initializes this ClMatrixTEB object using the power spectra cl
        Spectra which are not present in cl are assumed to be zero.
        """
        lmax = int(np.max(cl["L"]))
        zs = np.zeros(lmax + 1)

        clmat = np.zeros(
            (lmax + 1, 3, 3)
        )  # matrix of TEB correlations at each l.
        clmat[:, 0, 0] = cl.get("TT", zs.copy())
        clmat[:, 0, 1] = cl.get("TE", zs.copy())
        clmat[:, 1, 0] = clmat[:, 0, 1]
        clmat[:, 0, 2] = cl.get("TB", zs.copy())
        clmat[:, 2, 0] = clmat[:, 0, 2]
        clmat[:, 1, 1] = cl.get("EE", zs.copy())
        clmat[:, 1, 2] = cl.get("EB", zs.copy())
        clmat[:, 2, 1] = clmat[:, 1, 2]
        clmat[:, 2, 2] = cl.get("BB", zs.copy())

        self.lmax = lmax
        self.clmat = clmat

    def hashdict(self):
        return {
            "lmax": self.lmax,
            "clmat": hashlib.md5(self.clmat.view(np.uint8)).hexdigest(),
        }

    def compatible(self, other):
        """Test whether this object and the ClMatrixTEB object other
        can be added, subtracted, or multiplied.
        """
        return (self.lmax == other.lmax) and (
            self.clmat.shape == other.clmat.shape
        )

    def clone(self, lmax=None):
        if lmax == None:
            lmax = self.lmax
        ret = ClMatrixTEB({"L": np.zeros(lmax + 1)})
        ret.clmat[:, :, :] = self.clmat[0 : lmax + 1, :, :]

        return ret

    def __add__(self, other):
        if isinstance(other, ClMatrixTEB):
            assert self.compatible(other)
            ret = copy.deepcopy(self)
            ret.clmat += other.clmat
            return ret
        elif isinstance(other, MapSpectraTEB):
            teb = other
            ret = teb.copy()
            ell = teb.get_ell()

            ret["T"] += np.interp(
                ell.flatten(),
                np.arange(0, len(self.clmat[:, 0, 0])),
                self.clmat[:, 0, 0],
                right=0,
            ).reshape(ell.shape)
            ret["E"] += np.interp(
                ell.flatten(),
                np.arange(0, len(self.clmat[:, 1, 1])),
                self.clmat[:, 1, 1],
                right=0,
            ).reshape(ell.shape)
            ret["B"] += np.interp(
                ell.flatten(),
                np.arange(0, len(self.clmat[:, 2, 2])),
                self.clmat[:, 2, 2],
                right=0,
            ).reshape(ell.shape)

            return ret
        else:
            return NotImplemented

    def __mul__(self, other):
        if np.isscalar(other):
            ret = self.clone()
            ret.clmat *= other
            return ret
        elif isinstance(other, ClMatrixTEB):
            assert self.compatible(other)
            ret = self.clone()
            ret.clmat *= other.clmat
            return ret
        elif (getattr(other, "size", 0) > 1) and (
            len(getattr(other, "shape", ())) == 1
        ):
            lmax = self.lmax
            assert (lmax + 1) <= len(other)

            ret = self.clone()
            for i in range(0, 3):
                for j in range(0, 3):
                    ret.clmat[:, i, j] *= other[0 : lmax + 1]

            return ret

        elif isinstance(other, MapSpectraTEB):
            teb = other
            ret = teb.copy()
            ell = teb.get_ell()

            def fftxcl(fft, cl):
                return fft * np.interp(
                    ell.flatten(), np.arange(0, len(cl)), cl, right=0
                ).reshape(fft.shape)

            ret["T"] = (
                fftxcl(teb["T"], self.clmat[:, 0, 0])
                + fftxcl(teb["E"], self.clmat[:, 0, 1])
                + fftxcl(teb["B"], self.clmat[:, 0, 2])
            )
            ret["E"] = (
                fftxcl(teb["T"], self.clmat[:, 1, 0])
                + fftxcl(teb["E"], self.clmat[:, 1, 1])
                + fftxcl(teb["B"], self.clmat[:, 1, 2])
            )
            ret["B"] = (
                fftxcl(teb["T"], self.clmat[:, 2, 0])
                + fftxcl(teb["E"], self.clmat[:, 2, 1])
                + fftxcl(teb["B"], self.clmat[:, 2, 2])
            )

            return ret
        else:
            return NotImplemented

    def inverse(self):
        """ Return a new ClMatrixTEB object, which contains the matrix inverse
        of this one, multipole-by-multipole. """
        ret = copy.deepcopy(self)
        for l in range(0, self.lmax + 1):
            ret.clmat[l, :, :] = np.linalg.pinv(self.clmat[l])
        return ret

    def cholesky(self):
        """ Return a new ClMatrixTEB object, which contains the cholesky
        decomposition (or matrix square root) of this one,
        multipole-by-multipole. """
        ret = copy.deepcopy(self)
        for l in range(0, self.lmax + 1):
            u, t, v = np.linalg.svd(self.clmat[l])
            ret.clmat[l, :, :] = np.dot(u, np.dot(np.diag(np.sqrt(t)), v))
        return ret


def cl2tebfft(cl, map_nx, dx, map_ny=None, dy=None):
    """
    Interpolate TEB Cls into the 2D ell grid
    """
    ell = (
        MapSpectrum2D(map_nx, dx, map_ny=map_ny, dy=dy)
        .get_real()
        .get_ell()
        .flatten()
    )
    return MapSpectraTEB(
        ffts=[
            np.array(
                np.interp(ell, cl["L"], cl["TT"], right=0).reshape(
                    map_ny, map_nx // 2 + 1
                ),
                dtype=np.complex,
            ),
            np.array(
                np.interp(ell, cl["L"], cl["EE"], right=0).reshape(
                    map_ny, map_nx // 2 + 1
                ),
                dtype=np.complex,
            ),
            np.array(
                np.interp(ell, cl["L"], cl["BB"], right=0).reshape(
                    map_ny, map_nx // 2 + 1
                ),
                dtype=np.complex,
            ),
        ],
        map_nx=map_nx,
        map_ny=map_ny,
        dx=dx,
        dy=dy,
    )


def get_camb_scalcl(prefix=None, lmax=np.inf, lmin=0, fname=None, g3unit=True):
    """
    loads and returns a "scalar Cls" file produced by CAMB (camb.info).
    can either use a prefix indicating one of the Cl files included
    in lensing/data, or a full filename.
    lmax sets maximum multipole to load
    (all multipoles will be loaded by default).
    """
    if fname == None:
        basedir = os.path.dirname(cmb.__file__)
        if prefix == None:
            prefix = "planck15"
        fname = basedir + "/data/camb/" + prefix + "/*_scalCls.dat"
    tf = glob.glob(fname)
    assert len(tf) == 1
    cls = cmb.read_camb(tf[0], as_cls=True, lmax=lmax, lmin=lmin)
    if g3unit:
        for key in ["TT", "EE", "TE", "PP", "TP", "BB", "EP", "EB"]:
            if key in cls.keys():
                if "P" not in key:
                    cls[key] = cls[key] * uK ** 2
                elif "PP" in key:
                    pass
                elif "P" in key:
                    cls[key] = cls[key] * uK
    return cls


def get_camb_lensedcl(
    prefix=None, lmax=np.inf, lmin=0, fname=None, g3unit=True
):
    """
    loads and returns a "lensed Cls" file produced by CAMB (camb.info).
    can either use a prefix indicating one of the Cl files included
    in lensing/data, or a full filename.
    lmax sets maximum multipole to load
    (all multipoles will be loaded by default). """
    if fname == None:
        basedir = os.path.dirname(cmb.__file__)
        if prefix == None:
            prefix = "planck15"
        fname = basedir + "/data/camb/" + prefix + "/*_lensedCls.dat"
    tf = glob.glob(fname)
    assert len(tf) == 1
    cls = cmb.read_camb(tf[0], as_cls=True, lmax=lmax, lmin=lmin)
    if g3unit:
        for key in ["TT", "EE", "TE", "PP", "TP", "BB", "EP", "EB"]:
            if key in cls.keys():
                if "P" not in key:
                    cls[key] = cls[key] * uK ** 2
                elif "PP"in key:
                    pass
                elif "P" in key:
                    cls[key] = cls[key] * uK
    return cls


def get_camb_lenspotentialcl(
    prefix=None, lmax=np.inf, lmin=0, fname=None, g3unit=True
):
    """
    loads and returns a "lenspotential Cls" file produced by CAMB
    (camb.info).
    can either use a prefix indicating one of the Cl files included in
    lensing/data, or a full filename.
    lmax sets maximum multipole to load
    (all multipoles will be loaded by default). """
    if fname == None:
        basedir = os.path.dirname(cmb.__file__)
        if prefix == None:
            prefix = "planck15"
        fname = basedir + "/data/camb/" + prefix + "/*_lenspotentialCls.dat"
    tf = glob.glob(fname)
    assert len(tf) == 1
    cls = cmb.read_camb(tf[0], as_cls=True, lmax=lmax, lmin=lmin)
    if g3unit:
        for key in ["TT", "EE", "TE", "PP", "TP", "BB", "EP", "EB"]:
            if key in cls.keys():
                if "P" not in key:
                    cls[key] = cls[key] * uK ** 2
                elif "PP"in key:
                    pass
                elif "P" in key:
                    cls[key] = cls[key] * uK
    return cls


def is_cambcl(input):
    if not isinstance(input, dict):
        return False
    if "L" not in input.keys():
        return False
    return set(input.keys()).issubset(
        set(["L", "TT", "EE", "TE", "PP", "TP", "BB", "EP", "EB"])
    )
