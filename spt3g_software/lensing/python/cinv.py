import os
import sys
import numpy as np
import pickle as pk
import hashlib
from spt3g import core
from . import utils
from . import map_spec_utils
from . import cinv_utils
from .map_spec_utils import calculate_teb, MapSpectraTEB, mult_map

rad = core.G3Units.rad
arcmin = core.G3Units.arcmin
uK = core.G3Units.uK
K = core.G3Units.K

class CinvFiltFFT(object):

    """ simple library for inverse-variance filtered maps using ffts.
    There are four steps: 1) apply apodized mask, 2) fft tqu -> teb,
    3) deconvolve transfer function teb /= transf, 4) divide by
    (signal+noise) power teb /= (cl_theory + noise_theory)
    """

    def __init__(self, obs_lib, cl_len, transf, nlev_t, nlev_p, apod):
        """
        Parameters
        --------
        obs_lib: obs object
            observation library
        cl_len: dictionary
            lensed cls.
        transf: MapSpectraTEB object
            the 2D transfer function.
        nlev_t: float
            temperature noise level in uK.arcmin
        nlev_p: float
            polarization noise level in uK.arcmin
        apod: np.ndarray or a frame-like object
            the apodization mask(s).
        """
        self.obs_lib = obs_lib
        self.cl_len = cl_len
        self.transf = transf
        self.nlev_t = nlev_t
        self.nlev_p = nlev_p
        self.apod = apod

        nl2d = transf.inverse() * transf.inverse()
        nl2d["T"] *= (nlev_t * uK * arcmin) ** 2
        nl2d["E"] *= (nlev_p * uK * arcmin) ** 2
        nl2d["B"] *= (nlev_p * uK * arcmin) ** 2

        # check
        cl2d = map_spec_utils.cl2tebfft( cl_len, nl2d )

        self.fl = (cl2d + nl2d).inverse()

        self.apod = apod

    def hashdict(self):
        ret = {
            "obs_lib": self.obs_lib.hashdict(),
            "transf": self.transf.hashdict(),
            "nlev_t": self.nlev_t,
            "nlev_p": self.nlev_p,
            "apod": hashlib.md5(self.apod.view(np.uint8)).hexdigest(),
        }

    def get_fmask(self):
        return self.apod

    def get_fl(self):
        return self.fl

    def get_flt(self):
        return self.get_fl().get_complex()["T"]

    def get_fle(self):
        return self.get_fl().get_complex()["E"]

    def get_flb(self):
        return self.get_fl().get_complex()["B"]

    def calc_tebbar(self, tqu):
        teb = calculate_teb(mult_map(tqu, self.apod))
        teb *= self.transf.inverse()
        return teb * self.get_fl()

    def get_dat_teb(self):
        return self.calc_tebbar(self.obs_lib.get_dat_tqu())

    def get_sim_teb(self, idx):
        return self.calc_tebbar(self.obs_lib.get_sim_tqu(idx))

    def get_dat_t(self):
        return self.get_dat_teb().get_complex()["T"]

    def get_dat_e(self):
        return self.get_dat_teb().get_complex()["E"]

    def get_dat_b(self):
        return self.get_dat_teb().get_complex()["B"]

    def get_sim_t(self, idx):
        return self.get_sim_teb(idx).get_complex()["T"]

    def get_sim_e(self, idx):
        return self.get_sim_teb(idx).get_complex()["E"]

    def get_sim_b(self, idx):
        return self.get_sim_teb(idx).get_complex()["B"]


class CinvFilt(object):

    """ Library for C-inv filtered maps.
    """

    def __init__(self, obs_lib, sinv_filt, ninv_filt, lib_dir, eps_min=5.0e-3):
        """
        Parameters
        ----------
        obs_lib: obs object
            observation object
        sinv_filt: ClMatrixTEB object
            sinv filter
        ninv_filt: ninv_filt object
            ninv filter
        lib_dir: string
            location to save the caches
        eps_min: float
            cutoff threshold for iterative solve
        """
        self.obs_lib = obs_lib
        self.sinv_filt = sinv_filt
        self.ninv_filt = ninv_filt

        self.lib_dir = lib_dir
        self.eps_min = eps_min

        self.hash = self.hashdict()
        if True:

            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if not os.path.exists(lib_dir + "/sim_hash.pk"):
                pk.dump(
                    self.hash,
                    open(lib_dir + "/sim_hash.pk", "wb"),
                    protocol=pk.HIGHEST_PROTOCOL,
                )
        utils.hash_check(
            pk.load(open(lib_dir + "/sim_hash.pk", "rb")), self.hash
        )

    def hashdict(self):
        return {
            "obs_lib": self.obs_lib.hash,
            "sinv_filt": self.sinv_filt.hashdict(),
            "ninv_filt": self.ninv_filt.hash,
            "eps_min": self.eps_min,
        }

    def get_fmask(self):
        return self.ninv_filt.get_fmask()

    def get_fl(self):
        """
        Analytic approximation, assuming diagonal matrices.
        X_bar(k) = C^-1 [ C^-1 + P^+ N^-1 P]^-1 P^+ N^-1 P (S_k + N_k)
                 = [N/(P^2) + (C^-1)^-1]^-1 (S_k + N_k)
                 = get_fl() * (S_k + N_k)
        operates on S_k + N_k
        """
        return (
            self.sinv_filt.inverse() + self.ninv_filt.get_fl().inverse()
        ).inverse()

    def get_flt(self):
        return self.get_fl().get_complex()["T"]

    def get_fle(self):
        return self.get_fl().get_complex()["E"]

    def get_flb(self):
        return self.get_fl().get_complex()["B"]

    def cache_teb(self, tfname, tqu_obs):
        """
        Outputs teb_filt_cinv = C-inv filtered TEB: tbar, bbar, ebar
        in e.g. EBPhi notes
        cd_solve is conjugate direction method
        x      = [fwd_op]^-1 b, where fwd_op may not be invertible
        fwd_op = C^-1 + P^+ N^-1 P    (cf. eqn 31, EBPhi notes)
        b      = P^+ N^-1 d
        """
        assert not os.path.exists(tfname)
        sinv_filt, ninv_filt = self.sinv_filt, self.ninv_filt
        pre_op = cinv_utils.opfilt_teb.PreOperatorDiag(sinv_filt, ninv_filt)
        monitor = cinv_utils.cd_monitors.MonitorBasic(
            cinv_utils.opfilt_teb.DotOperator(), iter_max=np.inf, eps_min=self.eps_min
        )
        map_ny, map_nx = tqu_obs["T"].shape
        dx = tqu_obs["T"].x_res
        dy = tqu_obs["T"].y_res
        teb_filt_cinv = MapSpectraTEB(
            map_nx=map_nx, map_ny=map_ny, dx=dx, dy=dy
        )
        cinv_utils.cd_solve.cd_solve(
            x=teb_filt_cinv,
            b=cinv_utils.opfilt_teb.calc_prep(tqu_obs, sinv_filt, ninv_filt),
            fwd_op=cinv_utils.opfilt_teb.ForwardOperator(sinv_filt, ninv_filt),
            pre_ops=[pre_op],
            dot_op=cinv_utils.opfilt_teb.DotOperator(),
            criterion=monitor,
            tr=cinv_utils.cd_solve.tr_cg,
            cache=cinv_utils.cd_solve.CacheMemory(),
        )
        # the following returns C^-1 [C^-1 + P^+ N^-1 P]^-1 P^+ N^-1 d
        teb_filt_cinv = cinv_utils.opfilt_teb.calc_fini(
            teb_filt_cinv, sinv_filt, ninv_filt
        )

        pk.dump(teb_filt_cinv, open(tfname, "wb"), protocol=pk.HIGHEST_PROTOCOL)

    def get_dat_teb(self):
        tfname = self.lib_dir + "/dat_teb_bar.pk"

        if not os.path.exists(tfname):
            self.cache_teb(tfname, self.obs_lib.get_dat_tqu())
        return pk.load(open(tfname, "rb"))

    def get_sim_teb(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_teb_bar.pk"

        if not os.path.exists(tfname):
            self.cache_teb(tfname, self.obs_lib.get_sim_tqu(idx))
        return pk.load(open(tfname, "rb"))

    def get_dat_t(self):
        return self.get_dat_teb().get_complex()["T"]

    def get_dat_e(self):
        return self.get_dat_teb().get_complex()["E"]

    def get_dat_b(self):
        return self.get_dat_teb().get_complex()["B"]

    def get_sim_t(self, idx):
        return self.get_sim_teb(idx).get_complex()["T"]

    def get_sim_e(self, idx):
        return self.get_sim_teb(idx).get_complex()["E"]

    def get_sim_b(self, idx):
        return self.get_sim_teb(idx).get_complex()["B"]


class CinvFiltMasked(object):
    """
    CinvFilt with ell mask
    """
    def __init__(
        self,
        cinv,
        lmin=None,
        lmax=None,
        lxmin=None,
        lxmax=None,
        lymin=None,
        lymax=None,
    ):
        self.cinv = cinv
        self.lmin = lmin
        self.lmax = lmax
        self.lxmin = lxmin
        self.lxmax = lxmax
        self.lymin = lymin
        self.lymax = lymax
        self.hash = self.hashdict()

    def trim(self, ret):
        return ret.get_l_masked(
            lmin=self.lmin,
            lmax=self.lmax,
            lxmin=self.lxmin,
            lxmax=self.lxmax,
            lymin=self.lymin,
            lymax=self.lymax,
        )

    def hashdict(self):
        return {
            "cinv": self.cinv.hash,
            "lmin": self.lmin,
            "lmax": self.lmax,
            "lxmin": self.lxmin,
            "lymin": self.lymin,
            "lxmax": self.lxmax,
            "lymax": self.lymax,
        }

    def get_fmask(self):
        return self.cinv.get_fmask()

    def get_fl(self):
        return self.trim(self.cinv.get_fl())

    def get_flt(self):
        return self.get_fl().get_complex()["T"]

    def get_fle(self):
        return self.get_fl().get_complex()["E"]

    def get_flb(self):
        return self.get_fl().get_complex()["B"]

    def get_dat_teb(self):
        return self.trim(self.cinv.get_dat_teb())

    def get_dat_t(self):
        return self.get_dat_teb().get_complex()["T"]

    def get_dat_e(self):
        return self.get_dat_teb().get_complex()["E"]

    def get_dat_b(self):
        return self.get_dat_teb().get_complex()["B"]

    def get_sim_teb(self, idx):
        return self.trim(self.cinv.get_sim_teb(idx))

    def get_sim_t(self, idx):
        return self.get_sim_teb(idx).get_complex()["T"]

    def get_sim_e(self, idx):
        return self.get_sim_teb(idx).get_complex()["E"]

    def get_sim_b(self, idx):
        return self.get_sim_teb(idx).get_complex()["B"]


class CinvFiltData(object):
    """
    CinvFilt that only gets data
    """
    def __init__(self, cinv):
        self.cinv = cinv
        self.hash = self.hashdict()

    def hashdict(self):
        return {"cinv": self.cinv.hash, "dat": True}

    def get_fmask(self):
        return self.cinv.get_fmask()

    def get_fl(self):
        return self.cinv.get_fl()

    def get_flt(self):
        return self.get_fl().get_complex()["T"]

    def get_fle(self):
        return self.get_fl().get_complex()["E"]

    def get_flb(self):
        return self.get_fl().get_complex()["B"]

    def get_flp(self):
        return self.cinv.get_flp()

    def get_dat_teb(self):
        return self.cinv.get_dat_teb()

    def get_dat_t(self):
        return self.cinv.get_dat_t()

    def get_dat_e(self):
        return self.cinv.get_dat_e()

    def get_dat_b(self):
        return self.cinv.get_dat_b()

    def get_sim_teb(self, idx):
        return self.get_dat_teb()

    def get_sim_t(self, idx):
        return self.get_dat_t()

    def get_sim_e(self, idx):
        return self.get_dat_e()

    def get_sim_b(self, idx):
        return self.get_dat_b()

    def get_dat_p(self):
        return self.cinv.get_dat_p()

    def get_sim_p(self, idx):
        return self.get_dat_p()


class CinvFiltSim(object):
    """a remapped library of sims that can grab an simulation
    with idx rolled by a number
    """

    def __init__(self, cinv, roll=2):
        self.cinv = cinv
        self.roll = roll

        self.get_fmask = self.cinv.get_fmask

        self.get_fl = self.cinv.get_fl
        self.get_flt = self.cinv.get_flt
        self.get_fle = self.cinv.get_fle
        self.get_flb = self.cinv.get_flb

        self.get_dat_teb = self.cinv.get_dat_teb
        self.get_dat_t = self.cinv.get_dat_t
        self.get_dat_e = self.cinv.get_dat_e
        self.get_dat_b = self.cinv.get_dat_b
        self.hash = self.hashdict()

    def hashdict(self):
        return {"cinv": self.cinv.hash, "roll": self.roll}

    def roll_i(self, i):
        return i - np.mod(i, self.roll) + np.mod(i + 1, self.roll)

    def get_sim_teb(self, i):
        return self.cinv.get_sim_teb(self.roll_i(i))

    def get_sim_t(self, i):
        return self.cinv.get_sim_t(self.roll_i(i))

    def get_sim_e(self, i):
        return self.cinv.get_sim_e(self.roll_i(i))

    def get_sim_b(self, i):
        return self.cinv.get_sim_b(self.roll_i(i))


class CinvFiltRescale(object):

    """ Library which wraps another cinv library, rescaling its output. """

    def __init__(self, ivf_lib, rescaler):
        self.ivf_lib = ivf_lib
        self.rescaler = rescaler

    def hashdict(self):
        return {
            "ivf_lib": self.ivf_lib.hashdict(),
            "rescaler": self.rescaler.hashdict(),
        }

    def get_fmask(self):
        return self.ivf_lib.get_fmask()

    def get_fl(self):
        return self.rescaler * self.ivf_lib.get_fl()

    def get_flt(self):
        return self.get_fl().get_complex()["T"]

    def get_fle(self):
        return self.get_fl().get_complex()["E"]

    def get_flb(self):
        return self.get_fl().get_complex()["B"]

    def get_dat_teb(self):
        return self.rescaler * self.ivf_lib.get_dat_teb()

    def get_dat_t(self):
        return self.get_dat_t()

    def get_dat_e(self):
        return self.get_dat_e()

    def get_dat_b(self):
        return self.get_dat_b()

    def get_sim_teb(self, idx):
        return self.rescaler * self.ivf_lib.get_sim_teb(idx)

    def get_sim_t(self, idx):
        return self.get_sim_teb(idx).get_complex()["T"]

    def get_sim_e(self, idx):
        return self.get_sim_teb(idx).get_complex()["E"]

    def get_sim_b(self, idx):
        return self.get_sim_teb(idx).get_complex()["B"]
