import os, sys, hashlib, pdb
import numpy as np
import pickle as pk
from . import utils
from . import quadest_utils
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum1D, MapSpectrum2D
from spt3g.mapspectra.basicmaputils import map_to_ft
from .map_spec_utils import mult_map


class QuadEstLib(object):
    """ Class holding a simulated quadratic estimate of a field

    Attributes:
    ---------
    lib_dir: string
        directory of this library's data products
    cl_unl: dictionary
        unlensed spectrum
    cl_len: dictionary
        lensed spectrum
    ivfs1: cinv library
        cinv-filtered data library for field 1
    ivfs2: cinv library
        cinv-filtered data library for field 2
    """

    def __init__(self, cl_unl, cl_len, ivfs1=None, lib_dir=None, ivfs2=None):

        self.estimator_sets = {
            "Phi_set": (
                "Phi_TT",
                "Phi_EE",
                "Phi_TE",
                "Phi_ET",
                "Phi_TB",
                "Phi_BT",
                "Phi_EB",
                "Phi_BE",
            ),
            "Phi_curl_set": (
                "Phi_TT_curl",
                "Phi_EE_curl",
                "Phi_TE_curl",
                "Phi_ET_curl",
                "Phi_TB_curl",
                "Phi_BT_curl",
                "Phi_EB_curl",
                "Phi_BE_curl",
            ),
            "Phi_TE_set": ("Phi_TE", "Phi_ET"),
            "Phi_TE_curl_set": ("Phi_TE_curl", "Phi_ET_curl"),
            "Phi_TB_set": ("Phi_BT", "Phi_TB"),
            "Phi_TB_curl_set": ("Phi_TB_curl", "Phi_BT_curl"),
            "Phi_EB_set": ("Phi_EB", "Phi_BE"),
            "Phi_EB_curl_set": ("Phi_EB_curl", "Phi_BE_curl"),
            "Phi_pol_set": ("Phi_EB", "Phi_BE", "Phi_EE"),
            "Phi_pol_curl_set": ("Phi_EB_curl", "Phi_BE_curl", "Phi_EE_curl"),
        }
        # TODO: add back in "tfe" value, a multiplicative weighting factor i
        # associated with each estimator. previously they were all set to 1.
        # should it go here, or in quadest_utils?

        if ivfs2 == None:
            ivfs2 = ivfs1

        self.lib_dir = lib_dir

        self.cl_unl = cl_unl
        self.cl_len = cl_len
        self.ivfs1 = ivfs1
        self.ivfs2 = ivfs2
        self.quadest = quadest_utils.QuadEst(cl_len)

        self.hash = self.hashdict()

        if True:
            if not os.path.exists(lib_dir.format(prefix="temp")):
                os.makedirs(lib_dir.format(prefix="temp"))

            if not os.path.exists(lib_dir.format(prefix="outputs")):
                os.makedirs(lib_dir.format(prefix="outputs"))

            if not os.path.exists(
                lib_dir.format(prefix="outputs") + "/sim_hash.pk"
            ):
                pk.dump(
                    self.hash,
                    open(
                        lib_dir.format(prefix="outputs") + "/sim_hash.pk", "wb"
                    ),
                    protocol=pk.HIGHEST_PROTOCOL,
                )

        utils.hash_check(
            pk.load(
                open(lib_dir.format(prefix="outputs") + "/sim_hash.pk", "rb")
            ),
            self.hash,
        )

    def get_fmasks(self):
        return self.ivfs1.get_fmask(), self.ivfs2.get_fmask()

    def hashdict(self):
        ret = {
            # "cl_unl": self.cl_unl.hashdict(),
            # "cl_len": self.cl_len.hashdict(),
            "ivfs1": self.ivfs1.hash
        }

        if self.ivfs2 is not self.ivfs1:
            ret["ivfs2"] = self.ivfs2.hash

        return ret

    def get_response(self, id, npad=2):
        """
        Wrapper for quadest_utils.add_response()

        Parameters:
        ---------
        id: string
            The estimator to use
        npad:  int
            the padding factor in the FFT convolution calculation
        """
        teb_dict1 = self.ivfs1.get_fl().get_complex()
        teb_dict2 = (
            teb_dict1
            if self.ivfs1 is self.ivfs2
            else self.ivfs2.get_fl().get_complex()
        )

        ret = MapSpectrum2D(
            map_nx=teb_dict1["T"].map_nx,
            dx=teb_dict1["T"].dx,
            map_ny=teb_dict1["T"].map_ny,
            dy=teb_dict1["T"].dy,
        )

        estimator_ids = self.estimator_sets.get(id, [id])
        for estimator_id in estimator_ids:
            cfft = MapSpectrum2D(
                map_nx=teb_dict1["T"].map_nx,
                dx=teb_dict1["T"].dx,
                map_ny=teb_dict1["T"].map_ny,
                dy=teb_dict1["T"].dy,
            )
            teb_key1, teb_key2 = self.quadest.estimator_dict[estimator_id][
                "fls"
            ]
            fl1, fl2 = teb_dict1[teb_key1], teb_dict2[teb_key2]
            ret += self.quadest.fill_response(
                cfft, fl1, fl2, estimator_id, estimator_id
            )  # * tfs
        return ret
        # TODO: add functionality for two different estimators
        # (i.e add second id argument).

    def get_qft(self, id, teb_dict1, teb_dict2):
        cfft = MapSpectrum2D(
            map_nx=teb_dict1["T"].map_nx,
            dx=teb_dict1["T"].dx,
            map_ny=teb_dict1["T"].map_ny,
            dy=teb_dict1["T"].dy,
        )
        estimator_ids = self.estimator_sets.get(id, [id])
        for estimator_id in estimator_ids:
            teb_key1, teb_key2 = self.quadest.estimator_dict[estimator_id][
                "fls"
            ]
            f1, f2 = teb_dict1[teb_key1], teb_dict2[teb_key2]
            cfft += self.quadest.eval(f1, f2, estimator_id)  # * tfe
        return cfft

    def get_dat_tebfts(self):
        """
        Return the data T,E, and B FFT fields
        """
        teb_dict1 = self.ivfs1.get_dat_teb().get_complex()

        teb_dict2 = (
            teb_dict1
            if self.ivfs1 is self.ivfs2
            else self.ivfs2.get_dat_teb().get_complex()
        )

        return teb_dict1, teb_dict2

    def get_sim_tebfts(self, i):
        """
        Return the simulated T,E, and B FFT fields
        """
        teb_dict1 = self.ivfs1.get_sim_teb(i).get_complex()
        teb_dict2 = (
            teb_dict1
            if self.ivfs1 is self.ivfs2
            else self.ivfs2.get_sim_teb(i).get_complex()
        )
        return teb_dict1, teb_dict2

    def get_dat_qft(self, id):
        """
        Return the Fourier-space quadratic estimate of type k, from the data
        example: k='ptt'
        """
        teb_dict1, teb_dict2 = self.get_dat_tebfts()
        return self.get_qft(id, teb_dict1, teb_dict2)

    def get_sim_qft(self, id, i):
        """
        Return the Fourier-space quadratic estimate of type k, from the sims
        """
        if i == -1:
            return self.get_dat_qft(id)
        teb_dict1, teb_dict2 = self.get_sim_tebfts(i)
        return self.get_qft(id, teb_dict1, teb_dict2)

    def get_sim_qft_mf(self, k, idxs):
        """
        Return the mean-field Fourier-space quadratic estimate of type k,
        from the sims
        """
        tfname = self.lib_dir.format(
            prefix="temp"
        ) + "/sim_qft_mf_%s_%s.pk" % (
            k,
            hashlib.sha1(np.ascontiguousarray(idxs)).hexdigest(),
        )
        if not os.path.exists(tfname):
            qft_mf_avg = utils.AverageObjects()
            for i, idx in utils.enumerate_progress(idxs, "get_sim_qft_mf"):
                qft_mf_avg.add(self.get_sim_qft(k, idx))
            pk.dump(
                qft_mf_avg.get(),
                open(tfname, "wb"),
                protocol=pk.HIGHEST_PROTOCOL,
            )
        return pk.load(open(tfname, "rb"))


class QuadEstLibKappa(object):
    """
    Class which pretends to be a quadratic lensing estimator,
    but actually just returns lensing kappa for some cmb library.
    """

    def __init__(self, qe, cmbs):
        self.qe = qe
        self.cmbs = cmbs
        self.hash = self.hashdict()

    def get_fmasks(self):
        m1, m2 = self.qe.ivfs1.get_fmask(), self.qe.ivfs2.get_fmask()
        return np.ones(m1.shape), np.ones(m2.shape)

    def hashdict(self):
        return {"qe": self.qe.hash, "cmbs": self.cmbs.hash}

    def get_response(self, id, npad=2):
        assert id[0] == "P"  # this is a phi estimator.
        return self.qe.get_response(id, npad=npad)

    def get_sim_qft(self, k, idx):
        assert idx >= 0
        spec_scale_fac = MapSpectrum1D(
            lbins=np.arange(0.0, 100001.0) - 0.5,
            spec=np.nan_to_num(2.0 / np.arange(0.0, 100000.0) ** 2),
        )
        kap_ft = map_to_ft(self.cmbs.get_sim_kap(idx))
        return (
            kap_ft * spec_scale_fac.get_2d(kap_ft)* self.get_response(k)
        )  # convert kappa->phi, then convolve to look like qe.

    def get_sim_qft_mf(self, k, idxs):
        assert 0


class QuadEstLibKappaMasked(object):
    """
    Class which pretends to be a quadratic lensing estimator,
    but actually just returns lensing kappa for some cmb library.
    """

    def __init__(self, qe, cmbs, mask):
        self.qe = qe
        self.cmbs = cmbs
        self.mask = mask

    def get_fmasks(self):
        m1, m2 = self.qe.ivfs1.get_fmask(), self.qe.ivfs2.get_fmask()
        return np.ones(m1.shape), np.ones(m2.shape)

    def hashdict(self):
        return {
            "qe": self.qe.hashdict(),
            "cmbs": self.cmbs.hashdict(),
            "mask": hashlib.sha1(
                self.mask.flatten().view(np.uint8)
            ).hexdigest(),
        }

    def get_response(self, id, npad=2):
        assert id[0] == "p"  # this is a phi estimator.
        return self.qe.get_response(id, npad=npad)

    def get_sim_qft(self, k, idx):
        assert idx >= 0
        spec_scale_fac = MapSpectrum1D(
            lbins=np.arange(0.0, 100001.0) - 0.5,
            spec=np.nan_to_num(2.0 / np.arange(0.0, 100000.0) ** 2),
        )
        kap_ft = map_to_ft(self.cmbs.get_sim_kap(idx) * self.mask)
        return (
            kap_ft * spec_scale_fac.get_2d(kap_ft) * self.get_response(k)
        )  # convert kappa->phi, then convolve to look like qe.

    def get_sim_qft_mf(self, k, idxs):
        assert 0


class QuadEstLibData(object):
    """
    A remapped QuadEstLib that only gets data
    """
    def __init__(self, qest):
        self.qest = qest

    def hashdict(self):
        return {"qest": self.qest.hashdict(), "dat": True}

    def get_fmasks(self):
        return self.qest.get_fmasks()

    def get_response(self, *args, **kwargs):
        return self.qest.get_response(*args, **kwargs)

    def get_qft(self, *args, **kwargs):
        assert 0

    def get_dat_tebffts(self):
        assert 0

    def get_sim_tebffts(self, idx):
        assert 0

    def get_dat_qft(self, k):
        return self.qest.get_dat_qft(k)

    def get_sim_qft(self, k, idx):
        return self.qest.get_dat_qft(k)

    def get_sim_qft_mf(self, k, idxs):
        return self.qest.get_sim_qft_mf(k, idxs)
