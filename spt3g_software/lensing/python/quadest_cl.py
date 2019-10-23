import os
import sys
import hashlib
import numpy as np
import pickle as pk

from . import utils
from spt3g import core
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum2D
from spt3g.mapspectra.map_analysis import psd_2d_to_1d


class QuadEstCl(object):
    """Class to hold Cl's from the quadratic estimates.
    """

    def __init__(self, qeA, lib_dir, qeB=None, mc_sims_mf=None):
        """
        Arguments:
        ---------
        qeA: QuadEstLib
            quadratic estimate A
        lib_dir: string
            location to cache the files
        qeB: QuadEstLib
             quadratic estimate B
        mc_sims_mf: list
            mean field from the MonteCarlo sims
        """
        if not qeB:
            qeB = qeA

        self.qeA = qeA
        self.qeB = qeB

        self.lib_dir = lib_dir

        if isinstance(mc_sims_mf, tuple):
            self.mc_sims_mfA, self.mc_sims_mfB = mc_sims_mf
        else:
            if mc_sims_mf is None:
                self.mc_sims_mfA = None
                self.mc_sims_mfB = None
            else:
                self.mc_sims_mfA = mc_sims_mf[0::2]
                self.mc_sims_mfB = mc_sims_mf[1::2]

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

            if not os.path.exists(
                lib_dir.format(prefix="outputs") + "/qecl_fcut.dat"
            ):
                mask1, mask2 = [m.flatten() for m in self.qeA.get_fmasks()]
                mask3, mask4 = [m.flatten() for m in self.qeB.get_fmasks()]

                shape = mask1.shape
                assert shape == mask2.shape
                assert shape == mask3.shape
                assert shape == mask4.shape
                npix = mask1.size

                fcut11 = np.sum(mask1 ** 2) / npix
                fcut12 = np.sum(mask1 * mask2) / npix
                fcut13 = np.sum(mask1 * mask3) / npix
                fcut14 = np.sum(mask1 * mask4) / npix
                fcut22 = np.sum(mask2 ** 2) / npix
                fcut23 = np.sum(mask2 * mask3) / npix
                fcut24 = np.sum(mask2 * mask4) / npix
                fcut33 = np.sum(mask3 ** 2) / npix
                fcut34 = np.sum(mask3 * mask4) / npix
                fcut44 = np.sum(mask4 ** 2) / npix

                fcut1234 = np.sum(mask1 * mask2 * mask3 * mask4) / npix

                np.savetxt(
                    lib_dir.format(prefix="outputs") + "/qecl_fcut.dat",
                    [
                        fcut1234,
                        fcut11,
                        fcut12,
                        fcut13,
                        fcut14,
                        fcut22,
                        fcut23,
                        fcut24,
                        fcut33,
                        fcut34,
                        fcut44,
                    ],
                )

        utils.hash_check(
            pk.load(
                open(lib_dir.format(prefix="outputs") + "/sim_hash.pk", "rb")
            ),
            self.hash,
        )

        [
            self.fcut1234,
            self.fcut11,
            self.fcut12,
            self.fcut13,
            self.fcut14,
            self.fcut22,
            self.fcut23,
            self.fcut24,
            self.fcut33,
            self.fcut34,
            self.fcut44,
        ] = np.loadtxt(lib_dir.format(prefix="outputs") + "/qecl_fcut.dat")

    def hashdict(self):
        ret = {
            "qeA": self.qeA.hash,
            "mc_sims_mfA": self.mc_sims_mfA,
            "mc_sims_mfB": self.mc_sims_mfB,
        }

        if self.qeB is not self.qeA:
            ret["qeB"] = self.qeB.hash

        return ret

    def get_qcr(self, k1, k2=None):
        if k2 == None:
            k2 = k1
        return self.qeA.get_response(k1) * self.qeB.get_response(k2)

    def get_dat_qcl(self, k1, k2=None):
        """
        Return MapSpectrum2Ds from qeA and qeB
        Note: the mean field is subtracted at this step
        """
        if k2 == None:
            k2 = k1

        qeA_qft = self.qeA.get_dat_qft(k1)
        qeB_qft = self.qeB.get_dat_qft(k2)

        if "None" not in str(type(self.mc_sims_mfA)):
            # subtract the mean-field correction
            qeA_qft -= self.qeA.get_sim_qft_mf(k1, self.mc_sims_mfA)
        if "None" not in str(type(self.mc_sims_mfB)):
            # subtract the mean-field correction
            qeB_qft -= self.qeB.get_sim_qft_mf(k2, self.mc_sims_mfB)

        return MapSpectrum2D(
            map_nx=qeA_qft.map_nx,
            dx=qeA_qft.dx,
            spec=(qeA_qft * np.conj(qeB_qft)) / self.fcut1234,
            map_ny=qeA_qft.map_ny,
            dy=qeA_qft.dy,
        )

    def get_sim_qcl(self, k1, idx, k2=None):
        """
        Return MapSpectrum2Ds from qeA and qeB for sim number 'idx'
        Note: the mean field is subtracted at this step
        """
        if k2 == None:
            k2 = k1

        qeA_qft = self.qeA.get_sim_qft(k1, idx)
        qeB_qft = self.qeB.get_sim_qft(k2, idx)
        if "None" not in str(type(self.mc_sims_mfA)):
            assert idx not in self.mc_sims_mfA
            # subtract the mean-field correction
            qeA_qft -= self.qeA.get_sim_qft_mf(k1, self.mc_sims_mfA)

        if "None" not in str(type(self.mc_sims_mfB)):
            assert idx not in self.mc_sims_mfB
            # subtract the mean-field correction
            qeB_qft -= self.qeB.get_sim_qft_mf(k2, self.mc_sims_mfB)

        return MapSpectrum2D(
            map_nx=qeA_qft.map_nx,
            dx=qeA_qft.dx,
            spec=(qeA_qft * np.conj(qeB_qft)) / self.fcut1234,
            map_ny=qeA_qft.map_ny,
            dy=qeB_qft.dy,
        )

    def get_dat_ncl_lm(self, k1, k2=None):
        return self.get_sim_ncl_lm(k1, idx=-1, k2=k2)

    def get_sim_ncl_lm(self, k1, idx, k2=None):
        if k2 == None:
            k2 = k1
        if self.qeA is self.qeB:
            k1, k2 = sorted([k1, k2])
        tfname = self.lib_dir.format(prefix="temp") + (
            "/cache_lm_%s_%s.pk"
            % (
                hashlib.sha1(
                    np.ascontiguousarray(
                        (k1, idx) + list({"k2": k2}.items())[0]
                    )
                ).hexdigest(),
                "get_sim_ncl_lm",
            )
        )
        if os.path.exists(tfname):
            try:
                return pk.load(open(tfname, "rb"))
            except EOFError:
                core.log_notice(
                    "get_sim_ncl_lm caching lm for k1=%s, k2=%s, idx=%d : %s"
                    % (k1, k2, idx, tfname), unit="Quadest_cl"
                )
        else:
            core.log_notice(
                "get_sim_ncl_lm caching lm for k1=%s, k2=%s, idx=%d : %s"
                % (k1, k2, idx, tfname), unit="Quadest_cl"
            )
        if idx == -1:
            teb_dict1, teb_dict2 = self.qeA.get_dat_tebfts()
            if self.qeA is not self.qeB:
                teb_dict3, teb_dict4 = self.qeB.get_dat_tebfts()
            else:
                teb_dict3 = teb_dict1
                teb_dict4 = teb_dict2
        else:
            teb_dict1, teb_dict2 = self.qeA.get_sim_tebfts(idx)
            if self.qeA is not self.qeB:
                teb_dict3, teb_dict4 = self.qeB.get_sim_tebfts(idx)
            else:
                teb_dict3, teb_dict4 = (teb_dict1, teb_dict2)

        estimator_ids1 = self.qeA.estimator_sets.get(k1, [k1])
        estimator_ids2 = self.qeA.estimator_sets.get(k2, [k2])
        for id1 in estimator_ids1:
            for id2 in estimator_ids2:
                f1, f2 = self.qeA.quadest.estimator_dict[id1]["fls"]
                f3, f4 = self.qeA.quadest.estimator_dict[id2]["fls"]
                f1 = teb_dict1[f1]
                f2 = teb_dict2[f2]
                f3 = teb_dict3[f3]
                f4 = teb_dict4[f4]
                cfft = MapSpectrum2D(
                    map_nx=f1.map_nx, dx=f1.dx, map_ny=f1.map_ny, dy=f1.dy
                )
                cfft_added = cfft.copy()
                cfft_added += self.qeA.quadest.covariance(
                    cfft,
                    f1 * np.conj(f3) / self.fcut13,
                    f2 * np.conj(f4) / self.fcut24,
                    estimator_id1=id1,
                    estimator_id2=id2,
                    switch_34=False,
                    conj_34=True,
                )
                cfft[:] = 0
                cfft_added += self.qeA.quadest.covariance(
                    cfft,
                    f1 * np.conj(f4) / self.fcut14,
                    f2 * np.conj(f3) / self.fcut23,
                    estimator_id1=id1,
                    estimator_id2=id2,
                    switch_34=True,
                    conj_34=True,
                )
                try:
                    spec += psd_2d_to_1d(cfft_added, lmax=5000, delta_l=1)
                except NameError:
                    spec = psd_2d_to_1d(cfft_added, lmax=5000, delta_l=1)
        pk.dump(spec, open(tfname, "wb"), protocol=pk.HIGHEST_PROTOCOL)
        return spec

    @utils.cache_pk()
    def get_qcr_lm(self, *args, **kwargs):
        spec = psd_2d_to_1d(
            self.get_qcr(*args, **kwargs), lmax=5000, delta_l=1
        )
        return spec

    @utils.cache_pk()
    def get_dat_qcl_lm(self, *args, **kwargs):
        spec = psd_2d_to_1d(
            self.get_dat_qcl(*args, **kwargs), lmax=5000, delta_l=1
        )
        return spec

    @utils.cache_pk()
    def get_sim_qcl_lm(self, *args, **kwargs):
        spec = psd_2d_to_1d(
            self.get_sim_qcl(*args, **kwargs), lmax=5000, delta_l=1
        )
        return spec
