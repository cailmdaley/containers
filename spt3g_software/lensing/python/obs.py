"""
Purpose: Collect information about observations including noise,
cmb, foregrounds, etc.

Classes
   ObsHalfNoise()
      Purpose: Library for data and simulated observations,
               using half-observation differencing to estimate noise.
      get_ninv
      get_dat_tqu
      get_sim_tqu
      get_sim_tqu_sky

   ObsHomogeneousNoise()
      Purpose: Library for data and simulated observations,
               with a homogeneous noise model.

   ObsHomogeneousNoiseSim()
      Purpose: Library for testing on the simulation,
               with a homogeneous noise model.
"""

import os
import sys
import hashlib
import numpy as np
import pickle as pk
from builtins import str
from builtins import object
from spt3g import core, mapmaker, mapspectra
from spt3g.mapspectra import map_analysis
from . import map_spec_utils, utils

uK = core.G3Units.uK
arcmin = core.G3Units.arcmin

def threshold(maps, vmin, vmax=None, vcut=0.0):
    """ returns a new, thresholded version of the current map.
    threshold(v) -> set all pixels which don't satisfy
    (-|v| < val < |v|) equal to vcut.
    threshold(min,max) -> set all pixels which don't satisfy
    (vmin < val < vmax) equal to vcut.
    """
    if vmax is None:
        vmin = -np.abs(vmin)
        vmax = +np.abs(vmin)
    assert vmin < vmax
    for map_type in ["T", "Q", "U"]:
        m = np.array(maps[map_type])
        m[np.where(m < vmin)] = vcut
        m[np.where(m > vmax)] = vcut
        np.asarray(maps[map_type])[:] = m

    return maps


class ObsHalfNoise(object):
    """
    Library for data and simulated observations,
    using half-observation differencing to estimate noise.
    """

    def __init__(
        self,
        transf,
        sky_lib,
        half_dir,
        tq=0.0,
        tu=0.0,
        qu=0.0,
        fp=True,
        thresh=1000.0 * uK,
        tcal=1.0,
        pcal=1.0,
    ):
        """
        Parameters:
        ---------
        tq: float
            T-Q leakage term
        tu: float
            T-U leakage term
        fp: bool
            if fp=True, call flattenPol()
        thresh: float
            cut pixels in sim-maps with an absolute
            value greater than "thresh"
        tcal: float
            the absolute calibration of T,Q,and U in map units.
            Note: In some analyses, the tcal factor is included in the beam in
            "transf".
        pcal: float
            the calibration of Q and U maps.
            this is applied to Q&U in addition to tcal.
    """
        self.transf = transf
        self.sky_lib = sky_lib
        self.half_dir = half_dir

        self.tq = tq
        self.tu = tu
        self.qu = qu
        self.fp = fp
        self.tcal = tcal
        self.pcal = pcal

        self.thresh = thresh
        self.hash = self.hashdict()

    def hashdict(self):
        ret = {
            "sky_lib": self.sky_lib.hash,
            "half_dir": self.half_dir,
            "tq": self.tq,
            "tu": self.tu,
            "qu": self.qu,
            "thresh": self.thresh,
        }

        if not self.fp:
            ret["fp"] = self.fp

        if type(self.transf) == np.ndarray:
            ret["transf"] = hashlib.sha1(
                self.transf.view(np.uint8)
            ).hexdigest()
        else:
            ret["transf"] = self.transf.hashdict()

        if self.tcal != 1.0:
            ret["tcal"] = self.tcal
        if self.pcal != 1.0:
            ret["pcal"] = self.pcal

        return ret

    def fixmap(self, hmap):
        """
        Remove weight, flatten pol
        Remove T->P leakage,
        rotate qu
        """
        if hmap["T"].is_weighted:
            mapmaker.mapmakerutils.RemoveWeightModule(hmap)
        if self.fp == True:
            hmap = mapspectra.map_analysis.flatten_pol(hmap)
        else:
            core.log_notice("fixmap-- not flattening pol!", unit="Obs")

        if self.tq != 0.0:
            qmap = hmap["Q"] + self.tq * -hmap["T"]
            del hmap["Q"]
            hmap["Q"] = qmap
        if self.tu != 0.0:
            umap = hmap["U"] + self.tu * -hmap["T"]
            del hmap["U"]
            hmap["U"] = umap
        if self.qu != 0.0:
            core.log_info("rotating QU", unit="Obs")
            qiu = np.array(hmap["Q"]) + 1.0j * np.array(hmap["U"])
            qiu *= np.exp(-2.0j * np.pi * self.qu / 180.0)
            qmap = coordinateutils.FlatSkyMap(
                qiu.real,
                hmap["T"].res,
                proj=hmap["T"].proj,
                alpha_center=hmap["T"].alpha_center,
                delta_center=hmap["T"].delta_center,
                coord_ref=hmap["T"].coord_ref,
            )
            umap = coordinateutils.FlatSkyMap(
                qiu.imag,
                hmap["T"].res,
                proj=hmap["T"].proj,
                alpha_center=hmap["T"].alpha_center,
                delta_center=hmap["T"].delta_center,
                coord_ref=hmap["T"].coord_ref,
            )
            del hmap["Q"]
            hmap["Q"] = qmap
            del hmap["U"]
            hmap["U"] = umap

        return hmap

    def get_ninv(self):
        """
        Return N^-1 as the coadded weight map
        """
        halfA = self.half_dir + "/half_A_%d.g3" % (0)
        halfB = self.half_dir + "/half_B_%d.g3" % (0)
        hmap = map_analysis.add_two_maps(halfA, halfB)
        ninv = hmap["Wpol"]
        return ninv

    def get_dat_tqu(self):
        """
        Return the coadded data (observed sky)
        """
        halfA = self.half_dir + "/half_A_%d.g3" % (0)
        halfB = self.half_dir + "/half_B_%d.g3" % (0)
        hmap = self.fixmap(map_analysis.add_two_maps(halfA, halfB))

        # Apply calibration factors
        # hmap *= self.tcal   # apply tcal to all maps [T,Q,U]

        # hmap['Q'].map *= self.pcal   # additionally, apply pcal to Q and U
        # hmap['U'].map *= self.pcal
        tmap = hmap["T"] * self.tcal
        qmap = hmap["Q"] * self.tcal * self.pcal
        umap = hmap["U"] * self.tcal * self.pcal
        del hmap["T"]
        hmap["T"] = tmap
        del hmap["Q"]
        hmap["Q"] = qmap
        del hmap["U"]
        hmap["U"] = umap

        # Apply a threshold to the maps
        hmap = threshold(hmap, self.thresh)

        return hmap

    def get_sim_tqu(self, isim):
        """
        Return the simulated sky + noise
        Note: noise (but not sim cmb skies) should
        be scaled by the calibration factors.
        """
        halfA = self.half_dir + "/half_A_%d.g3" % isim
        halfB = self.half_dir + "/half_B_%d.g3" % isim

        hmap = self.fixmap(map_analysis.subtract_two_maps(halfA, halfB))

        # Apply calibration factors to the noise
        tmap = hmap["T"] * self.tcal
        qmap = hmap["Q"] * self.tcal * self.pcal
        umap = hmap["U"] * self.tcal * self.pcal
        del hmap["T"]
        # additional factor of two for noise estimation
        hmap["T"] = tmap / 2.0
        del hmap["Q"]
        hmap["Q"] = qmap / 2.0
        del hmap["U"]
        hmap["U"] = umap / 2.0

        cmb_tqu = self.sky_lib.get_sim_tqu(isim)

        # Apply a threshold to the maps
        hmap = threshold(hmap, self.thresh)
        cmb_tqu = threshold(cmb_tqu, self.thresh)
        cmb_tqu = (
            map_spec_utils.calculate_teb(cmb_tqu) * self.transf
        ).get_tqu()
        for k in ["T", "Q", "U"]:
            np.asarray(cmb_tqu[k])[:] = np.array(cmb_tqu[k]) + np.array(
                hmap[k]
            )
        return cmb_tqu

    def get_sim_tqu_sky(self, isim):
        """
        Return the simulated sky only, no noise
        """
        cmb_tqu = self.sky_lib.get_sim_tqu(isim)

        # Apply a threshold to the maps
        cmb_tqu = threshold(cmb_tqu, self.thresh)

        return (map_spec_utils.calculate_teb(cmb_tqu) * self.transf).get_tqu()


class ObsHomogeneousNoise(object):
    """
    Library for data and simulated observations,
    with a homogeneous noise model.
    """

    def __init__(
        self, transf, sky_lib, nlev_t=0.0, nlev_p=0.0, thresh=1000.0 * uK
    ):
        """
        Parameters:
        --------
        thresh: float
            cut pixels in sim-maps with an absolute value
             greater than "thresh"
        nlev_t: float
            temperature noise level in uK.arcmin
        nlev_p: float
            polarization (q and u) noise level in uK.arcmin
        """
        self.transf = transf
        self.sky_lib = sky_lib

        self.nlev_t = nlev_t
        self.nlev_p = nlev_p

        self.thresh = thresh
        self.hash = self.hashdict()

    def hashdict(self):
        ret = {
            "sky_lib": self.sky_lib.hash,
            "nlev_t": self.nlev_t,
            "nlev_p": self.nlev_p,
            "thresh": self.thresh,
        }

        if type(self.transf) == np.ndarray:
            ret["transf"] = hashlib.sha1(
                self.transf.view(np.uint8)
            ).hexdigest()
        else:
            ret["transf"] = self.transf.hashdict()

        return ret

    def get_dat_tqu(self):
        ret = self.get_sim_tqu(0)
        for k in ["T", "Q", "U"]:
            np.asarray(ret[k])[:] = 0
        return ret

    def get_sim_tqu(self, isim):
        """
        Return the simulated sky + noise
        """
        cmb_tqu = self.sky_lib.get_sim_tqu(isim)

        # Apply a threshold to the maps
        cmb_tqu = threshold(cmb_tqu, self.thresh)

        # Add the noise realizations
        hmap ={}
        dx = cmb_tqu["T"].x_res
        dy = cmb_tqu["T"].y_res

        hmap["T"] = (
            np.random.standard_normal(cmb_tqu["T"].shape)
            * self.nlev_t * uK * arcmin
            / np.sqrt(dx * dy)
        )
        hmap["Q"] = (
            np.random.standard_normal(cmb_tqu["Q"].shape)
            * self.nlev_p * uK * arcmin
            / np.sqrt(dx * dy)
        )
        hmap["U"] = (
            np.random.standard_normal(cmb_tqu["U"].shape)
            * self.nlev_p * uK * arcmin
            / np.sqrt(dx * dy)
        )

        cmb_tqu = (
            map_spec_utils.calculate_teb(cmb_tqu) * self.transf
        ).get_tqu()

        for k in ["T", "Q", "U"]:
            np.asarray(cmb_tqu[k])[:] = np.array(cmb_tqu[k]) + np.array(
                hmap[k]
            )
        return cmb_tqu


class ObsHomogeneousNoiseSim(object):
    """
    Library for data and simulated observations, with a homogeneous noise model.
    Treat sim with index 0 as the data
    """

    def __init__(
        self, transf, sky_lib, nlev_t=0.0, nlev_p=0.0, thresh=1000.0 * uK
    ):
        """
        Selected Parameters:
        ---------
        thresh: float
            cut pixels in sim-maps with an absolute value
            greater than "thresh"
        nlev_t: float
            temperature noise level in uK.arcmin
        nlev_p: float
            polarization (q and u) noise level in uK.arcmin
        """
        self.transf = transf
        self.sky_lib = sky_lib

        self.nlev_t = nlev_t
        self.nlev_p = nlev_p

        self.thresh = thresh
        self.hash = self.hashdict()

    def hashdict(self):
        ret = {
            "sky_lib": self.sky_lib.hash,
            "nlev_t": self.nlev_t,
            "nlev_p": self.nlev_p,
            "thresh": self.thresh,
        }

        if type(self.transf) == np.ndarray:
            ret["transf"] = hashlib.sha1(
                self.transf.view(np.uint8)
            ).hexdigest()
        else:
            ret["transf"] = self.transf.hashdict()

        return ret

    def get_dat_tqu(self):
        return self.get_sim_tqu(0)

    def get_sim_tqu(self, isim):
        """
        Return the simulated sky + noise
        """
        cmb_tqu = self.sky_lib.get_sim_tqu(isim)

        # Apply a threshold to the maps
        cmb_tqu = cmb_tqu.threshold(self.thresh)

        dx = cmb_tqu["T"].x_res
        dy = cmb_tqu["T"].y_res

        # Add the noise realizations
        hmap ={}

        hmap["T"] = (
            np.random.standard_normal(cmb_tqu["T"].shape)
            * self.nlev_t * uK * arcmin
            / np.sqrt(dx * dy)
        )
        hmap["Q"] = (
            np.random.standard_normal(cmb_tqu["Q"].shape)
            * self.nlev_p * uK * arcmin
            / np.sqrt(dx * dy)
        )
        hmap["U"] = (
            np.random.standard_normal(cmb_tqu["U"].shape)
            * self.nlev_p * uK * arcmin
            / np.sqrt(dx * dy)
        )
        cmb_tqu = (
            map_spec_utils.calculate_teb(cmb_tqu) * self.transf
        ).get_tqu()

        for k in ["T", "Q", "U"]:
            np.asarray(cmb_tqu[k])[:] = np.array(cmb_tqu[k]) + np.array(
                hmap[k]
            )
        return cmb_tqu


class SumTQU(object):
    """
    Library for summing multiple components.
    """

    def __init__(self, cpts):
        """
        Parameters:
        --------
        cpts: list of libraries.
            The first library should be the signal (data or sims)
        """
        assert len(cpts) > 0
        self.cpts = cpts
        self.hash = self.hashdict()

    def hashdict(self):
        return {i: cpt.hashdict() for i, cpt in enumerate(self.cpts)}

    def get_dat_tqu(self):
        return self.cpts[0].get_dat_tqu()

    def get_sim_tqu(self, idx):
        ret = self.cpts[0].get_sim_tqu(idx)
        for k in ["T", "Q", "U"]:
            temp = ret[k]
            for cpt in self.cpts[1:]:
                np.asarray(temp)[:] = temp + np.array(cpt.get_sim_tqu(idx)[k])
            del ret[k]
            ret[k] = temp
        return ret
