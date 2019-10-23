"""
Load simulated lensed and un-lensed CMB skies.
This is mostly for the tutorial demonstration with lens100d data
and will not be used for 3G
"""
import os
import hashlib
import glob
import numpy as np
import pickle as pk
from spt3g import core, coordinateutils, mapmaker, mapspectra
import warnings
from . import phas
from .. import utils

def deepcopy_dict(org):
    """
    faster than deepcopy, for a dict of the simple python types.
    """
    out = dict().fromkeys(org)
    for k, v in org.items():
        try:
            out[k] = v.copy()  # dicts, sets
        except AttributeError:
            try:
                out[k] = v[:]  # lists, tuples, strings, unicode
            except TypeError:
                out[k] = v  # ints

    return out

def radec_to_thetaphi(ras, decs):
    """ras,decs are arrays degrees. return (thts, phis) in radians"""
    thts = np.pi / 2.0 - np.deg2rad(
        decs
    )  # radians. 0 at north pole up to +pi at south pole
    phis = np.deg2rad(
        ras
    )  # for continuous (unwrapped) ra, this may go negative or above 2pi.
    return (thts, phis)


class SimECP(object):
    """
    Library for loading maps in an equidistant cylindrical projection (ECP).
    """

    def __init__(
        self,
        nra=1024,
        ndec=1024,  # number of ra, dec points.
        dra=25.0,
        ddec=25.0,  # throw in ra, dec (degrees)
        center=(352.5, -55),  # (ra,dec) in degrees
        lib_dir=None,
    ):

        self.nra = nra
        self.ndec = ndec
        self.dra = dra
        self.ddec = ddec
        self.center = center

        # pixel ra centers, not edges.  unwrapped (+ and -).
        self.ras = (np.arange(0, nra) - 0.5 * (nra - 1)) * dra / nra + center[
            0
        ]
        self.decs = (
            (np.arange(0, ndec) - 0.5 * (ndec - 1)) * ddec / ndec + center[1]
        )  # pixel dec centers.  can go -90 up to +90
        (self.thts, self.phis) = radec_to_thetaphi(
            self.ras, self.decs
        )  # deg2rad in addition to pi/2 shift for theta
        self.hash = self.hashdict()

        if lib_dir is not None:
            if True:

                if not os.path.exists(lib_dir):
                    os.makedirs(lib_dir)

                if not os.path.exists(os.path.join(lib_dir, "sim_hash.pk")):
                    pk.dump(
                        self.hash,
                        open(os.path.join(lib_dir, "sim_hash.pk"), "wb"),
                        protocol=pk.HIGHEST_PROTOCOL,
                    )
            utils.hash_check(
                pk.load(open(os.path.join(lib_dir, "sim_hash.pk"), "rb")),
                self.hash,
            )

    def hashdict(self):
        return {"thts": self.thts, "phis": self.phis}

    def get_ecp_tmap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_tmap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))

    def get_ecp_pmap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_pmap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))

    def get_ecp_emap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_emap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))

    def get_ecp_bmap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_bmap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))


class SimECPSum(SimECP):
    """
    Library for loading maps in an equidistant cylindrical projection (ECP).
    """

    def __init__(self, libs):
        nra = libs[0].nra
        ndec = libs[0].ndec
        dra = libs[0].dra
        ddec = libs[0].ddec
        center = libs[0].center
        super(SimECPSum, self).__init__(
            nra=nra, ndec=ndec, dra=dra, ddec=ddec, center=center, lib_dir=None
        )

        self.libs = libs

    def hashdict(self):
        ret = {"super": super(SimECPSum, self).hash}
        for i, lib in enumerate(self.libs):
            ret[i] = lib.hashdict()
        return ret

    def get_ecp_tmap(self, idx):
        ret = self.libs[0].get_ecp_tmap(idx)
        for lib in self.libs[1:]:
            ret += lib.get_ecp_tmap(idx)
        return ret

    def get_ecp_pmap(self, idx):
        ret = self.libs[0].get_ecp_pmap(idx, elm, blm)
        for lib in self.libs[1:]:
            ret += lib.get_ecp_pmap(idx, elm, blm)
        return ret

    def get_ecp_emap(self, idx):
        ret = self.libs[0].get_ecp_emap(idx)
        for lib in self.libs[1:]:
            ret += lib.get_ecp_emap(idx)
        return ret

    def get_ecp_bmap(self, idx):
        ret = self.libs[0].get_ecp_bmap(idx)
        for lib in self.libs[1:]:
            ret += lib.get_ecp_bmap(idx)
        return ret


class SimECPUnl(SimECP):
    """
    Library for loading unlensed ECP sims.
    """

    def __init__(
        self,
        nra=1024,
        ndec=1024,  # number of ra, dec points.
        dra=25.0,
        ddec=25.0,  # throw in ra, dec (degrees)
        center=(352.5, -55),
        phas=None,
        lib_dir=None,
    ):

        self.phas = phas
        self.lib_dir = lib_dir
        if self.phas is None:
            lmax = 6000
            self.phas = phas.RandomPhase(
                4 * ((lmax + 1) * (lmax + 2)), lib_dir=self.lib_dir + "/phase"
            )

        super(SimECPUnl, self).__init__(
            nra=nra,
            ndec=ndec,
            dra=dra,
            ddec=ddec,
            center=center,
            lib_dir=lib_dir,
        )
        self.hash = self.hashdict()

    def hashdict(self):
        return {
            "parent__phase": hashlib.sha512(
                pk.dumps(self.phas.hashdict())
            ).hexdigest(),
            "super": super(SimECPUnl, self).hashdict(),
        }


class SimECPLen(SimECP):
    """
    Library for loading lensed ECP sims.
    """

    def __init__(
        self,
        nra=1024,
        ndec=1024,  # number of ra, dec points.
        dra=25.0,
        ddec=25.0,  # throw in ra, dec (degrees)
        center=(352.5, -55),
        phas=None,
        phas_phi=None,
        lib_dir=None,
    ):

        self.phas = phas
        self.phas_phi = phas_phi
        self.lib_dir = lib_dir

        if self.phas is None:
            lmax = 6000
            self.phas = phas.RandomPhase(
                4 * ((lmax + 1) * (lmax + 2)), lib_dir=self.lib_dir + "/phase"
            )

        if self.phas_phi is None:
            self.phas_phi = self.phas
        else:
            if not (
                (self.phas_phi.size == self.phas.size)
                or (self.phas_phi.size == 4 * self.phas.size)
            ):
                raise ValueError(
                    "input size for phas_phi must be either (a) equal to or (b) 1/4 the size of phas."
                )

        super(SimECPLen, self).__init__(
            nra=nra,
            ndec=ndec,
            dra=dra,
            ddec=ddec,
            center=center,
            lib_dir=lib_dir,
        )
        self.hash = self.hashdict()

    def hashdict(self):

        return {
            "parent__phase": hashlib.sha512(
                pk.dumps(self.phas.hashdict())
            ).hexdigest(),
            "parent__phase_phi": hashlib.sha512(
                pk.dumps(self.phas_phi.hashdict())
            ).hexdigest(),
            "super": super(SimECPLen, self).hashdict(),
        }

    def get_ecp_kmap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_kmap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))

    def get_ecp_dmap(self, idx):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_ecp_dmap.npy"
        try:
            return np.load(tfname)
        except IOError:
            return pk.load(open(tfname[:-3] + "pk", "rb"))


class SimProj(object):
    """
    Library for loading maps interpolated to an SPTpol projection from ECP.
    """

    def __init__(
        self,
        nx,
        ny,
        reso_arcmin,
        lib_ecp,
        center=(352.5, -55),
        proj=5,
        lib_dir=None,
        cache=False,
    ):
        self.nx = nx  # ra-like
        self.ny = ny  # dec-like
        self.reso_arcmin = reso_arcmin
        self.lib_ecp = lib_ecp
        self.center = center
        self.proj = proj
        self.lib_dir = lib_dir
        self.cache = cache
        self.hash = self.hashdict()

        if True:
            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if not os.path.exists(os.path.join(lib_dir, "sim_hash.pk")):
                pk.dump(
                    self.hash,
                    open(os.path.join(lib_dir, "sim_hash.pk"), "wb"),
                    protocol=pk.HIGHEST_PROTOCOL,
                )

        utils.hash_check(
            pk.load(open(os.path.join(lib_dir, "sim_hash.pk"), "rb")),
            self.hash,
        )

    def hashdict(self):
        return {
            "nx": self.nx,
            "ny": self.ny,
            "reso_arcmin": self.reso_arcmin,
            "lib_ecp": self.lib_ecp.hash,
            "center": self.center,
            "proj": self.proj,
        }

    def get_sim_kap(self, idx, fl=1.0):
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_kap.g3"
        if (self.cache is False) or (not os.path.exists(tfname)):
            core.log_warn("kappa map does not exist")
        else:
            kmap = core.G3File(tfname).next()

        return kmap["none"]  # (kmap.get_rfft() * fl).get_rmap()

    # allow call SimProj.get_sim_kappa()
    get_sim_kappa = get_sim_kap


class SimScan(object):
    """
    Library for loading SPTpol mock-observed ECP maps. 
    """

    def __init__(
        self,
        nx,
        ny,
        lib_ecp,
        idfs=None,
        lib_dir=None,
        master_configfile=None,
        analysis_function="",
        analysis_function_kwargs={},
        preprocessing_function=[],
        preprocessing_function_kwargs=[{}],
        xtalk_files=None,
        use_leftgoing=None,
    ):

        self.nx = nx  # ra-like
        self.ny = ny  # dec-like
        self.lib_ecp = lib_ecp
        self.idfs = idfs
        self.lib_dir = lib_dir
        self.use_leftgoing = use_leftgoing

        # xtalk sims
        self.xtalk_files = xtalk_files
        if xtalk_files is None:
            self.do_xtalk = False
        else:
            self.do_xtalk = True
            self.setup_xtalk_TP_deprojection()

        self.master_configfile = master_configfile
        self.analysis_function = analysis_function
        self.analysis_function_kwargs = analysis_function_kwargs
        self.preprocessing_function = preprocessing_function
        self.preprocessing_function_kwargs = preprocessing_function_kwargs
        self.hash = self.hashdict()
        try:
            self.center = self.analysis_function_kwargs["map_center"]
        except KeyError:
            self.center = "source"

        if True:

            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if not os.path.exists(os.path.join(lib_dir, "sim_hash.pk")):
                pk.dump(
                    self.hash,
                    open(os.path.join(lib_dir, "sim_hash.pk"), "wb"),
                    protocol=pk.HIGHEST_PROTOCOL,
                )

        utils.hash_check(
            pk.load(open(os.path.join(lib_dir, "sim_hash.pk"), "rb")),
            self.hash,
        )

    def hashdict(self):
        # Create md5sum hash of master_config file
        if self.master_configfile is not None:
            try:
                master_configfile_hash = hashlib.md5(
                    open(self.master_configfile).read().encode("utf-8")
                ).hexdigest()
            except:
                core.log_warn(
                    "library_scan.hashdict(): Could not create md5 hash of master_config.  Quitting."
                )
                assert 0
        else:
            master_configfile_hash = ""

        kwargs_for_hash = deepcopy_dict(self.analysis_function_kwargs)
        ptsrc_md5hash = ""

        # build the hashdict
        ret = {
            "nx": self.nx,
            "ny": self.ny,
            "lib_ecp": self.lib_ecp.hash,
            "idfs": self.idfs,
            "use_leftgoing": self.use_leftgoing,
            "master_configfile_md5hash": master_configfile_hash,
            "ptsrc_md5hash": ptsrc_md5hash,
            "analysis_function": self.analysis_function,
            "analysis_function_kwargs": kwargs_for_hash,
            "preprocessing_function": self.preprocessing_function,
            "preprocessing_function_kwargs": self.preprocessing_function_kwargs,
        }

        # xtalk files
        if self.do_xtalk:
            ret["xtalk_files"] = self.xtalk_files
        return ret

    def get_sim_map(self, idx):
        """
        Return the simulated map number 'idx'.
        Units: uK
        Note: a) this map is weighted
              b) you need to call flattenPol() before taking a power spectrum
        """
        tfname = self.lib_dir + "/sim_" + ("%04d" % idx) + "_buf.npy"
        if not os.path.exists(tfname):
            if not os.path.exists(tfname[:-3] + "pk"):
                core.log_warn("The simulation files do not exist.")

        buf = np.load(tfname)
        weight = np.load(self.lib_dir + "/weight.npy")

        analysis_function_kwargs = deepcopy_dict(self.analysis_function_kwargs)
        if "scans" in analysis_function_kwargs:
            del analysis_function_kwargs["scans"]
        if "scan_poly_order" in analysis_function_kwargs:
            del analysis_function_kwargs["scan_poly_order"]

        smap = core.G3Frame(core.G3FrameType.Map)
        smap["Id"] = "sptpol_sim"
        smT = buf[(0 * self.nx * self.ny) : (1 * self.nx * self.ny)].reshape(
            (self.ny, self.nx)
        )
        smQ = buf[(1 * self.nx * self.ny) : (2 * self.nx * self.ny)].reshape(
            (self.ny, self.nx)
        )
        smU = buf[(2 * self.nx * self.ny) : (3 * self.nx * self.ny)].reshape(
            (self.ny, self.nx)
        )

        t_map = coordinateutils.FlatSkyMap(
            smT,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            is_weighted=True,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
            units=core.G3TimestreamUnits.Tcmb,
            pol_type=coordinateutils.MapPolType.T,
        )
        q_map = coordinateutils.FlatSkyMap(
            smQ,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            is_weighted=True,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
            units=core.G3TimestreamUnits.Tcmb,
            pol_type=coordinateutils.MapPolType.T,
        )
        u_map = coordinateutils.FlatSkyMap(
            smU,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            is_weighted=True,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
            units=core.G3TimestreamUnits.Tcmb,
            pol_type=coordinateutils.MapPolType.T,
        )
        WeightsMap = coordinateutils.G3SkyMapWeights(t_map)
        tt_w = weight[:, :, 0, 0].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        tt__w = coordinateutils.FlatSkyMap(
            tt_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.TT = tt__w
        tq_w = weight[:, :, 0, 1].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        tq__w = coordinateutils.FlatSkyMap(
            tq_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.TQ = tq__w
        tu_w = weight[:, :, 0, 2].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        tu__w = coordinateutils.FlatSkyMap(
            tu_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.TU = tu__w
        qq_w = weight[:, :, 1, 1].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        qq__w = coordinateutils.FlatSkyMap(
            qq_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.QQ = qq__w
        qu_w = weight[:, :, 1, 2].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        qu__w = coordinateutils.FlatSkyMap(
            qu_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.QU = qu__w
        uu_w = weight[:, :, 2, 2].copy(
            order="C"
        )  # makes contiguous so b.p. works with python slices
        uu__w = coordinateutils.FlatSkyMap(
            uu_w,
            self.analysis_function_kwargs["reso_arcmin"] * core.G3Units.arcmin,
            proj=coordinateutils.MapProjection(
                self.analysis_function_kwargs["proj"]
            ),
            alpha_center=self.analysis_function_kwargs["map_center"][0]
            * core.G3Units.deg,
            delta_center=self.analysis_function_kwargs["map_center"][1]
            * core.G3Units.deg,
            coord_ref=coordinateutils.MapCoordReference.Equatorial,
        )
        WeightsMap.UU = uu__w
        smap["T"] = t_map
        smap["Q"] = q_map
        smap["U"] = u_map
        smap["Wpol"] = WeightsMap

        return smap

    def get_sim_tqu(self, idx):
        smap = self.get_sim_map(idx)
        mapmaker.mapmakerutils.RemoveWeightModule(smap)
        smap = mapspectra.map_analysis.flatten_pol(smap)

        return smap
