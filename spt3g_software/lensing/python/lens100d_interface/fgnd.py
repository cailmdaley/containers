import os
import sys
import hashlib
import numpy as np
import pickle as pk
import healpy as hp
from spt3g import core, coordinateutils
from .. import utils, map_spec_utils
from ..map_spec_utils import MapSpectraTEB
from . import phas

uK = core.G3Units.uK


def tebfft(nx, ny, dx, dy, tcl):
    """A function that generates Gaussian realization of a power spectrum
    """
    tfft = (
        np.random.standard_normal((ny, nx // 2 + 1))
        + 1.0j * np.random.standard_normal((ny, nx // 2 + 1))
    ) / np.sqrt(2.0)
    efft = (
        np.random.standard_normal((ny, nx // 2 + 1))
        + 1.0j * np.random.standard_normal((ny, nx // 2 + 1))
    ) / np.sqrt(2.0)
    bfft = (
        np.random.standard_normal((ny, nx // 2 + 1))
        + 1.0j * np.random.standard_normal((ny, nx // 2 + 1))
    ) / np.sqrt(2.0)

    tfft[0, 0] = np.sqrt(2.0) * np.real(tfft[0, 0])
    efft[0, 0] = np.sqrt(2.0) * np.real(efft[0, 0])
    bfft[0, 0] = np.sqrt(2.0) * np.real(bfft[0, 0])

    teb = MapSpectraTEB(
        ffts=[tfft, efft, bfft], map_nx=nx, map_ny=ny, dx=dx, dy=dy
    )
    return map_spec_utils.ClMatrixTEB(tcl).cholesky() * teb


class ForegroundT(object):

    """
    Older T only foregrounds class for generating the Gaussian foregrounds.
    This class exists for interfacing with lens100d results (tutorial purpose).
    We shouldn't use this in 3G.
    """

    def __init__(
        self,
        nx,
        ny,
        dx,
        dy,
        lmax,
        lib_dir,
        sztmpl=None,
        asrc=0.0,
        acib=0.0,
        atsz=0.0,
        pcib=0.8,
        seed=None,
    ):
        self.lmax = lmax
        self.lib_dir = lib_dir
        self.nx, self.ny, self.dx, self.dy = nx, ny, dx, dy
        if sztmpl == None:
            shawcl = np.vstack(
                [
                    [0, 0],
                    [0, 0],
                    np.loadtxt(
                        os.path.dirname(__file__) + "/cl_tsz_148_shaw.dat"
                    ),
                ]
            )[:, 1]
            sztmpl = shawcl / shawcl[3000]

        ls = np.arange(0.0, lmax + 1.0)
        dl2cl = (
            2.0
            * np.pi
            * np.concatenate(
                [
                    [0],
                    1.0
                    / np.arange(1, lmax + 1)
                    / (np.arange(1, lmax + 1) + 1.0),
                ]
            )
        )

        self.ls = ls
        self.clsrc = {
            "L": ls,
            "TT": dl2cl * asrc * (ls / 3000.0) ** 2 * uK ** 2,
        }
        self.clcib = {
            "L": ls,
            "TT": dl2cl * acib * (ls / 3000.0) ** pcib * uK ** 2,
        }
        self.cltsz = {
            "L": ls,
            "TT": dl2cl * atsz * sztmpl[0 : lmax + 1] * uK ** 2,
        }

        self.phas = phas.RandomPhase(
            3 * 6 * nx * (ny // 2 + 1), lib_dir + "/phas", seed=seed
        )

        if True:
            if not os.path.exists(lib_dir):
                os.makedirs(lib_dir)

            if not os.path.exists(lib_dir + "/sim_hash.pk"):
                pk.dump(self.hashdict(), open(lib_dir + "/sim_hash.pk", "wb"))
        utils.hash_check(
            pk.load(open(lib_dir + "/sim_hash.pk", "rb")),
            self.hashdict(),
            ignore=["clsrc", "clcib", "cltsz"],
        )

    def get_clfg(self, lmax=None):
        if lmax == None:
            lmax = self.lmax
        assert lmax <= self.lmax

        return {
            "L": np.arange(lmax + 1),
            "TT": (self.clsrc["TT"] + self.clcib["TT"] + self.cltsz["TT"])[
                0 : lmax + 1
            ],
        }

    def hashdict(self):
        return {
            "clsrc": hashlib.md5(self.clsrc["TT"].view(np.uint8)).hexdigest(),
            "clcib": hashlib.md5(self.clcib["TT"].view(np.uint8)).hexdigest(),
            "cltsz": hashlib.md5(self.cltsz["TT"].view(np.uint8)).hexdigest(),
            "phas": self.phas.hashdict(),
        }

    def get_sim_tqu(self, idx):
        self.phas.set_state(idx)
        teb = tebfft(self.nx, self.ny, self.dx, self.dy, self.clsrc)
        teb += tebfft(self.nx, self.ny, self.dx, self.dy, self.clcib)
        teb += tebfft(self.nx, self.ny, self.dx, self.dy, self.cltsz)
        self.phas.check_state_final(idx)
        ret = teb.get_tqu()
        tmap = np.ascontiguousarray(np.array(ret["T"])[::-1, ::-1])
        np.asarray(ret["T"])[:] = tmap
        qmap = np.ascontiguousarray(np.array(ret["Q"])[::-1, ::-1])
        np.asarray(ret["Q"])[:] = qmap
        umap = np.ascontiguousarray(np.array(ret["U"])[::-1, ::-1])
        np.asarray(ret["U"])[:] = umap
        return ret
