import os
import sys
import numpy as np
import pylab as pl
import pickle as pk
import glob
import hashlib
import warnings
import healpy as hp
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum2D, MapSpectrum1D
from spt3g.mapspectra.map_analysis import mapffts_to_powerspectra
from . import utils

class BinnedClStats(object):
    """
    Gather spectrum statistics of MapSpectrum1D objects

    Parameters
    ------------
    lbins: numpy-array
        array definining the bin (lower-edge)
    w: function
        weight function
    wb: numpy array
        weight for each lbin
    bcl_av_sum: numpy array
        sum of binned Cl spectrum
    bcl_sq_sum: numpy array
        Sum of binned Cl^2 spectrum
    bcl_cv_sum: numpy array
        Sum of covariance
    """

    def __init__(self, lbins, w=lambda l: 1.0, wb=None, cov=False):
        if (wb == None).all():
            wb = np.ones(len(lbins) - 1)

        self.lbins = lbins
        self.w = w
        self.wb = wb
        self.docov = cov

        self.nsims = 0
        self.bcl_av_sum = np.zeros(len(lbins) - 1)
        self.bcl_sq_sum = np.zeros(len(lbins) - 1)

        if cov == True:
            self.bcl_cv_sum = np.zeros((len(lbins) - 1, len(lbins) - 1))

    def add(self, obj):
        tml = obj.rebin(self.lbins, w=self.w).real * self.wb

        self.nsims += 1
        self.bcl_av_sum += tml
        self.bcl_sq_sum += tml ** 2

        if self.docov == True:
            self.bcl_cv_sum += np.outer(tml, tml)

    def add_cl(self, obj):
        tml = next(iter(mapffts_to_powerspectra(obj,
            lbins=self.lbins,  ell_weights_2d=self.w).values())) * self.wb
        self.nsims += 1
        self.bcl_av_sum += tml
        self.bcl_sq_sum += tml ** 2

        if self.docov == True:
            self.bcl_cv_sum += np.outer(tml, tml)

    def avg(self):
        return self.bcl_av_sum / self.nsims

    def std(self):
        return np.sqrt((self.bcl_sq_sum - self.bcl_av_sum ** 2 / self.nsims)/ self.nsims)

    def var(self):
        return ((self.bcl_sq_sum - self.bcl_av_sum ** 2 / self.nsims)/ self.nsims)

    def cov(self):
        assert self.docov == True
        return (
            self.bcl_cv_sum
            - np.outer(self.bcl_av_sum, self.bcl_av_sum) / self.nsims
            ) / self.nsims

    def save(self, fname):
        """
        Turn this object into a dictionary so it can be Pickled.
        Note: lambda functions can't be pickled, so you have to
        pass in the same lambda function on read-in.
        """
        ml_stats_dict = {
            "lbins": self.lbins,
            "wb": self.wb,
            "docov": self.docov,
            "nsims": self.nsims,
            "bcl_av_sum": self.bcl_av_sum,
            "bcl_sq_sum": self.bcl_sq_sum,
        }
        if hasattr(self, "bcl_cv_sum"):
            ml_stats_dict["bcl_cv_sum"] = self.bcl_cv_sum

        if not os.path.exists(fname):
            pk.dump(
                ml_stats_dict, open(fname, "wb"), protocol=pk.HIGHEST_PROTOCOL
            )
        else:
            raise IOError("File already exists!")

    def load(self, fname):
        """
        Load an object from a saved pickle file.
        Note: you must pass in the same lambda function
        for w as was used for the original object
        """
        ml_stats_dict = pk.load(open(fname, "rb"))
        self.lbins = ml_stats_dict["lbins"]
        self.wb = ml_stats_dict["wb"]
        self.docov = ml_stats_dict["docov"]
        self.nsims = ml_stats_dict["nsims"]
        self.bcl_av_sum = ml_stats_dict["bcl_av_sum"]
        self.bcl_sq_sum = ml_stats_dict["bcl_sq_sum"]
        if "bcl_cv_sum" in list(ml_stats_dict.keys()):
            self.bcl_cv_sum = ml_stats_dict["bcl_cv_sum"]

    def plot_fill_between(self, t=lambda l: 1.0, m=None, **kwargs):
        """
        plot this spectrum as a histogram

        Parameters
        --------------
        t: function
            l-dependent coefficient.  e.g.,
            t = lambda l : (l*(l+1.))**2 / (2.*np.pi) * 1.e7
        m: np array
            the spectrum to be plotted
        """
        ls = 0.5 * (self.lbins[0:-1] + self.lbins[1:])  # bin centers
        if (m == None).all():
            m = self.avg()

        pl.fill_between(
            ls,
            t(ls) * (m - self.std()),
            t(ls) * (m + self.std()),
            **kwargs
        )

    def plot_error_bars(
        self, t=lambda l: 1.0, m=None, p=pl.errorbar, **kwargs
    ):
        ls = 0.5 * (self.lbins[0:-1] + self.lbins[1:])  # bin centers
        if (m == None).all():
            m = self.avg()

        p(ls, t(ls) * m, yerr=(t(ls) * self.std()), **kwargs)


class BinnedClCrossStats(object):
    """
    Class to gater stats from a cross-spectrum.
    This emulates BinnedClStats(), but starting from MapSpectrum2Ds
    of the two objects that you want to cross instead of MapSpectrum1Ds
    """

    def __init__(self, lbins, w=lambda l: 1.0, wb=None, cov=False):
        if (wb == None).all():
            wb = np.ones(len(lbins) - 1)

        self.lbins = lbins
        self.w = w
        self.wb = wb
        self.docov = cov

        self.nsims = 0
        self.bcl_av_sum = np.zeros(len(lbins) - 1)
        self.bcl_sq_sum = np.zeros(len(lbins) - 1)

        if cov == True:
            self.bcl_cv_sum = np.zeros((len(lbins) - 1, len(lbins) - 1))

    def add(self, obj1, obj2):
        tml = (
            next(iter(mapffts_to_powerspectra(input1=obj1, lbins=self.lbins,
                input2=obj2, ell_weights_2d=self.w).values()))
            * self.wb
        )

        self.nsims += 1
        self.bcl_av_sum += tml
        self.bcl_sq_sum += tml ** 2

        if self.docov == True:
            self.bcl_cv_sum += np.outer(tml, tml)

    def add_cl(self, obj1, obj2):
        tml = (
            next(iter(mapffts_to_powerspectra(input1=obj1, lbins=self.lbins,
                input2=obj2, ell_weights_2d=self.w).values()))
            * self.wb
        )

        self.nsims += 1
        self.bcl_av_sum += tml
        self.bcl_sq_sum += tml ** 2

        if self.docov == True:
            self.bcl_cv_sum += np.outer(tml, tml)

    def avg(self):
        return self.bcl_av_sum / self.nsims

    def std(self):
        return np.sqrt(
                    (self.bcl_sq_sum - self.bcl_av_sum ** 2 / self.nsims)
                    / self.nsims
                )


    def var(self):
        return ((self.bcl_sq_sum - self.bcl_av_sum ** 2 / self.nsims)
                / self.nsims)

    def save(self, fname):
        """
        Turn this object into a dictionary so it can be Pickled.
        Note: lambda functions can't be pickled,
        so you have to pass in the same lambda function on read-in.
        """
        ml_stats_dict = {
            "lbins": self.lbins,
            "wb": self.wb,
            "docov": self.docov,
            "nsims": self.nsims,
            "bcl_av_sum": self.bcl_av_sum,
            "bcl_sq_sum": self.bcl_sq_sum,
        }
        if hasattr(self, "bcl_cv_sum"):
            ml_stats_dict["bcl_cv_sum"] = self.bcl_cv_sum

        if not os.path.exists(fname):
            pk.dump(ml_stats_dict, open(fname, "wb"))
        else:
            raise IOError("File already exists!")

    def load(self, fname):
        """
        Load an object from a saved pickle file.
        Note: you must pass in the same lambda
        function for w as was used for the original object
        """
        ml_stats_dict = pk.load(open(fname, "rb"))
        self.lbins = ml_stats_dict["lbins"]
        self.wb = ml_stats_dict["wb"]
        self.docov = ml_stats_dict["docov"]
        self.nsims = ml_stats_dict["nsims"]
        self.bcl_av_sum = ml_stats_dict["bcl_av_sum"]
        self.bcl_sq_sum = ml_stats_dict["bcl_sq_sum"]
        if "bcl_cv_sum" in list(ml_stats_dict.keys()):
            self.bcl_cv_sum = ml_stats_dict["bcl_cv_sum"]

    def plot_fill_between(self, t=lambda l: 1.0, m=None, **kwargs):
        """
        plot this spectrum as a histogram

        Parameters
        --------------
        t: function
            l-dependent coefficient.  e.g.,
            t = lambda l : (l*(l+1.))**2 / (2.*np.pi) * 1.e7
        m: np array
            the spectrum to be plotted
        """
        ls = 0.5 * (self.lbins[0:-1] + self.lbins[1:])  # bin centers
        if (m == None).all():
            m = self.avg()

        pl.fill_between(
            ls,
            t(ls) * (m - self.std()),
            t(ls) * (m + self.std()),
            **kwargs
        )

    def plot_error_bars(
        self, t=lambda l: 1.0, m=None, p=pl.errorbar, **kwargs
    ):
        ls = 0.5 * (self.lbins[0:-1] + self.lbins[1:])  # bin centers
        if (m == None).all():
            m = self.avg()

        p(ls, t(ls) * m, yerr=(t(ls) * self.std()), **kwargs)


class ClStats(object):
    """
    Class to gater stats from a auto-spectrum.
    This emulates BinnedClStats(), but starting from MapSpectrum2Ds
    instead of Mapspectrum1Ds
    """
    def __init__(self, lbins, w=None, wb=None):
        self.lbins = lbins
        self.w = w
        self.wb = wb
        self.nsims = 0
        self.bcl_av_sum = utils.SumObjects()
        self.bcl_sq_sum = utils.SumObjects()

    def add(self, obj):
        tml = next(iter(mapffts_to_powerspectra(obj,
            lbins=self.lbins,  ell_weights_2d=self.w).values()))
        if self.wb != None:
            tml *= self.wb

        self.nsims += 1
        self.bcl_av_sum += tml
        self.bcl_sq_sum += tml * tml

    def avg(self):
        return self.bcl_av_sum.get() / self.nsims

    def std(self):
        return np.sqrt(
                    (self.bcl_sq_sum - self.bcl_av_sum ** 2 / self.nsims)
                    / self.nsims
                )

    def var(self):
        return (
            self.bcl_sq_sum.get()
            - self.bcl_av_sum.get() * self.bcl_av_sum.get() / self.nsims
        ) / self.nsims

    def std(self):
        ret = self.var()
        ret = np.sqrt(ret)
        return ret
