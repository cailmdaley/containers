"""
Script for calculating qcl spectra, amplitudes, and S/N-weighted spectra
Each qecl_analyzer instance handles ONE spectrum type, e.g., k1='Phi_TT', k2='Phi_EB_set'

See qecl_plotter.py for how to use an object from this class.

The parameter file is expected to define the following libraries:
  qecl_len_dd
  qecl_len_ds
  qecl_len_ss

  if considering unlensed sims:
    qecl_len_uu

  if use_n1:
    qecl_len_ss_nofg
    qecl_len_ss2_nofg

  if use_qcr_mc:
    qecl_len_dk

Output:
  None.  the function analyze_qecl() fills the qecl_analyzer object.

  The output of the qecl_analyzer object is now organized as follows:
      self.qcl   # contains un-binned 2D harmonic-space quantities
      self.clpp  # contains binned 1D quantities of phi-phi
      self.clkk  # contains binned 1D quantities of kappa-kappa
      self.amp   # contains single number (over full L-range) quantities


Example call:
  k1='Phi_TT'
  k2='Phi_EB'
  par_file = "/spt3g_software/lensing/specific_analyses/lens100d/par_run.py"
  par    = imp.load_source('par', par_file)
  n0eb   = True
  use_n1 = True
  use_qcr_mc = True
  ampweights = 'analytical_with_lensing_sample_variance' # empirical # analytical_with_lensing_sample_variance
  clear  = True # clear enumerate_progress indicators once finished.
  lmax   = 2000
  mv = anq.qecl_analyzer(par, k1, k2, lmax=lmax, n0eb=n0eb, ampweights=ampweights, use_n1=use_n1, use_qcr_mc=use_qcr_mc, clear=clear, lbins=lbins)
  mv.analyze_qecl()
"""

import os, sys, argparse, imp
import numpy as np
import pickle as pk
from spt3g import lensing as sl
from spt3g.lensing import utils, spectra_stat
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum1D
from spt3g.mapspectra.map_analysis import psd_2d_to_1d
#################################
#################################
class struct(dict):
    """
    Simple extension of the dict type which lets you access
    dictionary elements as if they were data elements of an object,
    instead of using a string inside brackets, as long as the
    dictionary element doesn't have the same name as a method
    of the class dict.
    """

    def __getattr__(self, name):
        """
        Maps instance.key calls into instance['key'], as long as the instance
        doesn't have an existing attribute called 'key'.
        """
        if name in self:
            return self[name]
        else:
            raise AttributeError(
                "Object "
                + typeName(self, short=True)
                + " has no attribute named '"
                + str(name)
                + "'"
            )

    def __getitem__(self, key):
        """
        This will let you use integers to access elements of
        the struct (assuming that the integer isn't already a key).
        NB: Order is arbitrary, and definitely not preserved from any non-dict object of origin!
        """
        try:
            return super(struct, self).__getitem__(key)
        except KeyError:
            try:
                return self.values()[key]
            except TypeError:
                raise KeyError(
                    "Object "
                    + typeName(self, short=True)
                    + " has no attribute named '"
                    + str(key)
                    + "'"
                )

    def __setattr__(self, name, value):
        """
        Inverse of the modifications to __getattr__: If we try to assign
        an attribute with the '.' operator, put it in the dictionary
        instead of making it a new instance attribute.
        """
        if name in self.__dict__:
            self.__dict__[name] = value
        else:
            self[name] = value

    def assignAttribute(self, name, value):
        """
        The normal way to create new attributes, obj.name = value now instead inserts
        "name" into the structure. If you actually want a new instance attribute, use
        this function instead.
        """
        self.__dict__[name] = value

    def setAttribute(self, name, value):
        return self.assignAttribute(name, value)

    def __str__(self, max_linesize=100, alphabetical_sort=True):
        """
        A structure might contain a lot of (possibly large) elements. We want some
        informative, small, display of the contents. This function is inspired by
        the output from IDL's "help" fuction, one of its more useful parts.
        """
        display_string = ""

        if len(self) != 0:
            keys, types, values = [], [], []
            for key in self:
                # Store the name and type of each field.
                keys += [str(key)]
                types += [typeName(self[key], short=True)]

                # Special case for the type field: If this is an arcfile.Register object, then we'd like
                # to know what kind of data it has.
                try:
                    additional_type_info = (
                        ": Holds "
                        + str(self[key].nFrames)
                        + " frames of "
                        + str(self[key].nAxes)
                        + "-axis "
                        + self[key].data_type_name
                        + " data."
                    )
                    types[-1] += additional_type_info
                except (AttributeError, KeyError, TypeError):
                    pass

                # Special case for arrays: Display the data type that they hold.
                if isinstance(self[key], np.ndarray):
                    types[-1] = "<Class 'ndarray'> dtype: " + str(
                        self[key].dtype
                    )

                # Now store some information about the contents, depending on what it is.
                val_string = ""
                try:
                    shape = self[key].shape
                    if len(shape) == 0:
                        val_string += str(
                            self[key]
                        )  # An item with an empty shape tuple is a NumPy scalar. Output the whole thing.
                    else:
                        val_string += "  shape=" + str(
                            self[key].shape
                        )  # If this is an array, show its shape.
                except (AttributeError, KeyError):
                    # except: This entry doesn't have a .shape.
                    try:
                        val_string += self[
                            key
                        ]  # If this entry is a string, display it.
                    except TypeError:
                        # except: This entry isn't a string.
                        try:
                            length = str(len(self[key]))
                            if isinstance(self[key], dict):
                                val_string += (
                                    "  #keys=(" + length + ")"
                                )  # If it's a dictionary, call the length "#keys".
                            else:
                                val_string += (
                                    " length=(" + length + ")"
                                )  # If it's not a string, display a length
                        except (AttributeError, TypeError):
                            # except: This entry doesn't have a __len__.
                            val_string += str(
                                self[key]
                            )  # If it doesn't have a length, just show the str(object)
                values += [val_string]

            # Sort alphabetically by key.
            if alphabetical_sort:
                i_sorted = np.argsort(keys)
            else:
                i_sorted = np.arange(len(keys))
            keys = np.asarray(keys)[i_sorted]
            types = np.asarray(types)[i_sorted]
            values = np.asarray(values)[i_sorted]

            # Figure out the size that we need to allow for the output fields
            key_maxsize = max([len(key) for key in keys])
            type_maxsize = max([len(this_type) for this_type in types])
            value_maxsize = max([len(val) for val in values])

            # Limit the maximum linesize, so that one large field doesn't screw up the whole screen.
            if (
                key_maxsize + type_maxsize + value_maxsize
            ) > max_linesize and max_linesize is not None:
                value_maxsize = np.max(
                    [max_linesize - key_maxsize - type_maxsize, 0]
                )

            # Now create the output string.
            format_string = "{0:>%d} | {1:<%d}|| {2:<%d}\n" % (
                key_maxsize + 1,
                value_maxsize + 5,
                type_maxsize,
            )
            scalars_string, lists_string = "", ""
            for key, this_type, val in zip(keys, types, values):
                # If the string representation of the value is too long, then truncate it.
                if len(val) > max_linesize:
                    val = (
                        val[: max_linesize / 3]
                        + " ... "
                        + val[-max_linesize / 3 :]
                    )

                # Sort the output string, with dictionaries, structures, lists, arrays, etc., last.
                if len(val) > 0 and val[0] == " ":
                    lists_string += format_string.format(key, val, this_type)
            else:
                scalars_string += format_string.format(key, val, this_type)
            display_string = scalars_string + lists_string

        return display_string

    def __repr__(self):
        return self.__str__()

    def __dir__(self):
        """
        Define a custom __dir__ to allow ipython's tab-completion to function.
        This new __dir__ shows the struct's keys as attributes.
        """
        attributes_and_functions = sorted(
            set(
                dir({}) + self.__class__.__dict__.keys() + self.__dict__.keys()
            )
        )
        keys = sorted(self.keys())
        return attributes_and_functions + keys

    def copy(self):
        """
        Override the "copy" method from dict so that it returns a struct.
        """
        return type(self)(super(struct, self).copy())



class qecl_analyzer(object):
    """
    Class for organizing the power spectrum products: clpp, amp
    One qecl_analyzer object handles power spectrum products of one type (k1, k2).

    Inputs:
    par,            pointer to parameter file.  You must imp.load_source() the par file first, then
                    hand the pointer to this class.
    cl_theory_amp,  dictionary with the theory spectrum to compare Amplitudes to
    """

    def __init__(
        self,
        par,
        k1,
        k2,
        lmax=2000,
        n0eb=True,
        simwt=None,
        use_n1=True,
        use_qcr_mc=True,
        clear=True,
        lbins=None,
        ampweights=None,
        cl_theory_amp=None,
    ):
        self.par = par  # pointer to a parameter file object
        self.k1 = k1
        self.k2 = k2
        self.lmax = lmax
        self.n0eb = n0eb
        self.use_n1 = use_n1
        self.use_qcr_mc = use_qcr_mc
        self.clear = clear

        # This object has three levels of "binning", which are separated into three structures.
        self.qcl = struct()  # contains un-binned 2D harmonic-space quantities
        self.clpp = struct()  # contains binned 1D quantities of phi-phi
        self.clkk = struct()  # contains binned 1D quantities of kappa-kappa
        self.amp = (
            struct()
        )  # contains single number (over full L-range) quantities

        ### Specify which sets of sims to use for different pieces of the analysis

        # mean-field
        self.mc_sims_mf = self.par.mc_sims_mf

        # variance
        self.mc_sims_var = self.par.mc_sims_var

        # unlensed variance (default to mc_sims_var)
        if hasattr(self.par, "mc_sims_unl"):
            self.mc_sims_unl = self.par.mc_sims_unl
        else:
            self.mc_sims_unl = self.mc_sims_var

        # phi amplitude correction (default to mc_sims_var)
        if hasattr(self.par, "mc_sims_qcr_mc"):
            self.mc_sims_qcr_mc = self.par.mc_sims_qcr_mc
        elif hasattr(self.par, "mc_sims_tf"):
            print(
                "found OBSOLETE var mc_sims_tf in par file.  change this to mc_sims_qcr_mc."
            )
            self.mc_sims_qcr_mc = self.par.mc_sims_tf
        else:
            self.mc_sims_qcr_mc = self.par.mc_sims_var

        # N0 bias (default to mc_sims_var)
        if hasattr(self.par, "mc_sims_n0"):
            self.mc_sims_n0 = self.par.mc_sims_n0
        else:
            self.mc_sims_n0 = self.par.mc_sims_var

        # N1 bias (default to mc_sims_var)
        if hasattr(self.par, "mc_sims_n1"):
            self.mc_sims_n1 = self.par.mc_sims_n1
        else:
            self.mc_sims_n1 = self.par.mc_sims_var

        if simwt != None:
            print(
                "QECL_ANALZER -- NOTE: simwt argument is deprecated, please switch to ampweights."
            )
            assert ampweights == None
            if simwt == True:
                ampweights = "empirical"
            if simwt == False:
                ampweights = "analytical_without_lensing_sample_variance"

        if ampweights == None:
            ampweights = (
                "analytical_with_lensing_sample_variance"
            )  # Default to this
        self.ampweights = ampweights

        if len(lbins) == 0:
            lbins = np.linspace(10, lmax, 51)
        self.lbins = lbins

        # The theory spectrum to compare amplitudes to
        if cl_theory_amp == None:
            self.cl_theory_amp = self.par.cl_unl
        else:
            self.cl_theory_amp = cl_theory_amp

        self.cache_qcl = {}
        self.cache_ncl = {}

        # Common spectra.
        self.t = lambda l: (l * (l + 1.0)) ** 2 / (2.0 * np.pi) * 1.0e7
        self.qcr = par.qecl_len_dd.get_qcr_lm(self.k1, k2=self.k2)
        self.clpp.cbins = 0.5 * (
            lbins[0:-1] + lbins[1:]
        )  # bin centers - FIX FOR LOG(L)

        # Get the binned theory spectrum.
        # Note, we use the input theory (par.cmbs.cl_unl.clpp), not what was passed in as an argument.
        self.clpp.theory = self.get_cl(
            max(self.par.cl_unl["L"]), self.par.cl_unl["PP"]
        ).rebin(lbins, w=self.t)
        self.clkk.theory = np.array(
            self.clpp.theory * self.t(self.clpp.cbins)
        ).real

        self.clpp.theory_amp = (
            self.get_cl(max(self.cl_theory_amp["L"]), self.cl_theory_amp["PP"])
            .rebin(lbins, w=self.t)
        )
        self.clkk.theory_amp = np.array(
            self.clpp.theory_amp * self.t(self.clpp.cbins)
        ).real

    # -------------------------------------
    # Helper Functions
    # -------------------------------------
    def get_cl(self, lmax, cl):
        cl = MapSpectrum1D(
            np.arange(len(cl) + 1) - 0.5,
            cl,
            dx=self.par.reso_rad,
            dy=self.par.reso_rad,
            map_nx=self.par.npix,
            map_ny=self.par.npix,
        )
        cl_2d = cl.get_2d()
        cl_ret = psd_2d_to_1d(cl_2d, lmax=lmax)
        return cl_ret



    def get_sim_qcl_lm(self, qcllib, idx):
        lib = self.cache_qcl.get(qcllib, None)
        if lib == None:
            lib = {}
            self.cache_qcl[qcllib] = lib
        qcl = lib.get(idx, None)
        if qcl == None:
            qcl = qcllib.get_sim_qcl_lm(self.k1, idx, k2=self.k2)
            lib[idx] = qcl
        return qcl

    def get_sim_ncl_lm(self, qcllib, idx):
        lib = self.cache_ncl.get(qcllib, None)
        if lib == None:
            lib = {}
            self.cache_ncl[qcllib] = lib
        ncl = lib.get(idx, None)
        if ncl == None:
            ncl = qcllib.get_sim_ncl_lm(self.k1, idx, k2=self.k2)
            lib[idx] = ncl
        return ncl

    def get_sim_stats(
        self,
        qcllib,
        mc_sims,
        rdn0=False,
        tscal=1.0,
        lbins=None,
        t=lambda l: 1.0,
        tb=None,
        cov=False,
    ):
        """
        Collect statistics of a collection of sims
        """
        qcr = self.qcr
        if len(lbins) == 0:
            lbins = self.lbins
        if tb == None:
            tb = np.ones(len(lbins) - 1)

        if rdn0 == True:
            ncl = utils.AverageObjects()
            for i, idx in utils.enumerate_progress(
                mc_sims,
                "accumulating ncl, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
                clear=self.clear,
            ):
                ncl.add(self.get_sim_ncl_lm(qcllib, idx))
            ncl_avg = ncl.get()

            ret = spectra_stat.BinnedClStats(lbins, w=t, wb=tb, cov=cov)
            for i, idx in utils.enumerate_progress(
                mc_sims,
                "accumulating qcl, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
                clear=self.clear,
            ):
                ret.add(
                    (
                        self.get_sim_qcl_lm(qcllib, idx)
                        - self.get_sim_ncl_lm(qcllib, idx)
                        + ncl_avg
                    )
                    / qcr
                    * tscal
                )
            return ret
        else:
            ret = spectra_stat.BinnedClStats(lbins, w=t, wb=tb, cov=cov)
            for i, idx in utils.enumerate_progress(
                mc_sims,
                "accumulating qcl, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
                clear=self.clear,
            ):
                ret.add(self.get_sim_qcl_lm(qcllib, idx) / qcr * tscal)
            return ret

    def analyze_qecl(self):
        # self.analyze_clpp()
        # self.analyze_clpp_uu()

        self.analyze_amps()

        # self.print_clpp_amps()
        self.print_amps()

    def analyze_amps(self):
        """ produce l-by-l weighted amplitude estimates for individual bandpower bins, as well as globally from lmin -> lmax. """

        self.amps = {}

        qcr = self.qcr
        lbins = self.lbins
        par = self.par
        t = self.t

        # =====================
        # Initialize this class structure for storing the output of this function
        # =====================

        # Set up the arrays I care about
        nbins = len(lbins) - 1
        nsims_var = len(self.par.mc_sims_var)
        cl_max = 5000

        # sim and data spectra: un-normalized, un-scaled
        self.qcl.sim = {}
        self.qcl.sim["dd"] = np.ndarray([nsims_var, cl_max], dtype=complex)
        self.qcl.sim["uu"] = np.ndarray([nsims_var, cl_max], dtype=complex)
        self.qcl.dat = np.zeros(cl_max, dtype=complex)

        # binned phi-phi spectra
        self.clpp.sim = {
            "dd": np.ndarray(
                [nsims_var, nbins], dtype=float
            ),  # all sim spectra
            "uu": np.ndarray([nsims_var, nbins], dtype=float),
        }
        self.clpp.dat = np.zeros(nbins)  # data spectrum
        self.clpp.n0 = np.zeros(nbins)  # N0 bias spectrum
        self.clpp.n1 = np.zeros(nbins)  # N1 bias spectrum
        self.clpp.rdn0 = np.zeros(nbins)  # RDN0 bias spectrum
        self.clpp.qcr_mc = np.zeros(
            nbins
        )  # quadratic response MC correction to the spectrum

        # binned kappa-kappa spectra.
        # this is what we plot and report
        self.clkk.sim_std = {
            "dd": np.zeros(nbins),  # std of sim spectra, used for error-bars
            "uu": np.zeros(nbins),
        }
        self.clkk.sim_avg = {
            "dd": np.zeros(nbins),  # average of sim spectra
            "uu": np.zeros(nbins),
        }
        self.clkk.dat = np.zeros(nbins)  # data spectrum
        self.clkk.n0 = np.zeros(nbins)  # N0 bias spectrum
        self.clkk.n1 = np.zeros(nbins)  # N1 bias spectrum
        self.clkk.rdn0 = np.zeros(nbins)  # RDN0 bias spectrum

        # amplitudes (over full L-range)
        self.amp.sim = {}
        self.amp.sim["dd"] = np.zeros(
            nsims_var, dtype=float
        )  # all sim amplitudes
        self.amp.sim["uu"] = np.zeros(nsims_var, dtype=float)
        self.amp.sim_std = {
            "dd": float(0.0),  # std of sim amps, used for error-bar
            "uu": float(0.0),
        }
        self.amp.sim_avg = {
            "dd": float(0.0),  # average of sim amps
            "uu": float(0.0),
        }
        self.amp.dat = float(0.0)  # data amplitude
        self.amp.n0 = float(0.0)  # N0 bias
        self.amp.n1 = float(0.0)  # N1 bias
        self.amp.rdn0 = float(0.0)  # RDN0 bias
        self.amp.qcr_mc = float(0.0)  # quadratic response MC correction

        # =====================
        # Configure binning weight functions
        # =====================
        if self.ampweights == "analytical_without_lensing_sample_variance":
            # Use analytical S/N weighting in limit that sample variance is small.
            clpp_lcl_vcl = self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"]) * qcr
            clpp_lcl_vcl_uu = (
                clpp_lcl_vcl
            )  # Treat unlensed sims in the same way

        elif self.ampweights == "analytical_with_lensing_sample_variance":
            # Assume 1/qcr is N0, and then add lensing sample variance.
            clpp_lcl_vcl = self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"])
            clpp_lcl_vcl /= (
                np.sqrt(1.0 / qcr.real)
                + self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"])
            ) ** 2
            clpp_lcl_vcl *= clpp_lcl_vcl.get_num_modes()

            clpp_lcl_vcl_uu = (
                self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"]) * qcr
            )
            clpp_lcl_vcl_uu *= clpp_lcl_vcl_uu.get_num_modes()

        elif self.ampweights == "empirical":
            print("applying simwt")
            # Monte-Carlo based S/N weighting -- determine average variance of Fourier modes in L bins from sims.
            sim_stats_dd = spectra_stat.BinnedClStats(
                np.arange(0, self.qcr.lmax + 1), w=lambda l: 1.0
            )
            for i, idx in utils.enumerate_progress(
                self.mc_sims_var, "accumulating dd var", clear=self.clear
            ):
                sim_stats_dd.add(
                    self.get_sim_qcl_lm(self.par.qecl_len_dd, idx) / qcr
                )
            vcl = sim_stats_dd.var()

            clpp_lcl_vcl = self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"])
            clpp_lcl_vcl[np.nonzero(clpp_lcl_vcl.get_num_modes())] /= vcl[
                np.nonzero(clpp_lcl_vcl.get_num_modes())
            ]

            # Same thing, but for unlensed sims
            sim_stats_uu = spectra_stat.BinnedClStats(
                np.arange(0, self.qcr.lmax + 1), w=lambda l: 1.0
            )
            for i, idx in utils.enumerate_progress(
                self.mc_sims_unl, "accumulating uu var", clear=self.clear
            ):
                sim_stats_uu.add(
                    self.get_sim_qcl_lm(self.par.qecl_len_uu, idx) / qcr
                )
            vcl_uu = sim_stats_uu.var()

            clpp_lcl_vcl_uu = self.get_cl(self.qcr.lmax, self.par.cl_unl["PP"])
            clpp_lcl_vcl_uu[np.nonzero(clpp_lcl_vcl_uu.get_num_modes())] /= vcl_uu[
                np.nonzero(clpp_lcl_vcl_uu.get_num_modes())
            ]

        # Store the weights in the structure for easy access
        self.clpp_lcl_vcl = clpp_lcl_vcl
        self.clpp_lcl_vcl_uu = clpp_lcl_vcl_uu

        # =====================
        # Calculate bias terms
        # =====================
        # get average sim spectra
        self.sim_stats_ds = self.get_sim_stats(
            self.par.qecl_len_ds,
            self.mc_sims_n0,
            rdn0=False,
            lbins=np.arange(0, cl_max + 1),
        )
        self.sim_stats_ss = self.get_sim_stats(
            self.par.qecl_len_ss,
            self.mc_sims_n0,
            rdn0=False,
            lbins=np.arange(0, cl_max + 1),
        )
        if self.use_n1:
            self.sim_stats_ss_nofg = self.get_sim_stats(
                self.par.qecl_len_ss_nofg,
                self.mc_sims_n1,
                rdn0=False,
                lbins=np.arange(0, cl_max + 1),
            )
            self.sim_stats_ss2_nofg = self.get_sim_stats(
                self.par.qecl_len_ss2_nofg,
                self.mc_sims_n1,
                rdn0=False,
                lbins=np.arange(0, cl_max + 1),
            )

        # estimate n0 biases
        self.n0_bias = self.sim_stats_ss.avg() * 2.0  # simulation only n0 bias
        self.rdn0 = (
            self.sim_stats_ds.avg() * 4.0 - self.n0_bias
        )  # data-dependent n0 bias.

        # estimate n1 biases
        if self.use_n1:
            self.n1_bias = (
                self.sim_stats_ss2_nofg.avg() - self.sim_stats_ss_nofg.avg()
            ) * 2.0
        else:
            self.n1_bias = self.sim_stats_ss.avg()
            self.n1_bias *= 0.0  # if this is not used, set N1 to zero

        # for qcr_mc
        if self.use_qcr_mc:
            self.sim_stats_dk = self.get_sim_stats(
                par.qecl_len_dk,
                self.mc_sims_qcr_mc,
                lbins=np.arange(0, cl_max + 1),
                rdn0=False,
            )

        # =====================
        # Calculate lensing amplitudes, over full l range as well as individual bins.
        # Note: "zip(lbins[0:-1], lbins[1:])" -> individual bins
        #       "[[lbins[0], lbins[-1]]" -> full l range
        # =====================

        # Pre-calculate the sim spectra over the full 2D L-plane

        # estimate average analytical noise spectrum
        if self.n0eb:
            ncl = utils.AverageObjects()
            for i, isim in utils.enumerate_progress(
                self.mc_sims_var,
                "accumulating ncl, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
                clear=self.clear,
            ):
                ncl.add(self.get_sim_ncl_lm(self.par.qecl_len_dd, isim))
            ncl_avg = ncl.get()

        # estimate dd sim amplitudes
        for i, isim in utils.enumerate_progress(
            self.mc_sims_var,
            "accumulating amp_dd, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
            clear=self.clear,
        ):
            tmp_dd = (
                self.get_sim_qcl_lm(self.par.qecl_len_dd, isim) / self.qcr.real
                - self.n0_bias
            )

            if (
                self.n0eb
            ):  # subtract analytical estimate of sim's rdn0, if requested.
                tmp_dd -= (
                    self.get_sim_ncl_lm(self.par.qecl_len_dd, isim) - ncl_avg
                ) / self.qcr.real

            if self.use_n1:  # subtract N1, if requested
                tmp_dd = tmp_dd - self.n1_bias

            self.qcl.sim["dd"][i, :] = tmp_dd.copy()

        for i, isim in utils.enumerate_progress(
            self.mc_sims_unl,
            "accumulating amp_uu, (self.k1,k2)=(%s,%s)" % (self.k1, self.k2),
            clear=self.clear,
        ):
            self.qcl.sim["uu"][i, :] = (
                self.get_sim_qcl_lm(self.par.qecl_len_uu, isim) / self.qcr.real
                - self.n0_bias
            )

        # estimate data amplitude
        self.qcl.dat = (
            par.qecl_len_dd.get_dat_qcl_lm(self.k1, k2=self.k2) / qcr.real
            - self.rdn0
        )

        if self.use_n1:  # subtract N1, if requested
            self.qcl.dat = self.qcl.dat - self.n1_bias

        # Get amps for each bin
        for ibin in range(0, nbins + 1):

            if ibin == nbins:  # full L-range
                ibin = -1

            self.calculate_amp_in_one_bin(ibin)

        # Calculate std and avg of sims, for convenience
        self.amp.sim_avg["dd"] = np.real(np.average(self.amp.sim["dd"]))
        self.amp.sim_std["dd"] = np.real(np.std(self.amp.sim["dd"]))
        self.amp.sim_avg["uu"] = np.real(np.average(self.amp.sim["uu"]))
        self.amp.sim_std["uu"] = np.real(np.std(self.amp.sim["uu"]))

        # Now convert spectra to kappa: clkk
        # These are the spectra we plot and report
        for i in range(0, len(lbins) - 1):
            tl = self.clpp.cbins[i]
            self.clkk.sim_std["dd"][i] = (
                np.std(self.clpp.sim["dd"][:, i])
                * self.clpp.theory_amp[i].real
                * t(tl)
            )
            self.clkk.sim_avg["dd"][i] = (
                np.average(self.clpp.sim["dd"][:, i])
                * self.clpp.theory_amp[i].real
                * t(tl)
            )
            self.clkk.sim_std["uu"][i] = (
                np.std(self.clpp.sim["uu"][:, i])
                * self.clpp.theory_amp[i].real
                * t(tl)
            )
            self.clkk.sim_avg["uu"][i] = (
                np.average(self.clpp.sim["uu"][:, i])
                * self.clpp.theory_amp[i].real
                * t(tl)
            )

        self.clkk.dat = (
            self.clpp.dat * np.real(self.clpp.theory_amp) * t(self.clpp.cbins)
        )
        self.clkk.n0 = (
            self.clpp.n0 * np.real(self.clpp.theory_amp) * t(self.clpp.cbins)
        )
        self.clkk.n1 = (
            self.clpp.n1 * np.real(self.clpp.theory_amp) * t(self.clpp.cbins)
        )
        self.clkk.rdn0 = (
            self.clpp.rdn0 * np.real(self.clpp.theory_amp) * t(self.clpp.cbins)
        )
        self.clkk.cbins = self.clpp.cbins  # lame, I know, but less confusing

    def calculate_amp_in_one_bin(self, ibin):
        """
        Helper function to calculate the amplitude within one bin.
        INPUTS:
          ibin         : [Int] bin number as defined in self.lbins.
                               if ibin=-1, then consider the full L-range
        """

        # --------------------
        # Set the data pointers to the correct output structure
        # --------------------

        sim_amps = {
            "dd": np.zeros(len(self.par.mc_sims_var), dtype=float),
            "uu": np.zeros(len(self.par.mc_sims_var), dtype=float),
        }
        # Use this case to calculate amplitudes (over full L-range)
        if ibin == -1:
            lmin_bin = self.lbins[0]
            lmax_bin = self.lbins[-1]

        # Use this case to calculate binned spectra
        else:
            lmin_bin = self.lbins[ibin]
            lmax_bin = self.lbins[ibin + 1]

        # DEBUGGING - use self.clpp.tl[imin]
        # tl =  0.5*(lmin_bin + lmax_bin) # bin centers

        print("estimating amps for bin [%4d, %4d]" % (lmin_bin, lmax_bin))

        # prepare to estimate bin amplitudes.
        print(lmin_bin, lmax_bin)
        lmin_bin = int(lmin_bin)
        lmax_bin = int(lmax_bin)
        nmd = self.qcr.get_num_modes()[lmin_bin:lmax_bin]  # number of modes

        # Denominators for the amplitude calculation.
        # den    = np.sum( (self.cl_theory_amp.clpp[lmin_bin:lmax_bin] * self.clpp_lcl_vcl[lmin_bin:lmax_bin])[np.nonzero(nmd)] )
        # den_uu = np.sum( (self.cl_theory_amp.clpp[lmin_bin:lmax_bin] * self.clpp_lcl_vcl_uu[lmin_bin:lmax_bin])[np.nonzero(nmd)] )

        # To avoid floor problem, cl_theory_amp.clpp needs the same processing as sims
        clpp_theory = self.get_cl(self.lmax, self.cl_theory_amp["PP"])
        den = np.sum(
            (
                clpp_theory[lmin_bin:lmax_bin]
                * self.clpp_lcl_vcl[lmin_bin:lmax_bin]
            )[np.nonzero(nmd)]
        )
        den_uu = np.sum(
            (
                clpp_theory[lmin_bin:lmax_bin]
                * self.clpp_lcl_vcl_uu[lmin_bin:lmax_bin]
            )[np.nonzero(nmd)]
        )

        def get_amp(qcl):
            num = np.sum(
                (
                    qcl[lmin_bin:lmax_bin]
                    * self.clpp_lcl_vcl[lmin_bin:lmax_bin]
                )[np.nonzero(nmd)]
            )
            return np.real(num / den)

        def get_amp_uu(qcl):
            num = np.sum(
                (
                    qcl[lmin_bin:lmax_bin]
                    * self.clpp_lcl_vcl_uu[lmin_bin:lmax_bin]
                )[np.nonzero(nmd)]
            )
            return np.real(num / den_uu)

        # --

        # estimate phi transfer function
        # re-name this to "qcr_mc"
        if self.use_qcr_mc:
            tmp_dd = self.sim_stats_dk.avg()

            den_qcr_mc = np.sum(
                (
                    self.par.cl_unl["PP"][lmin_bin:lmax_bin]
                    * self.clpp_lcl_vcl[lmin_bin:lmax_bin]
                )[np.nonzero(nmd)]
            )
            num = np.sum(
                (
                    tmp_dd[lmin_bin:lmax_bin]
                    * self.clpp_lcl_vcl[lmin_bin:lmax_bin]
                )[np.nonzero(nmd)]
            )
            qcr_mc_amp = num / den_qcr_mc

            qcr_mc = (
                1.0 / qcr_mc_amp.real ** 2
            )  # store qcr_mc in output structure
        else:
            qcr_mc = 1.0

        # estimate dd sim amplitudes
        for i, isim in utils.enumerate_progress(
            self.mc_sims_var,
            "accumulating amp, (k1,k2)=(%s,%s)" % (self.k1, self.k2),
            clear=self.clear,
        ):
            tmp_dd = self.qcl.sim["dd"][i].copy()

            if self.use_qcr_mc:  # apply bin-specific QCR_MC, if requested
                tmp_dd = tmp_dd * qcr_mc

            sim_amps["dd"][i] = get_amp(tmp_dd)

        # estimate uu sim amplitudes
        for i, isim in utils.enumerate_progress(
            self.mc_sims_unl,
            "accumulating amp, (self.k1,k2)=(%s,%s)" % (self.k1, self.k2),
            clear=self.clear,
        ):
            tmp_uu = self.qcl.sim["uu"][i].copy()

            if self.use_qcr_mc:  # apply bin-specific QCR_MC, if requested
                tmp_uu = tmp_uu * qcr_mc

            sim_amps["uu"][i] = get_amp_uu(tmp_uu)

        tmp_dat = self.qcl.dat.copy()
        if self.use_qcr_mc:  # apply bin-specific QCR_MC, if requested
            tmp_dat = tmp_dat * qcr_mc

        dat_amps = get_amp(tmp_dat)

        n1 = get_amp(self.n1_bias * qcr_mc)
        n0 = get_amp(self.n0_bias * qcr_mc)
        rdn0 = get_amp(self.rdn0 * qcr_mc)

        # Store output to the correct output structure
        # Use this case to calculate amplitudes (over full L-range)
        if ibin == -1:
            self.amp.sim = sim_amps
            self.amp.qcr_mc = qcr_mc
            self.amp.dat = dat_amps
            self.amp.n1 = n1
            self.amp.n0 = n0
            self.amp.rdn0 = rdn0

        # Use this case to calculate binned spectra
        else:
            self.clpp.sim["dd"][:, ibin] = sim_amps["dd"]
            self.clpp.sim["uu"][:, ibin] = sim_amps["uu"]
            self.clpp.qcr_mc[ibin] = qcr_mc
            self.clpp.dat[ibin] = dat_amps
            self.clpp.n1[ibin] = n1
            self.clpp.n0[ibin] = n0
            self.clpp.rdn0[ibin] = rdn0

    def print_amps(self):
        """ print overall amplitude estimate obtained with l-by-l weighting. """

        if not hasattr(self, "amp"):
            self.analyze_amps()

        print("(k1,k2)=(%s,%s), amps" % (self.k1, self.k2))
        print(
            "        sim_unl: %5.3f +- %5.3f"
            % (self.amp.sim_avg["uu"], self.amp.sim_std["uu"])
        )
        print(
            "        sim_len: %5.4f +- %5.4f (expect %5.2f sigma)"
            % (
                self.amp.sim_avg["dd"],
                self.amp.sim_std["dd"],
                self.amp.sim_avg["dd"] / self.amp.sim_std["dd"],
            )
        )
        print(
            "            dat: %5.4f (%5.4f sigma)"
            % (self.amp.dat, (self.amp.dat) / self.amp.sim_std["dd"])
        )

        # print "(k1,k2)=(%s,%s), amps"%(self.k1,self.k2)
        # print "        sim_unl: %5.3f +- %5.3f"%(np.real(np.average(self.amp.sim['uu'])), np.std(self.amp.sim['uu']))
        # print "        sim_len: %5.4f +- %5.4f (expect %5.2f sigma)"%(np.real(np.average(self.amp.sim['dd'])),
        #                                                                   np.std(self.amp.sim['dd']),
        #                                                                   (np.real((np.average(self.amp.sim['dd']))))/np.std(self.amp.sim['dd']))
        # print "            dat: %5.4f (%5.4f sigma)"%(self.amp.dat, (self.amp.dat)/np.std(self.amp.sim['dd']))

    def print_chisq(self):
        """
        print the chisq of the spectrum relative to sims
        """
        from scipy.stats import chisqprob

        if hasattr(self, "amp") == False:
            self.analyze_amps()

        # scale the theory spectrum?
        best_fit = True
        if best_fit:
            scaled_amp = self.amp.dat
            scaled_amp_sim = self.amp.sim_avg["dd"]
        else:
            scaled_amp = 1
            scaled_amp_sim = 1

        bth_cl = (
            self.clkk.theory_amp
        )  # np.array(self.clpp.theory * self.t(tl)).real

        # Get data chisq
        chisq_dat = np.sum(
            ((self.clkk.dat - bth_cl * scaled_amp) / self.clkk.sim_std["dd"])
            ** 2.0
        )
        # chisq_dat = np.sum((((np.array(self.amp_clpp_dat)).real - bth_cl*scaled_amp) / (np.array(self.amp_clpp_std)).real)**2.)

        # --------------------
        # PTE method 1
        pte_1 = chisqprob(chisq_dat, len(self.lbins) - 1)

        # --------------------
        # PTE method 2
        # We use this one!
        chisq_sims = []
        # bins = self.amps.keys()[1:]
        nsims_var = len(self.par.mc_sims_var)
        for isim in range(nsims_var):

            sim_spec = (
                self.clpp.sim["dd"][isim, :]
                * np.real(self.clpp.theory_amp)
                * self.t(self.clpp.cbins)
            )
            # sim_spec = (self.amp.sim['dd'][isim,:] * self.clpp.theory_amp * self.t(self.clpp.cbins)).real

            # calculate chisq
            chsq = np.sum(
                (
                    (sim_spec - bth_cl * scaled_amp_sim)
                    / (np.array(self.clkk.sim_std["dd"])).real
                )
                ** 2.0
            )
            chisq_sims.append(chsq)

        # calculate the PTE
        pte_dat = len(np.where(chisq_sims > chisq_dat)[0]) / float(
            len(chisq_sims)
        )

        # Store the value
        self.clkk.chisq_dat = chisq_dat
        self.clkk.pte_dat = pte_dat
        self.clkk.chisq_sims = {"dd": chisq_sims}

        # --------------------
        # Print results
        print("PTE (sim hist)  = ", pte_dat)
        print("PTE (chisqprob) = ", pte_1)
