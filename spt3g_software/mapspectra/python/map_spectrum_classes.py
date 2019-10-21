import copy, hashlib
import numpy as np
from spt3g import core, mapmaker, coordinateutils
import warnings

rad = core.G3Units.rad

class MapSpectrum2D(np.ndarray):
    """ 2D map spectrum object.

    Selected attributes
    ---------
    map_nx: int
        number of pixels in the x direction
    dx: float
        size of pixels in the x direction [units of G3Unit (radians)]
    is_real: boolean
        Whether the mapspectrum corresponds to np.fft.rfft2 on a real map
    """
    def __new__(cls, map_nx, dx, spec=None, map_ny=None, dy=None,
                units=getattr(core.G3TimestreamUnits, "None"),
                proj=getattr(coordinateutils.MapProjection, "ProjNone")):
        if map_ny is None:
            map_ny = map_nx
        if dy is None:
            dy = dx
        if spec is None:
            spec = np.zeros((map_ny, map_nx), dtype=np.complex)
            is_real = False
        obj = np.asarray(spec).view(cls)
        obj.map_nx, obj.map_ny, obj.dx, obj.dy = map_nx, map_ny, dx, dy
        obj.units, obj.proj = units, proj
        if (map_ny, map_nx) == spec.shape:
            obj.is_real = False
        elif (map_ny, map_nx // 2 + 1) == spec.shape:
            obj.is_real = True
        else:
            raise TypeError("Input spectrum does not match the map size")
        return obj

    def hashdict(self):
        """ Uniquely characterize the object.

        Returns a dictionary which should uniquely
        characterize the contents of this object.
        """
        return {
            "pix": {
                "map_nx": self.map_nx,
                "dx": self.dx,
                "map_ny": self.map_ny,
                "dy": self.dy,
            },
            "spec": hashlib.sha1(self.view(np.uint8)).hexdigest(),
        }

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        outputs = kwargs.get('out', ())
        # check that data type is compatible
        for x in inputs + outputs:
            if isinstance(x, self.__class__):
                if not self.compatible(x):
                    raise TypeError("Shape or resolution not compatible! ")
            elif np.isscalar(x):
                pass
            elif isinstance(x, (np.ndarray, list)) and np.shape(x) == self.shape:
                pass
            elif isinstance(x, (np.ndarray, list)) and np.shape(x) != self.shape:
                raise ValueError("operands could not be broadcast together with shapes"
                                 + str(np.shape(x)) + str(self.shape))
            else:
                return NotImplemented
        # gather input arguments for ufunc
        args = []
        for i, input in enumerate(inputs):
            if isinstance(input, self.__class__):
                args.append(input.view(np.ndarray))
            else:
                args.append(input)
        if outputs:
            kwargs['out'] = tuple(
                x.view(np.ndarray) if isinstance(x, self.__class__) else x
                for x in outputs)
        # do the calculation
        results = super().__array_ufunc__(ufunc, method,
                                          *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented
        if method == 'at':
            return None
        if ufunc.nout == 1:
            results = (results,)
        results = (np.nan_to_num(result) for result in results)
        # return the object with the same attributes 
        results = [result.view(self.__class__) for result in results]
        for result in results:
            if isinstance(result, self.__class__):
                for attr in ["map_nx", "dx", "map_ny",
                             "dy", "units", "proj", "is_real"]:
                    setattr(result, attr, getattr(self, attr, None))
        results=tuple(results)
        return results[0] if len(results) == 1 else results

    def __array_finalize__(self, obj):
        if obj is None: return
        # same attributes when create new from template
        for attr in ["map_nx", "dx", "map_ny",
                     "dy", "units", "proj", "is_real"]:
            if not hasattr(self, attr):
                setattr(self, attr, getattr(obj, attr, None))

    def compatible(self, other):
        """ Check that two object are compatible

        Check whether this map can be added, subtracted, etc.
        to the map 'other'.
        """
        return (
            isinstance(other, self.__class__)
            and other.shape == self.shape
            and (self.map_nx == other.map_nx)
            and (self.map_ny == other.map_ny)
            and np.allclose(self.dx, other.dx, rtol=1.0e-8, atol=1.0e-8)
            and np.allclose(self.dy, other.dy, rtol=1.0e-8, atol=1.0e-8)
        )

    def copy(self):
        """ return a clone of this object. """
        return MapSpectrum2D(
            self.map_nx,
            self.dx,
            self.view(np.ndarray).copy(),
            map_ny=self.map_ny,
            dy=self.dy,
            units=self.units,
            proj=self.proj
        )

    def get_lxly(self):
        """ returns the (lx, ly) pair for each Fourier mode. """
        if self.is_real:
            return np.meshgrid(
                np.fft.fftfreq(self.map_nx, self.dx / rad)[
                    0 : self.map_nx // 2 + 1
                ]
                * 2.0
                * np.pi,
                np.fft.fftfreq(self.map_ny, self.dy / rad) * 2.0 * np.pi,
            )
        else:
            return np.meshgrid(
                np.fft.fftfreq(self.map_nx, self.dx / rad) * 2.0 * np.pi,
                np.fft.fftfreq(self.map_ny, self.dy / rad) * 2.0 * np.pi,
            )

    def get_ell(self):
        """ returns l = sqrt(lx**2 + ly**2) for each Fourier mode
        """
        lx, ly = self.get_lxly()
        return np.hypot(lx, ly)

    def get_pixel_window(self):
        """ return the pixel window transfer function of this object.
        """
        lx, ly = self.get_lxly()
        transf = self.copy()
        transf[0, 0] = 1.0
        transf[0, 1:] = np.sin(self.dx / rad * lx[0, 1:] / 2.0) / (
            self.dx / rad * lx[0, 1:] / 2.0
        )
        transf[1:, 0] = np.sin(self.dy / rad * ly[1:, 0] / 2.0) / (
            self.dy / rad * ly[1:, 0] / 2.0
        )
        transf[1:, 1:] = (
            np.sin(self.dx / rad * lx[1:, 1:] / 2.0)
            * np.sin(self.dy / rad * ly[1:, 1:] / 2.0)
            / (self.dx / rad * self.dy / rad * lx[1:, 1:] * ly[1:, 1:] / 4.0)
        )
        return transf

    def get_real(self):
        """ Shrink the size of self by a factor of two for real maps

        This converts the return of np.fft.fft2 to the return of np.fft.rfft2.

        Return:
        ------
        MapSpectrum2D object
            The return corresponds to the real part of the map associated
            with this fft, which is half the size of the original map.
        """
        if self.is_real:
            return self
        else:
            if self.map_nx%2 == 0:
                half_ind = [self.map_nx // 2 + 1, self.map_nx // 2]
            else:
                half_ind = [self.map_nx // 2 + 1, self.map_nx // 2+1]
            if not np.allclose(
                np.array(self[1:, half_ind[0] :]),
                np.conj(np.array(self[1:, 1 : half_ind[1]][::-1, ::-1])),
                rtol=1e-8, atol=1e-8
            ) or not np.allclose(
                np.array(self[0, half_ind[0] :]),
                np.conj(np.array(self[0, 1 : half_ind[1]][::-1])),
                rtol=1e-8, atol=1e-8
            ):
                warnings.warn("Does not correspond to a real map!")
            # return the part that corresponds to the real part of the map
            return MapSpectrum2D(
                self.map_nx, self.dx,
                spec=self[:, 0 : half_ind[0]],
                map_ny=self.map_ny, dy=self.dy,
                units=self.units,
                proj=self.proj
            )

    def get_complex(self):
        """ Expand the size of self by a factor of two for real maps
        This converts the return of np.fft.rfft2 to the return of np.fft.fft2.

        Return:
        ------
        MapSpectrum2D object
            The return is the same size as the original map
            instead of half the size (np.fft.rfft2).
        """
        if self.is_real:
            if self.map_nx%2 == 0:
                half_ind = [self.map_nx // 2 + 1, self.map_nx // 2]
            else:
                half_ind = [self.map_nx // 2 + 1, self.map_nx // 2+1]
            spec = np.zeros((self.map_ny, self.map_nx), dtype=np.complex)
            spec[:, 0 : half_ind[0]] = self
            spec[0, half_ind[0] :] = np.conj(
                np.array(self[0, 1 : half_ind[1]][::-1])
            )
            spec[1:, half_ind[0] :] = np.conj(
                np.array(self[1:, 1 : half_ind[1]][::-1, ::-1])
            )
            return MapSpectrum2D(
                self.map_nx, self.dx, spec=spec,
                map_ny=self.map_ny, dy=self.dy,
                units=self.units,
                proj=self.proj
            )
        else:
            return self

    def get_rmap(self, pol_type="T"):
        """ return the real map given by this FFT. """
        tfac = np.sqrt((self.map_nx * self.map_ny)
                       / (self.dx / rad * self.dy / rad))
        if self.is_real:
            rmap = np.fft.irfft2(self, s=(self.map_ny, self.map_nx)) * tfac
        else:
            cmap = np.fft.ifft2(self) * tfac
            assert np.allclose(
                cmap.imag, 0, rtol=1e-8, atol=1e-8
            ), "Does not correspond to a real map!"
            rmap = np.ascontiguousarray(cmap.real)
        rmap = coordinateutils.FlatSkyMap(
            rmap,
            self.dx,
            is_weighted=False,
            units=self.units,
            pol_type=getattr(coordinateutils.MapPolType, pol_type),
            proj=self.proj
        )
        return {pol_type: rmap}

    def get_l_masked(
        self,
        lmin=None,
        lmax=None,
        lxmin=None,
        lxmax=None,
        lymin=None,
        lymax=None,
    ):
        """ returns a copy of this object which has been
            masked to zero in a customizable range of Fourier space.
        """
        ret = self.copy()
        lx, ly = ret.get_lxly()
        ell = np.hypot(lx, ly)
        if lmin is not None:
            ret[ell < lmin] = 0.0
        if lmax is not None:
            ret[ell >= lmax] = 0.0
        if lxmin is not None:
            ret[np.abs(lx) < lxmin] = 0.0
        if lymin is not None:
            ret[np.abs(ly) < lymin] = 0.0
        if lxmax is not None:
            ret[np.abs(lx) >= lxmax] = 0.0
        if lymax is not None:
            ret[np.abs(ly) >= lymax] = 0.0
        return ret

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
        ret = self.copy()
        ret[:, :] = 1.0
        lx, ly = ret.get_lxly()
        ell = np.hypot(lx, ly)
        if lmin is not None:
            ret[ell < lmin] = 0.0
        if lmax is not None:
            ret[ell >= lmax] = 0.0
        if lxmin is not None:
            ret[np.abs(lx) < lxmin] = 0.0
        if lymin is not None:
            ret[np.abs(ly) < lymin] = 0.0
        if lxmax is not None:
            ret[np.abs(lx) >= lxmax] = 0.0
        if lymax is not None:
            ret[np.abs(ly) >= lymax] = 0.0
        return ret

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(MapSpectrum2D, self).__reduce__()
        outdict = self.__dict__.copy()
        # Convert these unpickleable enum things into strings
        outdict['units'] = str(outdict['units'])
        outdict['proj'] = str(outdict['proj'])
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (outdict,)
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        # Fill out __dict__
        self.__dict__.update(state[-1])
        self.units = core.G3TimestreamUnits.names[self.units]
        self.proj = coordinateutils.MapProjection.names[self.proj]
        # Call the parent's __setstate__ with the other tuple elements.
        super(MapSpectrum2D, self).__setstate__(state[0:-1])


class MapSpectrum1D(np.ndarray):
    """Class to hold 1-dimensional power spectrum.

    Selected attributes:
    -----------
    spec_type: string or None
        Options are: 'cl', 'dl', or None (unspecified).
    lbins: list or 1d numpy array
        the bin edges.
    bin_centers: 1d numpy array
        the bin centers.
    delta_l: int or float
        the size of binning, default is one.
    lmax: int or float
        the maximum ell to keep track of.
    res: float
        the resolution of the map.
    map_nx: int
        map pixel number in x
    map_ny: int
        map pixel number in y
    num_modes: 1d numpy array
        number of modes (number or pixels in each annulus).
    """

    def __new__(cls, lbins, spec, spec_type=None, dx=None,
                dy=None, map_nx=None, map_ny=None,
                units=getattr(core.G3TimestreamUnits, "None"),
                proj=getattr(coordinateutils.MapProjection, "ProjNone"),
                is_asd=False):
        """
        Selected arguments:
        ----------
        lbins: list, 1d numpy array, or list of tuples
            the bin edges.
        spec: list or numpy array
            the power spectrum value (cl or dl).
        spec_type: string or None
            Options are: 'cl', 'dl', or None (unspecified).
        dx, dy: float
            the x/y resolution of the map.
        map_nx: int
            map pixel number in x
        map_ny: int
            map pixel number in y
        """

        # initialize
        obj = np.asarray(spec).view(cls)
        obj.spec_type = spec_type
        # Convert lbins to a 1d list if it's a list of tuple. For example,
        # if it's [(0,100),(100,200)], convert it to [0, 100, 200]
        if isinstance(lbins[0], (list, tuple)):
            temp = [seq[0] for seq in lbins]
            temp.append(lbins[-1][1])
            lbins = np.array(temp)
        lbins = np.array(lbins)
        obj.lbins = lbins
        obj.bin_centers = (obj.lbins[1:] + obj.lbins[:-1]) / 2.0
        # if uniform binning
        if np.allclose(
            lbins[1:] - lbins[:-1], lbins[1] - lbins[0],
            rtol=1.0e-8, atol=1.0e-8
        ):
            obj.delta_l = lbins[1] - lbins[0]
        else:
            obj.delta_l = None
        obj.lmax = np.max(lbins)
        if map_ny is None:
            map_ny = map_nx
        if dy is None:
            dy = dx
        obj.map_nx, obj.map_ny, obj.dx, obj.dy = map_nx, map_ny, dx, dy
        obj.units, obj.proj, obj.is_asd = units, proj, is_asd
        return obj

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        outputs = kwargs.get('out', ())
        # check that data type is compatible
        for x in inputs + outputs:
            if isinstance(x, self.__class__):
                if not self.compatible(x):
                    return NotImplemented
            elif np.isscalar(x):
                pass
            elif isinstance(x, (np.ndarray, list)) and np.shape(x) == self.shape:
                pass
            elif isinstance(x, (np.ndarray, list)) and np.shape(x) != self.shape:
                raise ValueError("operands could not be broadcast together with shapes"
                                 + str(np.shape(x)) + str(self.shape))
            else:
                return NotImplemented
        # gather input arguments for ufunc
        args = []
        for i, input in enumerate(inputs):
            if isinstance(input, self.__class__):
                args.append(input.view(np.ndarray))
            else:
                args.append(input)
        if outputs:
            kwargs['out'] = tuple(
                x.view(np.ndarray) if isinstance(x, self.__class__) else x
                for x in outputs)
        # do the calculation
        results = super().__array_ufunc__(ufunc, method,
                                          *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented
        if method == 'at':
            return None
        if ufunc.nout == 1:
            results = (results,)
        results = (np.nan_to_num(result) for result in results)
        # return the object with the same attributes 
        results = [result.view(self.__class__) for result in results]
        for result in results:
            if isinstance(result, self.__class__):
                for attr in ["lbins", "spec_type", "dx", "dy",
                             "map_nx", "map_ny", "units", "proj",
                             "is_asd", "delta_l", "bin_centers",
                             "lmax",  "num_modes"]:
                    setattr(result, attr, getattr(self, attr, None))
        results=tuple(results)
        return results[0] if len(results) == 1 else results

    def __array_finalize__(self, obj):
        if obj is None: return
        # same attributes when create new from template
        for attr in ["lbins", "spec_type", "dx", "dy",
                     "map_nx", "map_ny", "units", "proj", 
                     "is_asd", "delta_l", "bin_centers", 
                     "lmax",  "num_modes"]:
            if not hasattr(self, attr):
                setattr(self, attr, getattr(obj, attr, None))


    def compatible(self, other):
        """ Check if this object can be added, subtracted, etc. with other.
        """

        return (
            isinstance(other, MapSpectrum1D)
            and (self.delta_l == other.delta_l or 
                 np.allclose(self.delta_l, other.delta_l, rtol=1.0e-8, atol=1.0e-8))
            and self.spec_type == other.spec_type
            and self.map_nx == other.map_nx
            and self.map_ny == other.map_ny
            and np.allclose(self.dx, other.dx, rtol=1.0e-8, atol=1.0e-8)
            and np.allclose(self.dy, other.dy, rtol=1.0e-8, atol=1.0e-8)
        )

    def get_num_modes(self):
        if hasattr(self, 'num_modes') and self.num_modes is not None:
            return self.num_modes
        if (self.dx is not None) and (self.map_nx is not None) and (
            self.map_ny is not None
        ):
            from .basicmaputils import make_ellgrid
            ell = make_ellgrid(self.dx, (self.map_ny, self.map_nx))
            self.num_modes, temp = np.histogram(ell.ravel(), bins=self.lbins)
        else:
            self.num_modes = None
        return self.num_modes

    def rebin(self, lbins, w=lambda l: 1.0):
        """
        Rebins this spectrum with non-uniform binning

        Arguments
        ----------
        lbins: 1D array or an array of tuples
            The bin edges.
        w: a function
            weight function to apply when accumulating into bins.
        """
        num_modes = self.get_num_modes()
        if num_modes is None:
            raise ValueError("Cannot rebin if number of modes is None!")
        if self.delta_l != 1:
            warnings.warn("Default Delta l for rebinning is one!")
        # Convert lbins to a 1d list if it's a list of tuple. For example,
        # if it's [(0,100),(100,200)], convert it to [0, 100, 200]
        if isinstance(lbins[0], (list, tuple)):
            temp = [seq[0] for seq in lbins]
            temp.append(lbins[-1][1])
            lbins = np.array(temp)
        lbins = np.array(lbins)
        # lb is the new bin centers
        lb = 0.5 * (lbins[:-1] + lbins[1:])
        wb = w(lb)
        # l is the original bin centers
        l = self.lbins
        l = 0.5 * (l[:-1] + l[1:])
        w = w(l)

        # get the weights in each l bin, proportional to # modes
        norm, temp = np.histogram(
            l, bins=lbins, weights=np.nan_to_num(num_modes)
        )
        if self.is_asd:
            unbinned_psd = np.nan_to_num(np.array(self))**2
        else:
            unbinned_psd = np.nan_to_num(np.array(self))
        # rebin the power spectrum
        binned_psd, temp = np.histogram(
            l,
            bins=lbins,
            weights=w
            * np.nan_to_num(num_modes)
            * unbinned_psd,
        )
        # divide by the weight
        binned_psd[np.nonzero(norm)] /= norm[np.nonzero(norm)] * wb
        if self.is_asd:
            binned_spec = np.sqrt(binned_psd)
        else:
            binned_spec = binned_psd
        binned = self.__class__(
            lbins,
            binned_spec,
            spec_type=self.spec_type,
            dx=self.dx,
            dy=self.dy,
            map_nx=self.map_nx,
            map_ny=self.map_ny,
            units=self.units,
            proj=self.proj,
            is_asd=self.is_asd,
        )
        return binned

    def get_dl(self):
        """ Convert Cl into Dl, or return self if already in Dl

        Not very accurate. Better to calculate from 2D!
        """
        if self.is_asd:
            raise NotImplementedError('Convert to psd first.')
        if self.spec_type == "dl":
            return self
        elif self.spec_type == "cl":
            dl = (
                np.array(self)
                * self.bin_centers
                * (self.bin_centers + 1)
                / 2
                / np.pi
            )
            return self.__class__(
                self.lbins,
                dl,
                spec_type="dl",
                dx=self.dx,
                dy=self.dy,
                map_nx=self.map_nx,
                map_ny=self.map_ny,
                units=self.units,
                proj=self.proj,
                is_asd=self.is_asd,
            )
        else:
            raise ValueError("spec_type not specified!")

    def get_cl(self):
        """ Convert Dl into Cl, or return self if already in Cl

        Not very accurate. Better to calculate from 2D!
        """
        if self.is_asd:
            raise NotImplementedError('Convert to psd first.')
        if self.spec_type == "cl":
            return self
        elif self.spec_type == "dl":
            cl = (
                np.array(self)
                / self.bin_centers
                / (self.bin_centers + 1)
                * 2
                * np.pi
            )
            return self.__class__(
                self.lbins,
                cl,
                spec_type="cl",
                dx=self.dx,
                dy=self.dy,
                map_nx=self.map_nx,
                map_ny=self.map_ny,
                units=self.units,
                proj=self.proj,
                is_asd=self.is_asd,
            )
        else:
            raise ValueError("spec_type not specified!")

    def get_2d(self, spec_2d = None):
        """Interpolate to 2D grid and return MapSpectrum2D

        If a MapSpectrum2D object is given as an argument, it will
        interpolate to the ell grid corresponding to it. Otherwise
        it will interpolate to the ell grid corresponding to self's
        shape and resolution
        """
        #if no arguments given
        if spec_2d is None and (self.map_nx is not None) and (
            self.map_ny is not None
        ) and (self.dx is not None):
            spec_2d = MapSpectrum2D(
                self.map_nx, self.dx, spec=None,
                map_ny=self.map_ny, dy=self.dy,
                units=self.units, proj=self.proj)
        #if a spec_2d object is given, interpolate to its ell grid
        elif isinstance(spec_2d, MapSpectrum2D):
            spec_2d = spec_2d.copy()
        else:
            return NotImplemented
        ell = spec_2d.get_ell()
        spec_2d[:,:] = np.interp(
            ell.ravel(), self.bin_centers, np.array(self), right=0
        ).reshape(spec_2d.shape)
        return spec_2d

    def copy(self):
        """ return a clone of this object. """
        return MapSpectrum1D(
            self.lbins,
            self.view(np.ndarray).copy(),
            spec_type=self.spec_type,
            dx=self.dx, dy=self.dy,
            map_nx=self.map_nx,
            map_ny=self.map_ny,
            units=self.units,
            proj=self.proj,
            is_asd=self.is_asd)

    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(MapSpectrum1D, self).__reduce__()
        outdict = self.__dict__.copy()
        # Convert these unpickleable enum things into strings
        outdict['units'] = str(outdict['units'])
        outdict['proj'] = str(outdict['proj'])
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (outdict,)
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        # Fill out __dict__
        self.__dict__.update(state[-1])
        self.units = core.G3TimestreamUnits.names[self.units]
        self.proj = coordinateutils.MapProjection.names[self.proj]
        # Call the parent's __setstate__ with the other tuple elements.
        super(MapSpectrum1D, self).__setstate__(state[0:-1])


class MapSpectrum1DDict(dict):
    """ A collection of MapSpectrum1D objects.

    This class iterates an operation to all objects it holds.

    Attributes
    -----------
    specs: a dictionary containing MapSpectrum1D objects
    """

    def __init__(self, ps):
        """ Initialize with a dictionary containing MapSpectrum1D objects

        Arguments
        ---------
        ps: a dictionary
            A dictionary in the format of {'TT': MapSpectrum1D object, ...}
        """
        for k, v in ps.items():
            if not isinstance(v, MapSpectrum1D):
                raise TypeError("Argument {} is not a MapSpectrum1D".format(k))
        super().__init__(ps)

    def copy(self):
        return MapSpectrum1DDict({k:self[k].copy() for k in self.keys()})

    def __add__(self, other):
        ret = self.copy()
        ret += other
        return ret

    def __iadd__(self, other):
        if isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] += other[spec]
            return self
        elif isinstance(other, MapSpectrum1D) or np.isscalar(other):
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
        elif isinstance(other, MapSpectrum1D) or np.isscalar(other):
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
        if isinstance(other, self.__class__) and other.keys() == self.keys():
            for spec in self.keys():
                self[spec] *= other[spec]
            return self
        elif isinstance(other, MapSpectrum1D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] *= other
            return self
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
        elif isinstance(other, MapSpectrum1D) or np.isscalar(other):
            for spec in self.keys():
                self[spec] /= other
            return self
        else:
            return NotImplemented

    def rebin(self, lbins, w=lambda l: 1.0):
        rebinned = {}
        for spec in self.keys():
            rebinned[spec] = self[spec].rebin(lbins, w)
        return self.__class__(rebinned)

    def get_dl(self):
        dl = {}
        for spec in self.keys():
            dl[spec] = self[spec].get_dl()
        return self.__class__(dl)

    def get_cl(self):
        cl = {}
        for spec in self.keys():
            cl[spec] = self[spec].get_cl()
        return self.__class__(cl)
