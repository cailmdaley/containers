#!/usr/bin/env python
import numpy as np
from spt3g import core, coordinateutils
from spt3g.mapspectra import basicmaputils, map_analysis
from spt3g.mapspectra.map_spectrum_classes import MapSpectrum1DDict

# generate 1uK arcmin noise TQU maps
dx = 40 * core.G3Units.arcmin
dy = dx
tmap = (
    np.random.standard_normal([75, 113])
    * 1
    * core.G3Units.uK
    * core.G3Units.arcmin
    / (np.sqrt(dx * dy))
)
tmap = coordinateutils.FlatSkyMap(
    tmap,
    dx,
    proj=coordinateutils.MapProjection.Proj5,
    alpha_center=0.0,
    delta_center=-57.5 * core.G3Units.deg,
    coord_ref=coordinateutils.MapCoordReference.Equatorial,
    is_weighted=False,
)
maps = {"T": tmap, "Q": tmap, "U": tmap}

# generate a power spectrum with delta_l=1
cls = map_analysis.calculate_powerspectra(maps, flatten=True, apod_mask=None)
cls = MapSpectrum1DDict(cls)
ell_bins = np.arange(25, 426, 200)

# generate a power spectrum with customized ell_bins
cls_binned = map_analysis.calculate_powerspectra(
    maps, lbins=ell_bins, flatten=True, apod_mask=None
)

# rebin from delta_l = 1
cls_rebinned = cls.rebin(ell_bins)

# rebin should generate the same results as defining ell bins earlier
assert (
    np.max(abs(cls_binned["TT"] - cls_rebinned["TT"]) / abs(cls_binned["TT"]))
    < 1e-12
), "Rebinning fails."

# the map is 1uK arcmin noise. check the power spectrum is 1uK arcmin
assert all(
    abs(cls_binned["TT"] / (core.G3Units.arcmin * core.G3Units.uK) ** 2 - 1)
    < 0.1
), "Power spectrum values are incorrect."

# check division, multiplication, and subtraction operations
operation_test = cls_binned["TT"] / 2 * 2 - cls_binned["TT"]
assert all(operation_test < 1e-12), "Operation fails for MapSpectrum1D."

# check that operation preserve attributes
assert (
    cls_binned["TT"].__dict__.keys() == operation_test.__dict__.keys()
    and operation_test.units == cls_binned["TT"].units
), "Operation fails for MapSpectrum1D."

# check multiplication and subtraction with objects of the same type
cls_square = cls * cls
nonzero = np.where(cls["TT"] != 0)
assert abs(
    np.nanmean(
        np.array(cls_square["TT"][nonzero] / cls["TT"][nonzero] ** 2) - 1
    )
    < 1e-12
), "Operation fails for MapSpectrum1D"

# Test get_dl()
dl = cls.get_dl()
ells = np.arange(50, 1000 + 1, 50)
dl_rebinned = dl.rebin(ell_bins)
cls_from_dl = dl_rebinned.get_cl()
assert dl["TT"].spec_type == "dl", "dl has incorrect spec_type"
assert cls_from_dl["TT"].spec_type == "cl", "cl has incorrect spec_type"

# the more accurate way to calculate dl is in 2d space instead of converting from cl
dl_accurate = map_analysis.calculate_powerspectra(
    maps, lbins=ell_bins, flatten=True, apod_mask=None, calculate_dls=True
)

# compare the dl calculated from 2D and the dl converted from 1D cl
diff = dl_accurate["TT"] - dl_rebinned["TT"]
assert all(
    (diff) < 1e-8
), "dl calculated from 2D is different from the one calculated from 1D cl"

# QU power spectrum
cls_queb = map_analysis.calculate_powerspectra(
    maps,
    flatten=True,
    b_mode_method="chi",
    apod_mask=None,
    symmetric_cross_spectra=True,
    qu_eb="both",
)
cls_binned_queb = MapSpectrum1DDict(cls_queb).rebin(ell_bins)

# Q and U should be 1 uK arcmin noise
assert all(
    abs(
        cls_binned_queb["QQ"] / (core.G3Units.arcmin * core.G3Units.uK) ** 2 - 1
    )
    < 0.1
), "Power spectrum values are incorrect."
assert all(
    abs(
        cls_binned_queb["UU"] / (core.G3Units.arcmin * core.G3Units.uK) ** 2 - 1
    )
    < 0.1
), "Power spectrum values are incorrect."

# test 2D power spectrum
cls_2d = map_analysis.calculate_powerspectra(
    maps,
    flatten=True,
    b_mode_method="chi",
    apod_mask=None,
    symmetric_cross_spectra=True,
    return_2d=True,
)
bb = cls_2d["BB"]

# test operation
assert (
    np.max(abs(bb * bb - bb ** 2)) < 1e-12
), "Operation fails for MapSpectrum2D"

# test the masking
masked = bb.get_l_mask(lmax=400)

# masked area should be zero
assert all(masked[masked.get_ell() > 400] < 1e-12), "Ell pasking does not work"

# test get_real
mapspectra = basicmaputils.map_to_ft(maps["T"])
real = mapspectra.get_real()

# from fft to map
tmap = real.get_rmap()

# the recovered map from real fft should be the same as the input map
assert np.max(
    np.abs(tmap["T"] - maps["T"]) < 1e-12
), "get_rmap or get_real does not work."
