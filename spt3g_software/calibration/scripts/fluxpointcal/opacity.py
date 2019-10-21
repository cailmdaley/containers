#!/usr/bin/env python

from __future__ import print_function

import numpy, sys, re
import scipy.ndimage
from spt3g import core, mapmaker
from spt3g.std_processing.sourcemaps import get_source_ra_dec
import argparse as ap

# Usage: opacity.py -i <files.g3> -o output.g3

P = ap.ArgumentParser(description='Sky opacity from a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-s', '--source', action='store', default='RCW38',
               help='Name of source for output data')
args = P.parse_args()

# Load maps from disk
map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
integrals = None
fluxcals = None
bp = None
def getintegralfluxes(fr):
    global integrals
    global fluxcals
    global bp
    if args.source + 'IntegralFlux' in fr:
        integrals = fr[args.source + 'IntegralFlux']
    if args.source + 'FluxCalibration' in fr:
        fluxcals = fr[args.source + 'FluxCalibration']
    if 'BolometerProperties' in fr:
        bp = fr['BolometerProperties']
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(getintegralfluxes)
pipe.Add(map_extractor)
pipe.Run()

data = map_extractor.maps

# Write output

outframe = core.G3Frame(core.G3FrameType.Calibration)

trans = core.G3MapDouble()
x_offset = core.G3MapDouble()
y_offset = core.G3MapDouble()
dpmodel = core.G3MapVectorDouble()
for m in data.keys():
    b = re.match(args.source + '-([0-9]+)GHz', m)
    if b is not None:
        band = int(b.groups()[0])

        if numpy.max(data[m]['Wunpol'].TT) == 0:
            core.log_warn(str(band)+" GHz map for "+args.output+" is all zeros.")
            continue

        # Get an appropriate detector for the integral
        fluxintegral = None
        for bolo, props in bp.iteritems():
            if props.band == band*core.G3Units.GHz and bolo in integrals:
                fluxintegral = integrals[bolo]
                break

        # Compute integral within 4x4 arcminute box for normalization purposes.
        # This is hardcoded because it *must* match the values in fit_fluxandpointing.py
        # and the values in apply_t_calibration.py.
        box = 2*core.G3Units.arcmin # +/- 2 arcmin
        integrand = numpy.asarray(data[m]['T'])*(data[m]['T'].res**2)
        mapres = data[m]['T'].res

        # Find the brightest pixel
        with numpy.errstate(divide='ignore', invalid='ignore'):
            sn = numpy.asarray(data[m]['T']) / data[m]['Wunpol'].TT**0.5
        sn[numpy.logical_not(numpy.isfinite(sn))] = 0
        # Smooth to avoid hot pixels
        sn = scipy.ndimage.gaussian_filter(sn, sigma=1*core.G3Units.arcmin/mapres)
        max = numpy.unravel_index(numpy.argmax(numpy.abs(sn)), sn.shape)

        # Compute the weighted integral of the source in a box around the peak pixel
        mask = (
          slice(max[0] - int(box/mapres), max[0] + int(box/mapres)),
          slice(max[1] - int(box/mapres), max[1] + int(box/mapres))
        )
        norm = numpy.sum(integrand[mask])/numpy.sum(numpy.asarray(data[m]['Wunpol'].TT)[mask])*numpy.sum(numpy.ones(integrand[mask].shape))

        trans[str(band)] = norm/fluxintegral

        # Compute the flux-weighted centroid
        with numpy.errstate(divide='ignore', invalid='ignore'):
            map_unw = numpy.asarray(data[m]['T']) / data[m]['Wunpol'].TT
        print('%d GHz: %f (peak %f)' % (band, trans[str(band)], numpy.nanmax(map_unw)))
        xybox = numpy.indices(sn.shape)
        ycen = numpy.sum(map_unw[mask]*(xybox[0,:,:])[mask])/numpy.sum(map_unw[mask])
        y_offset[str(band)] = (ycen - (numpy.float(map_unw.shape[0]//2)-0.5))*mapres
        xcen = numpy.sum(map_unw[mask]*(xybox[1,:,:])[mask])/numpy.sum(map_unw[mask])
        x_offset[str(band)] = (xcen - (numpy.float(map_unw.shape[1]//2)-0.5))*mapres

assert(len(trans) > 0)
outframe[args.source + 'SkyTransmission'] = trans
outframe[args.source + 'OffsetX'] = x_offset
outframe[args.source + 'OffsetY'] = y_offset
outframe[args.source + 'IntegralFlux'] = integrals
outframe[args.source + 'FluxCalibration'] = fluxcals

# translate weighted offsets into pointing parameters using
# first-order approximation to the pointing model
if args.source == 'RCW38' or args.source == 'MAT5A':
    el_source = -1 * get_source_ra_dec(args.source)[1]
    xoffs = numpy.asarray(x_offset.values())
    yoffs = numpy.asarray(y_offset.values())
    if len(numpy.where(numpy.isfinite(xoffs))[0]) > 0 and len(numpy.where(numpy.isfinite(yoffs))[0]) > 0:
        xoff_mean = numpy.nanmean(x_offset.values())
        yoff_mean = numpy.nanmean(y_offset.values())
        dpmodel['fixedCollimation'] = numpy.asarray([0.,-yoff_mean])
        dpmodel['flexure'] = numpy.asarray([0.,0.])
        dpmodel['tilts'] = numpy.asarray([0.,0.,xoff_mean/numpy.sin(el_source)])
    outframe[args.source + 'PointingModelCorrection'] = dpmodel

writer = core.G3Writer(args.output)
writer(outframe)
writer(core.G3Frame(core.G3FrameType.EndProcessing)) # Make sure data gets flushed

