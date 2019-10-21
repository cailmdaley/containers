#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import argparse as ap

# Usage: fit_fluxandpointing.py <files.g3> -o output.g3
#
# Computes best-fit relative detector pointing and flux calibration
# for all detectors in the focal plane. Flux results are normalized to
# a 4 arcminute x 4 arcminute box centered on the brightest point
# in the input maps.

P = ap.ArgumentParser(description='Pointing and calibration off of a point source',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-s', '--source', action='store', default='RCW38',
               help='Source name for output data')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-v', '--verbose', action='store_true', default=False,
               help='Print internal likelihoods')
args = P.parse_args()

# Load source maps from disk
map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
bp = None
calresponse = None
calresponsesn = None
def extract_bp(fr):
    global bp, calresponse, calresponsesn
    if 'BolometerProperties' in fr and bp is None:
        bp = fr['BolometerProperties']
    if 'CalibratorResponse' in fr:
        calresponse = fr['CalibratorResponse']
        calresponsesn = fr['CalibratorResponseSN']

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(extract_bp)
# Cut maps that are all zero
pipe.Add(lambda fr: 'Id' not in fr or 'Wunpol' in fr or ((numpy.asarray(fr['T']) != 0).any() and numpy.isfinite(fr['T']).any()))
pipe.Add(map_extractor)
pipe.Run()

data = map_extractor.maps
assert(calresponse is not None)
assert(calresponsesn is not None)

# We want 5% or better flux calibration, so ignore any detector with
# calibrator S/N < 20. Since the results of this code are calibrator
# relative, we will never do better than 5% in such circumstances no
# matter how good the fit is.
allbolos = list(filter(lambda b: b in calresponsesn and calresponsesn[b] > 20, data.keys()))

assert(len(allbolos) > 0)

# Build out our fit groups
if bp is not None:
    # Default is to return lists of bolos in each band, keyed by observing band
    groups = calibration.template_groups.get_template_groups(bp, include_keys=True, exclude_empty=True)
    bolodict = {}
    for group in groups.keys():
        groupbolos = [b for b in groups[group] if b in allbolos]
        if len(groupbolos) > 0:
            bolodict[group] = groupbolos
    allbolos = bolodict
else:
    allbolos = {'All': allbolos}
print('Fitting %d bolometer groups' % len(allbolos))

# Start fit

# NB: numpy and the mapmaker use opposite definitions of "X" and "Y". "X", unless
# noted otherwise, here follows the numpy convention (first, slowly varying axis)
# rather than the sane convention used elsewhere.

params = {} # bolo: (amp, delta X, delta Y)

# Get global map parameters
w = None
# Only the one weight map since all the hit maps are the same
for b in data.keys():
    if 'Wunpol' in data[b] and b != 'bsmap' and b != 'map':
        w = numpy.asarray(data[b]['Wunpol'].TT)[:]
        break
assert(w is not None)
mapres = list(data.values())[0]['T'].res

# Apodize the weight map a little in the scan direction
for row in range(len(w)):
    # Find where weight map transitions to or from 0
    trans = numpy.where(numpy.diff((w[row] == 0).astype('float')))[0]
    if len(trans) >= 2:
        w[row][:trans[0] + 5] = 0
        w[row][trans[-1] - 5:] = 0
    elif len(trans) == 1:
        if trans[0] > len(w[row])/2:
            w[row][trans[0] - 5:] = 0
        else:
            w[row][:trans[0] + 5] = 0

# Utility function for building a template by coadding
# shifted maps.
def recalctemplate(maps, weights, params, maxmaps=0):
    template = numpy.zeros(list(maps.values())[0].shape)
    template_weights = numpy.zeros(list(maps.values())[0].shape)
    x,y = numpy.indices(template.shape)
    i = 0
    ampls = numpy.asarray([params[bolo][0] for bolo in maps.keys()])
    avgampls = numpy.nanmedian(ampls)

    for bolo in maps.keys():
        # Drop wonky maps (zero or NaN amplitude, amplitude very
        # different from average) from the coadd
        if params[bolo][0] == 0 or not numpy.isfinite(params[bolo][0]):
            continue
        if numpy.abs((params[bolo][0] - avgampls)/avgampls) > .2:
            continue

        # Also avoid non-responsive bolometers
        # XXX: Should already be excluded by allbolos construction
        if calresponsesn[bolo] < 30:
            continue

        # Shift the bolometer map by the current best-fit pointing
        # offset and construct coadd of normalized maps
        template += scipy.ndimage.map_coordinates(maps[bolo]/params[bolo][0], [x + params[bolo][1], y + params[bolo][2]], order=1)
        template_weights += scipy.ndimage.map_coordinates(weights, [x + params[bolo][1], y + params[bolo][2]], order=1)

        # We may not need a huge number of maps to get high S/N in
        # the template (just need errors << the bolo map errors)
        i += 1
        if maxmaps != 0 and i > maxmaps:
            break

    # Zero everything outside of centered 30x30 pixel box 
    template[:template.shape[0]//2 - 30,:] = 0
    template[template.shape[0]//2 + 30:,:] = 0
    template[:,:template.shape[1]//2 - 30] = 0
    template[:,template.shape[1]//2 + 30:] = 0

    # Check the brightest pixel is actually in the middle.
    # So long as the algorithm isn't diverging, it has to
    # stay there.
    maxpix = numpy.unravel_index(numpy.nanargmax(template/template_weights), template.shape)
    assert(abs(maxpix[0] - template.shape[0]//2) < 2)
    assert(abs(maxpix[1] - template.shape[1]//2) < 2)

    return template, template_weights

# Iterate over all subgroups of detectors, with independent dynamic templates
final_templates = {}
renormalization = {}
for group, bolos in allbolos.items():
    rawmaps = {b: numpy.asarray(data[b]['T'])[:] for b in bolos}

    # First, pick the brightest (= highest S/N) pixel and use that as
    # a seed in both amplitude and position
    for bolo in rawmaps.keys():
        m = rawmaps[bolo]

        sn = m / w**0.5
        sn[numpy.logical_not(numpy.isfinite(sn))] = 0
        
        # Smooth to avoid hot pixels
        sn = scipy.ndimage.gaussian_filter(sn, sigma=1*core.G3Units.arcmin/mapres)

        max = numpy.unravel_index(numpy.argmax(numpy.abs(sn)), m.shape)
        amp = (m/w)[max]
        # Get offset of brightest pixel from map center
        max = (max[0] - m.shape[0]/2, max[1] - m.shape[1]/2)
        params[bolo] = (amp, max[0], max[1])
        if args.verbose:
            print('%s: %s' % (bolo, params[bolo]))

    # Second, iteratively fit the maps of each bolometer to a shifted
    # and scaled copy of the current co-add. This is an approximation
    # to the full MLH solution in which the likelihood is fully
    # seperable (i.e. we assume the covariance between the co-add template 
    # and the parameters for the bolometers is zero).

    deltachisq = {} # Dictionary of change in (approximate) chi-squared between subsequent iterations

    for iter in range(0,3):
        print('Starting iteration %d of group %s' % (iter, group))
        template, template_weights = recalctemplate(rawmaps, w, params, maxmaps=60)
        template = template/template_weights
        template[numpy.logical_not(numpy.isfinite(template))] = 0

        # Normalize center peak to 1
        # XXX: template is already peaked at 1?
        template /= template[template.shape[0]//2, template.shape[1]//2]

        # On the first iteration, renormalize the first guess, based on
        # the brightest pixel, to this convention
        if iter == 0:
            renorm = numpy.max(template)
            for bolo in rawmaps.keys():
                params[bolo] = (params[bolo][0]/renorm, params[bolo][1], params[bolo][2])

        template_spl = scipy.interpolate.RectBivariateSpline(
          numpy.arange(template.shape[0]),numpy.arange(template.shape[1]),template)
        for bolo in rawmaps.keys():
            unwmap = rawmaps[bolo]/w
            unwmap[numpy.logical_not(numpy.isfinite(unwmap))] = 0
            invnoise = w**0.5
            invnoise /= numpy.std(unwmap[unwmap != 0]*invnoise[unwmap != 0])
            errorfunction = lambda p: numpy.ravel(
            (p[0]*template_spl(
                    numpy.arange(template.shape[0]) - p[1],
                    numpy.arange(template.shape[1]) - p[2]) 
             - unwmap)*invnoise)
            p, success = scipy.optimize.leastsq(errorfunction, params[bolo])
            deltachisq[bolo] = numpy.sum(errorfunction(params[bolo])**2) - \
            numpy.sum(errorfunction(p)**2)
            if args.verbose:
                print('\tMoved %s from %s to %s: %f' % (bolo, params[bolo], p, deltachisq[bolo]))
            params[bolo] = p

    # Compute integral within 4x4 arcminute box for normalization purposes
    box = 2*core.G3Units.arcmin # +/- 2 arcmin
    integrand = template*(mapres**2)
    norm = numpy.sum(
      integrand[template.shape[0]//2 - int(box/mapres):template.shape[0]//2 + int(box/mapres),
      template.shape[1]//2 - int(box/mapres):template.shape[1]//2 + int(box/mapres)]
    )

    renormalization[group] = norm
    final_templates[group] = template

# Write output

outframe = core.G3Frame(core.G3FrameType.Calibration)

fluxcal = core.G3MapDouble()
peaktointegral = core.G3MapDouble()
xoffset = core.G3MapDouble()
yoffset = core.G3MapDouble()
# XXX: also error bars on these quantities?

# Get Az/el offsets by using lookup maps to transform from (fractional) pixel
# indices to the offsets from the source at that moment
for bolo, p in params.items():
    fluxcal[bolo] = p[0]/calresponse[bolo]
    xoffset[bolo] = p[2]*mapres
    yoffset[bolo] = p[1]*mapres
        

for group in allbolos:
    for bolo in allbolos[group]:
        peaktointegral[bolo] = renormalization[group]
# Peak brightness of source, relative to calibrator response
outframe[args.source + 'FluxCalibration'] = fluxcal
# Flux normalized to peak at 1, then integrated over 4x4 arcminute box
outframe[args.source + 'IntegralFlux'] = peaktointegral

assert(len(fluxcal) > 0)
assert(len(peaktointegral) > 0)

outframe['PointingOffsetX'] = xoffset
outframe['PointingOffsetY'] = yoffset

print(outframe)

writer = core.G3Writer(args.output)
writer(outframe)
writer(core.G3Frame(core.G3FrameType.EndProcessing)) # Make sure data gets flushed

