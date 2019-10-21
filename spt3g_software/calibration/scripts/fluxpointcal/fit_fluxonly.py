from __future__ import print_function

import numpy, sys
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import argparse as ap

# Usage: fit_fluxonly.py -i <files.g3> -o output.g3

P = ap.ArgumentParser(description='Flux calibration from a point source',
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
bp = None
calresponse = None
def extract_bp(fr):
    global bp, calresponse
    if 'BolometerProperties' in fr:
        bp = fr['BolometerProperties']
    if 'CalibratorResponse' in fr:
        calresponse = fr['CalibratorResponse']

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(extract_bp)
pipe.Add(map_extractor)
pipe.Run()

data = map_extractor.maps
assert(calresponse is not None)

# Build out our fit groups
if bp is not None:
    groups = calibration.template_groups.get_template_groups(bp)
    allbolos = [[b for b in group if b in data.keys()] for group in groups]
else:
    allbolos = [data.keys()]
print('Fitting %d bolometer groups' % len(allbolos))

# Start fit

# NB: numpy and the mapmaker use opposite definitions of "X" and "Y". "X", unless
# noted otherwise, here follows the numpy convention (first, slowly varying axis)
# rather than the sane convention used elsewhere.

fluxcal = core.G3MapDouble()

# Utility function for building a template by coadding maps.
def recalctemplate(maps, weights, params):
    template = numpy.zeros(maps.values()[0].shape)
    template_weights = numpy.zeros(maps.values()[0].shape)
    x,y = numpy.indices(template.shape)
    ampls = numpy.asarray([params[a] for a in maps.keys()])
    avgampls = numpy.nanmedian(ampls)

    for bolo in params.keys():
        # Drop wonky maps (zero or NaN amplitude, amplitude very
        # different from average) from the coadd
        if params[bolo] == 0 or not numpy.isfinite(params[bolo]):
            continue
        if numpy.abs((params[bolo] - avgampls)/avgampls) > .2:
            continue

        # Add to coadd
        template += maps[bolo]/params[bolo]
        template_weights += weights[bolo]

    return template, template_weights

def mapextremum(m):
    max = numpy.nanmax(m)
    min = numpy.nanmin(m)
    if numpy.abs(max) > numpy.abs(min):
        return max
    else:
        return min

# Iterate over all subgroups of detectors, with independent dynamic templates
for bolos in allbolos:
    rawmaps = {b: data[b]['T'][:] for b in bolos}
    weights = {b: data[b]['Wunpol'].TT[:] for b in bolos}

    # First guess to assist in template building
    for bolo in bolos:
        fluxcal[bolo] = mapextremum(rawmaps[bolo]/weights[bolo])

    template, template_weights = recalctemplate(rawmaps, weights, fluxcal)
    template = template/template_weights
    template[numpy.logical_not(numpy.isfinite(template))] = 0

    # Normalize brightest point to 1
    template /= mapextremum(template)

    # Solve analytically for the amplitude of each bolometer
    for bolo in bolos:
        w = weights[bolo]
        m = rawmaps[bolo]
        m[numpy.logical_not(numpy.isfinite(m))] = 0

        fluxcal[bolo] = numpy.sum(m*template)/numpy.sum(template**2 * w) / calresponse[bolo]
        print(bolo, fluxcal[bolo])

# Write output

outframe = core.G3Frame(core.G3FrameType.Calibration)

# XXX: also error bars on these quantities?

outframe[args.source + 'FluxCalibration'] = fluxcal

writer = core.G3Writer(args.output)
writer(outframe)
writer(core.G3Frame(core.G3FrameType.EndProcessing)) # Make sure data gets flushed

