'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
import scipy.stats
import os, numpy

# Usage: makecoadd.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
           help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
P.add_argument('-s', '--source', action='store', 
           default=None, help='name of source')
P.add_argument('-k', '--source-relative', action='store_true',
           default=False, help='calibrate in source-relative units rather than K_cmb')
P.add_argument('-r', '--res', action='store', 
           default=0.5, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', 
           default=3, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', 
           default=3, help='map height [deg]')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# Grab the observation time in case we need it for planets
starttime = None
for fname in args.input_files:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            starttime = frame['RawTimestreams_I'].start
            if args.source is None:
                # Grab source from data, removing any '-pixelraster' from
                # the end that indicates a fast point (every pixel rastered
                # over the source)
                args.source = frame['SourceName'].replace('-pixelraster', '')
                break
    if starttime is not None:
        break

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

# Deal with partially-complete calibration frames made during autoprocessing
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Cut uncalibratable detectors early
pipe.Add(calibration.elnod_analysis.RotateIQ, i_rotate_key='RawTimestreamsRotated')
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, input='RawTimestreamsRotated', output='FilteredTimestreams')

# Invert crosstalk and convert to Watts.
pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='FilteredTimestreams', output='TimestreamsWatts', units=core.G3TimestreamUnits.Power, ignore_missing=True)

# The equivalent thing that doesn't deal with crosstalk:
#pipe.Add(dfmux.ConvertTimestreamUnits, Input='FilteredTimestreams', Output='TimestreamsWatts', Units=core.G3TimestreamUnits.Power)

# Next to source-relative units (XXX: hardcode RCW38 here), without the opacity
# correction this is likely meant to find (XXX: what about other uses?)
pipe.Add(calibration.ApplyTCalibration, Source='RCW38', InKCMB=not args.source_relative,
    OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')

# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreams',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = True,
    poly_order = 4)

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams'])

# Quick and dirty weight calculation
class weights(object):
    def __init__(self):
        pass
    def __call__(self, fr):
        if 'CalibratorResponseSN' in fr:
            self.calsn = fr['CalibratorResponseSN']
        if 'PolyFilteredTimestreams' not in fr:
            return
        w = core.G3MapDouble()
        for k,ts in fr['PolyFilteredTimestreams'].iteritems():
            if not numpy.isfinite(ts).all():
                w[k] = 0
            elif (ts == 0).all():
                w[k] = 0
            elif k not in self.calsn or self.calsn[k] < 20:
                w[k] = 0
            else:
                w[k] = \
                  1./numpy.var(scipy.stats.sigmaclip(ts, 2.5, 2.5).clipped)
                if not numpy.isfinite(w[k]):
                    w[k] = 0
        fr['TodWeights'] = w
        # print(numpy.sum(numpy.asarray(w.values()) != 0))
pipe.Add(weights)

# pipe.Add(core.Dump)

# Kick off maps
pipe.Add(mapmaker.MapInjector, map_id = 'stub',
         maps_lst = [smstub], is_stub = True, 
         make_polarized = False, do_weight = True)
pipe.Add(mapmaker.CalculatePointing, map_id = 'stub',
  is_healpix = False, pointing_store_key = 'PixelPointing',
  ts_map_key = 'PolyFilteredTimestreams',
  boresight_ra_key='OnlineBoresightRa', boresight_dec_key='OnlineBoresightDec')
# Guess the list of bolos to use and other metadata for configuration
bolos = None
for fname in args.input_files:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            bolos = frame['RawTimestreams_I'].keys()
            starttime = frame['RawTimestreams_I'].start
            break
    if bolos is not None:
        break

pipe.Add(mapmaker.BinMap, map_id='stub',
      ts_map_key='PolyFilteredTimestreams',
         individual_bolos_to_map = bolos,
      pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

