'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging, pointing
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
P.add_argument('-p', '--polarized', action='store_true',
           default=False, help='make polarized maps')
P.add_argument('-r', '--res', action='store', type=float,
           default=0.5, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float,
           default=3, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
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

# Generate map stub and mask
mapkwargs = {'x_len': float(args.xlen)*core.G3Units.deg/res,
             'y_len': float(args.ylen)*core.G3Units.deg/res, 
             'res': res,
             'proj': coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
             'at_time': starttime}
smstub = std_processing.CreateSourceMapStub(args.source, **mapkwargs)
mask = std_processing.CreateSourceMask(args.source, **mapkwargs)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

# Deal with partially-complete calibration frames made during autoprocessing
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Add Online pointing to scans.
# Clean up pre-existing timestreams
pipe.Add(
    core.Delete,
    keys=['OnlineBoresightAz', 'OnlineBoresightEl',
          'OnlineBoresightRa', 'OnlineBoresightDec', 'OnlineRaDecRotation']
)
transform_key = 'OnlineRaDecRotation'
pipe.Add(
    pointing.CalculateCoordTransRotations,
    raw_az_key='RawBoresightAz',
    raw_el_key='RawBoresightEl',
    output='OnlineBoresight',
    transform_store_key=transform_key,
    model='OnlinePointingModel',
    flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
)

# do some very nice early flagging
pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
         ts_key='RawTimestreams_I', flag_on_delta_rfrac=True)

# Next to source-relative units (XXX: hardcode list of sources here), without the opacity
# correction this is likely meant to find (XXX: what about other uses?)
# Stop CalibrateRawTimestreams at Watts to (optionally) avoid the opacity correction.
pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Power,
    output='TimestreamsWatts')
pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
    OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')

# Basic timestream filtering with source filter to handle bright point sources
# The mask is fairly large to avoid issues with online pointing jitter
pipe.Add(mapmaker.MapInjector, map_id = 'mask', maps_lst = [mask],
         is_stub = False)
pipe.Add(mapmaker.mapmakerutils.CalculatePointing, 
         map_id = 'mask', 
         pointing_store_key = 'PixelPointing',
         trans_key = transform_key,
         ts_map_key = 'CalTimestreams')

pipe.Add(
    mapmaker.TodFiltering, ts_in_key = 'CalTimestreams',
    ts_out_key = 'PolyFilteredTimestreams', 
    point_source_mask_id = 'mask',
    point_source_pointing_store_key = 'PixelPointing', poly_order = 4)

# Standard bolometer weighting
pipe.Add(std_processing.weighting.AddMaskedVarWeight)

# do some very nice flagging
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
pipe.Add(timestreamflagging.FlagNaNs, ts_key = 'PolyFilteredTimestreams', 
         flag_key = 'Flags')
pipe.Add(timestreamflagging.flaggingutils.SigmaclipFlagGroupG3MapValue,
         m_key = 'TodWeights', low = 3.0, high = 3.0, flag_reason = 'BadWeight')

stats = timestreamflagging.GenerateFlagStats(flag_key='Flags')
pipe.Add(stats)

# remove flagged detectors
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')

# Split by band
pipe.Add(calibration.SplitByBand, input='DeflaggedTimestreams',
    output_root='DeflaggedTimestreams')

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams', 'PolyFilteredTimestreams'])

pipe.Add(std_processing.AddMetaData, stats = stats, args = vars(args))

pipe.Add(core.Dump)

# Kick off maps
for band in ['90', '150', '220']: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id=args.source + '-%sGHz' % band,
      maps_lst=[smstub], is_stub=True, make_polarized=args.polarized, do_weight=True)

for band in ['90', '150', '220']: # XXX should be automatic
    pipe.Add(mapmaker.BinMap, map_id=args.source + '-%sGHz' % band,
             ts_map_key='DeflaggedTimestreams%sGHz' % band,
             pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights',
             trans_key=transform_key)

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map
  or fr.type == core.G3FrameType.PipelineInfo) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()
