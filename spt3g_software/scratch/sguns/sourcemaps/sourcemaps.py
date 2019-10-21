'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging
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
P.add_argument('-p', '--poly_order', type=int, action='store',
           default=4, help='poly order to use for filtering')
P.add_argument('-sm', '--static_mask', action='store_true', 
           default=False, help='use static mask')
P.add_argument('-dm', '--dynamic_mask', action='store_true',
           default=False, help='use dynamic mask')
P.add_argument('-r', '--res', action='store', 
           default=0.5, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', 
           default=3, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', 
           default=3, help='map height [deg]')
P.add_argument('-ra', '--right_ascension', action='store',
           default = None, help = '(custom) right ascension [deg]')
P.add_argument('-dec', '--declination', action='store',
           default = None, help = '(custom) source declination [deg]')
P.add_argument('-rad', '--radius', action='store',
           default = 2.0, help = 'source radius [arcmin]')
P.add_argument('-pw', '--per_wafer', action='store_true',
           default=False, help = 'make per-wafer maps')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin
wafer_lst = ['w203', 'w201', 'w188', 'w174', 'w177', 'w176', 'w172', 'w180', 'w181', 'w187']

if args.static_mask and args.dynamic_mask:
    print "cannot use both static and dynamic mask. Use either -sm, -dm or none."
    sys.exit()

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

# Generate point source mask
psmask = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

ra, dec = std_processing.sourcemaps.get_source_ra_dec(args.source, starttime)
if args.right_ascension is not None:
    ra = float(args.right_ascension)*core.G3Units.deg
if args.declination is not None:
    dec = float(args.declination)*core.G3Units.deg
radius = float(args.radius)*core.G3Units.arcmin

mapmaker.fill_point_source_mask_flat_map([ra], [dec], [radius], False, False, psmask)

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

#only needed for data before feb 2018
pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OnlineRaDecRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
         bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
         do_bad_transform = True)


# Next to source-relative units (XXX: hardcode RCW38 here), without the opacity
# correction this is likely meant to find (XXX: what about other uses?)
# Stop CalibrateRawTimestreams at Watts to (optionally) avoid the opacity correction.
pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Power,
    output='TimestreamsWatts')
pipe.Add(calibration.ApplyTCalibration, Source='RCW38', InKCMB=not args.source_relative,
    OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')

# do some very nice flagging
pipe.Add(std_processing.flagsegments.FlagInvalidData, flag_key = 'Flags', ts_key = 'CalTimestreams')
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
pipe.Add(std_processing.flagsegments.FlagUncalibratable, ts_key = 'CalTimestreams', flag_key = 'Flags')
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'CalTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedCalTimestreams')

pipe.Add(mapmaker.MapInjector, map_id = 'PointSourceMask', maps_lst = [psmask], is_stub=False)
pipe.Add(mapmaker.CalculatePointing, map_id='PointSourceMask',
         pointing_store_key = 'PixelPointing',
         ts_map_key = 'DeflaggedCalTimestreams', trans_key='OnlineRaDecRotation')

if args.static_mask:
    ps_mask_id = 'PointSourceMask'
else:
    ps_mask_id = None

# Basic timestream filtering with dynamic source filter to handle bright point sources with
# unreliable pointing
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'DeflaggedCalTimestreams',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = args.dynamic_mask,
    poly_order = args.poly_order, point_source_mask_id = ps_mask_id, point_source_pointing_store_key = 'PixelPointing',
    filter_mask_key = 'FilterMask')

# Standard bolometer weighting
pipe.Add(std_processing.weighting.AddMaskedVarWeight, input='PolyFilteredTimestreams', output='TodWeights')
pipe.Add(timestreamflagging.noiseflagging.FlagUnphysicallyLowVariance, ts_key = 'PolyFilteredTimestreams') 
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedPFTimestreams')
pipe.Add(timestreamflagging.GenerateFlagStats, flag_key='Flags')

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'TimestreamsWatts', 'CalTimestreams', 'DeflaggedCalTimestreams', 'PolyFilteredTimestreams'])

# Split timestreams by wafer
if args.per_wafer:
    pipe.Add(calibration.SplitByWafer, in_ts_key='DeflaggedPFTimestreams')

    # Kick off maps
    for wafer in wafer_lst:
        pipe.Add(mapmaker.MapInjector, map_id = wafer,
                 maps_lst = [smstub], is_stub = True, 
                 make_polarized = False, do_weight = True)
        pipe.Add(mapmaker.BinMap, map_id=wafer,
                 ts_map_key='DeflaggedPFTimestreams'+wafer,
                 pointing_store_key='PixelPointing',
                 timestream_weight_key = 'TodWeights',trans_key = 'OnlineRaDecRotation')

        pipe.Add(calibration.SplitTimestreamsByBand, input='DeflaggedPFTimestreams'+wafer,
                 output_root='DeflaggedPFTimestreams'+wafer)
        for band in ['90GHz', '150GHz', '220GHz']:
            pipe.Add(mapmaker.MapInjector, map_id = wafer+band,
                     maps_lst = [smstub], is_stub = True, 
                     make_polarized = False, do_weight = True)
            pipe.Add(mapmaker.BinMap, map_id=wafer+band,
                     ts_map_key='DeflaggedPFTimestreams'+wafer+band,
                     pointing_store_key='PixelPointing',
                     timestream_weight_key = 'TodWeights',trans_key = 'OnlineRaDecRotation')
# Split by band
pipe.Add(calibration.SplitTimestreamsByBand, input='DeflaggedPFTimestreams',
    output_root='DeflaggedPFTimestreams')

# Kick off maps
for band in ['90', '150', '220']: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id='%sGHz' % band,
      maps_lst=[smstub], is_stub=True, make_polarized=False, do_weight=True)
    pipe.Add(mapmaker.BinMap, map_id='%sGHz' % band,
             ts_map_key='DeflaggedPFTimestreams%sGHz' % band,
             pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights',
             trans_key='OnlineRaDecRotation')

pipe.Add(mapmaker.MapInjector, map_id = 'SumMap',
             maps_lst = [smstub], is_stub = True, 
             make_polarized = False, do_weight = True)
pipe.Add(mapmaker.BinMap, map_id='SumMap',
             ts_map_key='DeflaggedPFTimestreams',
             pointing_store_key='PixelPointing',
             timestream_weight_key = 'TodWeights',trans_key = 'OnlineRaDecRotation')
pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
#pipe.Add(core.Dump)
pipe.Run()
