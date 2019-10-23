'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils, timestreamflagging
import scipy.stats
import os, numpy, sys

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
P.add_argument('-p', '--poly_order', type=int, default=4)
P.add_argument('-d', '--dynamic_mask', action='store_true')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin
wafer_lst = ['W148', 'W158', 'W157', 'W147', 'W153', 'W152', 'W142', 'W139', 'W162', 'W136']

#std_processing.weighting.AddPSDWeights(self, low_f=3.0000000000000004e-08, high_f=5e-08, input='PolyFilteredTimestreams')

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

source_list = ['0537-441']
source_ra = {'0537-441': 84.703*core.G3Units.deg}
source_dec = {'0537-441': -44.09*core.G3Units.deg}
source_radius = {'0537-441': 2.0*core.G3Units.arcmin}

#[84.703, -44.0977]
#84.6984, -44.091
# Generate point source mask

psmask = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)
if args.source in source_list:
    ra = source_ra[args.source]
    dec = source_dec[args.source]
    radius = source_radius[args.source]
else:
    print("supported sources for static masking:", source_list)
    sys.exit()
mapmaker.fill_point_source_mask_flat_map([ra], [dec], [radius], True, False, psmask)
# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files, n_frames_to_read=10)
#pipe.Add(core.G3Reader, filename=args.input_files,n_frames_to_read=10)

# Deal with partially-complete calibration frames made during autoprocessing
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Cut uncalibratable detectors early
pipe.Add(calibration.elnod_analysis.RotateIQ, i_rotate_key='RawTimestreamsRotated')

pipe.Add(std_processing.flagsegments.FlagInvalidDataStrict, 
         overbiased_rfrac_cutoff = 0.97,
         flag_key='Flags',
         ts_key='RawTimestreamsRotated')
#pipe.Add(timestreamflagging.RemoveFlaggedTimestreams, 
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key='RawTimestreamsRotated', 
         input_flag_key='Flags', 
         output_ts_key='FilteredTimestreams')

pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='FilteredTimestreams', 
         output='TimestreamsWatts', units=core.G3TimestreamUnits.Power, 
         ignore_missing=True)

pipe.Add(calibration.ApplyTCalibration, Source='RCW38', 
         InKCMB=not args.source_relative,
         OpacityCorrection=not args.source_relative, 
         Input='TimestreamsWatts', Output='CalTimestreams')

pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key='Flags')
pipe.Add(std_processing.flagsegments.FlagNoisyOrNonGaussianData, 
         flag_key='Flags',ts_key='CalTimestreams',
         max_derivative_val=2e-7,
         min_variance=1e-4,
         max_variance=1e-2,
         glitch_num_above = [1],
         glitch_thresholds = [10],
         plot_statistics=True, 
         plot_directory='/home/nlharr/spt3g_software/scratch/nlharr/focus/plots/')

pipe.Add(timestreamflagging.RemoveFlagged, input_ts_key='CalTimestreams', 
         input_flag_key='Flags', output_ts_key='FlaggedCalTimestreams')
pipe.Add(timestreamflagging.GenerateFlagStats, flag_key='Flags')
# Inject static point source mask and calculate pointing


pipe.Add(mapmaker.MapInjector, map_id = 'PointSourceMask', 
         maps_lst = [psmask], is_stub=False)
pipe.Add(mapmaker.CalculatePointing, map_id = 'PointSourceMask',
         is_healpix = False, pointing_store_key = 'PixelPointing',
         ts_map_key = 'FlaggedCalTimestreams',
         boresight_ra_key='OnlineBoresightRa', 
         boresight_dec_key='OnlineBoresightDec')

if args.dynamic_mask:
    ps_mask_id = None
else:
    ps_mask_id = 'PointSourceMask'
# Timestream filtering with point source mask
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'FlaggedCalTimestreams',
         ts_out_key = 'PolyFilteredTimestreams', 
         use_dynamic_source_filter = args.dynamic_mask,
         poly_order = args.poly_order, 
         point_source_mask_id = ps_mask_id,
         #point_source_mask_id = None,
         point_source_pointing_store_key = 'PixelPointing')

#pipe.Add(core.Dump)

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 
                            'FilteredTimestreams', 'TimestreamsWatts', 
                            'CalTimestreams', 'FlaggedCalTimestreams'])

pipe.Add(std_processing.weighting.AddPSDWeights, input='PolyFilteredTimestreams')
#pipe.Add(std_processing.weighting.AddMadVarWeight, input='PolyFilteredTimestreams')

# pipe.Add(core.Dump)
pipe.Add(calibration.SplitByWafer, in_ts_key='PolyFilteredTimestreams')

# Kick off maps
for wafer in wafer_lst:
    pipe.Add(mapmaker.MapInjector, map_id = wafer,
             maps_lst = [smstub], is_stub = True, 
             make_polarized = False, do_weight = True)
    pipe.Add(mapmaker.BinMap, map_id=wafer,
             ts_map_key='PolyFilteredTimestreams'+wafer,
             pointing_store_key='PixelPointing',
             timestream_weight_key = 'TodWeights')

    pipe.Add(calibration.SplitTimestreamsByBand, input='PolyFilteredTimestreams'+wafer,
             output_root='PolyFilteredTimestreams'+wafer)
    for band in ['90GHz', '150GHz', '220GHz']:
        pipe.Add(mapmaker.MapInjector, map_id = wafer+band,
                 maps_lst = [smstub], is_stub = True, 
                 make_polarized = False, do_weight = True)
        pipe.Add(mapmaker.BinMap, map_id=wafer+band,
                 ts_map_key='PolyFilteredTimestreams'+wafer+band,
                 pointing_store_key='PixelPointing',
                 timestream_weight_key = 'TodWeights')

pipe.Add(calibration.SplitTimestreamsByBand, input='PolyFilteredTimestreams',
         output_root='PolyFilteredTimestreams')
for band in ['90GHz', '150GHz', '220GHz']:
    pipe.Add(mapmaker.MapInjector, map_id = band,
             maps_lst = [smstub], is_stub = True, 
             make_polarized = False, do_weight = True)
    pipe.Add(mapmaker.BinMap, map_id=band,
             ts_map_key='PolyFilteredTimestreams'+band,
             pointing_store_key='PixelPointing',
             timestream_weight_key = 'TodWeights')



pipe.Add(core.Dump)
pipe.Add(mapmaker.MapInjector, map_id = 'SumMap',
             maps_lst = [smstub], is_stub = True, 
             make_polarized = False, do_weight = True)
pipe.Add(mapmaker.BinMap, map_id='SumMap',
             ts_map_key='PolyFilteredTimestreams',
             pointing_store_key='PixelPointing',
             timestream_weight_key = 'TodWeights')


pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()

