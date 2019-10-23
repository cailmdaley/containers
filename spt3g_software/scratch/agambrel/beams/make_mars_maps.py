import argparse as ap
import numpy as np
import scipy.stats
import os, numpy
from glob import glob
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g import core, std_processing, mapmaker, calibration, todfilter, \
coordinateutils, timestreamflagging

'''
This script is mostly copied from make_coadd.py 
(in spt3g_software/calibration/scripts/fluxpointcal), which does the
auto-processed point source maps. It adds more options than that script.
Can run locally on scott/amundsen for one or two Mars observations, but 
if making several, submit to the grid using condor_submit_mars_maps.py
in this directory.

If running locally, input_files can just be a directory that contains 
the observation and cal frame.
'''

P = ap.ArgumentParser(description='Make Mars maps with special options',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('--cal-file', action='store', default=None,
               help='Calibration file if want to override fist input file')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-s', '--source', action='store', 
               default='mars', help='name of source')
P.add_argument('--poly-order', action='store', type=int,
               default=3, help='polynomial filter order')
P.add_argument('-k', '--source-relative', action='store_true',
               default=False, help='calibrate in source-relative units rather '+
               'than K_cmb')
P.add_argument('-r', '--res', action='store', 
               default=0.2, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', 
               default=4, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', 
               default=4, help='map height [deg]')
P.add_argument('--wafer-split', action='store_true', default=False,
               help='Split maps by wafer')
P.add_argument('--lr', action='store_true', default=False,
               help='Split maps by scan direction')
P.add_argument('--poly1-map', action='store', default=None,
               help='Map to use to make point source mask for filter/weight')
P.add_argument('--flag-saturated', action='store_true', default=False,
               help='Flag detectors that saturate on Mars.')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# copied from Sam's script
waflist = ['w203', 'w201', 'w188', 'w174', 'w177', 'w176', 'w172', 'w180', 
           'w181', 'w187']

print('input files: {}'.format(args.input_files))

# If input_files is a directory, search it for data
# Use cal file in that directory if cal_file arg isn't specified.
if os.path.isdir(args.input_files[0]):
    print(args.input_files[0])
    obsid = os.path.basename(args.input_files[0])
    print('obsid', obsid)
    # calibration must be first argument!
    if args.cal_file is None:
        cal = glob(os.path.join(args.input_files[0], '*offline_calibration.g3'))
    else:
        cal = [args.cal_file]
    dat = sorted(glob(os.path.join(args.input_files[0], '[0-9]*.g3')))
    args.input_files = cal + dat
else:
    obsid = os.path.basename(args.input_files[0]).split('.')[0]

# Grab the observation time -- need it for planets
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

if args.poly1_map is not None:
    # Find poly1 map's brightest point to mask it
    poly1_map = core.G3File(os.path.join(args.poly1_map))
    xvals = []
    yvals = []
    for fr in auto_map:
        if fr.type == core.G3FrameType.Map:
            tmap = np.asarray(fr['T'])
            max_y, max_x = np.unravel_index(tmap.argmax(), tmap.shape)
            xvals.append(max_x)
            yvals.append(max_y)
    # median brightest point for three fequency maps
    xloc = int(np.median(xvals))
    yloc = int(np.median(yvals))
    print(xloc, yloc)
    mask_center = fr['T'].pixel_to_angle(xloc, yloc)
    print(mask_center)
    radius = 4.*core.G3Units.arcmin
else:
    mask_center = [0, 0]
    radius = 0

# Generate map stub
# This sets up a grid of pixels around where the planet is expected to be
# at the beginning of the observation
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len=float(args.xlen)*core.G3Units.deg/res,
    y_len=float(args.ylen)*core.G3Units.deg/res, res=res,
    proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time=starttime)

# Generate point source mask
# this is used for 1) fitting the polynomial coefficients in the filtering, 
# and 2) computing the inverse variance of the filtered timestream for 
# weighting
if args.poly1_map is not None:
    psmask = std_processing.CreateSourceMapStub(
        args.source, x_len=float(args.xlen)*core.G3Units.deg/res,
        y_len=float(args.ylen)*core.G3Units.deg/res, res=res,
        proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        at_time=starttime)

    ra, dec = list(mask_center) 
    print('masking ra {} dec {}'.format(ra, dec))
    mapmaker.fill_point_source_mask_flat_map([ra], [dec], [radius], False, 
                                             False, psmask)
    ps_map_id = 'PointSourceMap'
else:
    # there is not masking if there is no poly1 map
    ps_map_id = None

# Begin pipeline
pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)
# Deal with partially-complete calibration frames made during autoprocessing
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(FillCoordTransRotations,
         transform_store_key='OnlineRaDecRotation',
         bs_az_key='RawBoresightAz', bs_el_key='RawBoresightEl',
         bs_ra_key='OnlineBoresightRa', bs_dec_key='OnlineBoresightDec',
         do_bad_transform=True)

# do pre-calibration flagging
pipe.Add(std_processing.flagsegments.FlagInvalidData, flag_key='Flags', 
         ts_key='RawTimestreams_I')
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key='Flags')
pipe.Add(std_processing.flagsegments.FlagUncalibratable, flag_key='Flags',
         ts_key='RawTimestreams_I')
if args.flag_saturated:
    # Remove already flagged scans
    pipe.Add(timestreamflagging.RemoveFlagged, 
             input_ts_key='RawTimestreams_I',
             input_flag_key='Flags', output_ts_key='DeflaggedTimestreamsFirst')
    # On remaining good scans, flag saturated bolos
    pipe.Add(std_processing.flagsegments.FlagSaturatedBolosMars, 
             ts_key='RawTimestreams_I', flag_key='Flags')
    # Remove the newly flagged scans
    pipe.Add(timestreamflagging.RemoveFlagged, 
             input_ts_key='DeflaggedTimestreamsFirst',
             input_flag_key='Flags', output_ts_key='DeflaggedTimestreams')
else:
    pipe.Add(timestreamflagging.RemoveFlagged, 
             input_ts_key='RawTimestreams_I',
             input_flag_key='Flags', output_ts_key='DeflaggedTimestreams')

# Apply calibration
pipe.Add(std_processing.CalibrateRawTimestreams, 
         units=core.G3TimestreamUnits.Power, output='TimestreamsWatts')
pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
         OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', 
         Output='DeflaggedCalTimestreams')

if args.poly1_map is not None:
    pipe.Add(mapmaker.MapInjector, map_id=ps_map_id, maps_lst=[psmask], 
             is_stub=False)

# This is a blank map used for boresight pointing
pipe.Add(mapmaker.MapInjector, map_id='bsmap', maps_lst=[smstub,], 
         is_stub=False)
pipe.Add(mapmaker.CalculatePointing, map_id='bsmap',
         pointing_store_key='PixelPointing', 
         ts_map_key='DeflaggedCalTimestreams', 
         trans_key='OnlineRaDecRotation')

# Basic timestream filtering
# Would use dynamic_source_filter to flag point sources, but S. Guns 
# says it doesn't work so use fixed poly1 point source mask instead
pipe.Add(mapmaker.TodFiltering, ts_in_key='DeflaggedCalTimestreams',
         ts_out_key='PolyFilteredTimestreams', 
         use_dynamic_source_filter=False, poly_order=args.poly_order, 
         point_source_mask_id=ps_map_id, 
         point_source_pointing_store_key='PixelPointing',
         filter_mask_key='FilterMask')

# Standard bolometer weighting
if args.poly1_map is not None:
    # Mask the point source when calculating weight
    pipe.Add(std_processing.weighting.AddMaskedVarWeight, 
             input='PolyFilteredTimestreams', output='TodWeights')
else:
    pipe.Add(std_processing.weighting.AddSigmaClippedWeight,
             input='PolyFilteredTimestreams', output='TodWeights')

# Post-calibration flags
pipe.Add(timestreamflagging.noiseflagging.FlagUnphysicallyLowVariance, 
         ts_key='PolyFilteredTimestreams') 
pipe.Add(timestreamflagging.RemoveFlagged,  
         input_ts_key='PolyFilteredTimestreams', input_flag_key='Flags',
         output_ts_key='DeflaggedPFTimestreams')
pipe.Add(timestreamflagging.GenerateFlagStats, flag_key='Flags')

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 
                            'TimestreamsWatts', 'CalTimestreams', 
                            'DeflaggedCalTimestreams', 
                            'PolyFilteredTimestreams'])
pipe.Add(core.Dump)

# Ugly code to make any combination of splits possible
if args.lr and args.wafer_split:
    pipe.Add(std_processing.pointing.split_left_right_scans, 
             ts_in_key='DeflaggedPFTimestreams')
    for direction in ['Left', 'Right']:
        pipe.Add(calibration.SplitByWafer, 
                 input='DeflaggedPFTimestreams{}'.format(direction))
        for wafer in waflist:
            wafer = wafer.upper()
            pipe.Add(calibration.SplitByBand, 
                     input='DeflaggedPFTimestreams'+direction+wafer,
                     output_root='DeflaggedPFTimestreams'+direction+wafer)
                
    for direction in ['Left', 'Right']:
        for wafer in waflist:
            wafer = wafer.upper()
            for band in ['90GHz', '150GHz', '220GHz']:
                pipe.Add(
                    mapmaker.MapInjector, 
                    map_id=args.source+'{}-{}-{}'.format(
                        direction, wafer, band),
                    maps_lst=[smstub], is_stub=True, make_polarized=False, 
                    do_weight=True)
                pipe.Add(
                    mapmaker.BinMap, 
                    map_id=args.source+'{}-{}-{}'.format(direction, wafer, 
                                                         band),
                    ts_map_key='DeflaggedPFTimestreams{}{}{}'.format(direction,
                                                                     wafer,
                                                                     band),
                    pointing_store_key='PixelPointing', 
                    timestream_weight_key='TodWeights',
                    trans_key='OnlineRaDecRotation')

elif args.lr:
    pipe.Add(std_processing.pointing.split_left_right_scans, 
             ts_in_key='DeflaggedPFTimestreams')
    for direction in ['Left', 'Right']:
        pipe.Add(calibration.SplitByBand, 
                 input='DeflaggedPFTimestreams'+direction,
                 output_root='DeflaggedPFTimestreams'+direction)
                
    for direction in ['Left', 'Right']:
        for band in ['90GHz', '150GHz', '220GHz']:
            pipe.Add(
                mapmaker.MapInjector, 
                map_id=args.source+'{}-{}'.format(
                    direction, band),
                maps_lst=[smstub], is_stub=True, make_polarized=False, 
                do_weight=True)
            pipe.Add(
                mapmaker.BinMap, 
                map_id=args.source+'{}-{}'.format(direction, band),
                ts_map_key='DeflaggedPFTimestreams{}{}'.format(direction,
                                                               band),
                pointing_store_key='PixelPointing', 
                timestream_weight_key='TodWeights',
                trans_key='OnlineRaDecRotation')

elif args.wafer_split:
    pipe.Add(calibration.SplitByWafer, 
             input='DeflaggedPFTimestreams')
    for wafer in waflist:
        wafer = wafer.upper()
        pipe.Add(calibration.SplitByBand, 
                 input='DeflaggedPFTimestreams'+wafer,
                 output_root='DeflaggedPFTimestreams'+wafer)
                
    for wafer in waflist:
        wafer = wafer.upper()
        for band in ['90GHz', '150GHz', '220GHz']:
            pipe.Add(
                mapmaker.MapInjector, 
                map_id=args.source+'{}-{}'.format(wafer, band),
                maps_lst=[smstub], is_stub=True, make_polarized=False, 
                do_weight=True)
            pipe.Add(
                mapmaker.BinMap, 
                map_id=args.source+'{}-{}'.format(wafer, band),
                ts_map_key='DeflaggedPFTimestreams{}{}'.format(wafer, band),
                pointing_store_key='PixelPointing', 
                timestream_weight_key='TodWeights',
                trans_key='OnlineRaDecRotation')

else:
    pipe.Add(calibration.SplitByBand, 
             input='DeflaggedPFTimestreams',
             output_root='DeflaggedPFTimestreams')
    for band in ['90GHz', '150GHz', '220GHz']:
        pipe.Add(
            mapmaker.MapInjector, 
            map_id=args.source+'{}'.format(band),
            maps_lst=[smstub], is_stub=True, make_polarized=False, 
            do_weight=True)
        pipe.Add(
            mapmaker.BinMap, 
            map_id=args.source+'{}'.format(band),
            ts_map_key='DeflaggedPFTimestreams{}'.format(band),
            pointing_store_key='PixelPointing', 
            timestream_weight_key='TodWeights',
            trans_key='OnlineRaDecRotation')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()


