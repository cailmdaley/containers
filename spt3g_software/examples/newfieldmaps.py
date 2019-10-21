'''
Script to make CMB field maps.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, todfilter, coordinateutils, timestreamflagging
from spt3g.pointing import offline_pointing as op
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations

import scipy.stats
import pdb
import os, numpy, copy

# Usage: python newfieldmaps.py <input files.g3> -o outputmaps.g3

P = ap.ArgumentParser(description='Maps for a CMB field (SPTpol 500d by default)',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[], help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
P.add_argument('-r', '--res', action='store', 
           default=1.0, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', 
           default=45, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', 
           default=25, help='map height [deg]')
P.add_argument('-s', '--source', action='store', 
           default=None, help='Name of source')

args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# Grab the observation time in case we need it for planets
starttime = None
for fname in args.input_files:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            starttime = frame['RawTimestreams_I'].start
            if args.source is None:
                args.source = frame['SourceName']
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

#pipe.Add(core.G3Reader, filename=args.input_files, n_frames_to_read=10)
pipe.Add(core.G3Reader, filename=args.input_files)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Flag uncalibratable, invalid data
pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
         flag_key = 'Flags', ts_key = 'RawTimestreams_I')

# Apply calibrations
pipe.Add(std_processing.CalibrateRawTimestreams,
         output = 'CalTimestreams')

# Add point source mask and calculate detector pointing

# Make an empty flatsky map for the source filtering
fs_stub = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)
# XXX example file, use one appropriate for your field!
ps_map = mapmaker.pointsourceutils.make_point_source_map(fs_stub, '/spt/user/javva/lowell/ptsrc_config_ra0hdec-57p5_both_50mJy.txt') 
pipe.Add(mapmaker.MapInjector, map_id='PointSourceMask',maps_lst=[fs_stub,], is_stub=False)

pipe.Add(FillCoordTransRotations, transform_store_key = 'RaDecTransform',
         do_bad_transform = True)
pipe.Add(mapmaker.mapmakerutils.CalculatePointing, map_id = 'PointSourceMask', 
         pointing_store_key = 'PixelPointing', trans_key='RaDecTransform',
         ts_map_key = 'CalTimestreams')

# Basic timestream filtering
# DPD: I've commented out the point source mask bits
#      in order to make this example run faster.
pipe.Add(mapmaker.TodFiltering, ts_in_key='CalTimestreams',
    ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=False,
    poly_order=4, #point_source_mask_id = 'PointSourceMask',
         point_source_pointing_store_key = 'PixelPointing',
         filters_are_ell_based = True, lpf_filter_frequency=6600,
         mhpf_cutoff=20, boresight_az_key='OnlineBoresightAz',
         boresight_el_key='OnlineBoresightEl')

# Calculate Weights
pipe.Add(std_processing.weighting.AddPSDWeights,
         input = 'PolyFilteredTimestreams', output = 'TodWeights')

# More flags for calibrator SN, glitches, etc
pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,
         flag_key = 'Flags', ts_key = 'PolyFilteredTimestreams')

# Remove flagged timestreams
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'MapTimestreams')

# Split by band (could do other things: wafer, bolos, etc.)
pipe.Add(calibration.SplitByBand, input='MapTimestreams',
    output_root='MapTimestreams')

pipe.Add(core.Dump)

# Kick off maps
for band in ['90', '150', '220']: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id=args.source + '-%sGHz' % band,
             maps_lst=[smstub], is_stub=True, make_polarized=True, do_weight=True)
    pipe.Add(mapmaker.mapmakerutils.BinMap, map_id=args.source + '-%sGHz' % band,
             ts_map_key='MapTimestreams%sGHz' % band,
             trans_key='RaDecTransform',
             pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights')
    
pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop everything but maps

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run(profile=True)


