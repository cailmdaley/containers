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
P.add_argument('-r', '--res', action='store', type=float,
           default=0.5, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float,
               default=1, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
               default=1, help='map height [deg]')
args = P.parse_args()

res = .2 * core.G3Units.arcmin

class fix_bp(object):
    def __init__(self):
        self.buffer = []
        pass
    
    def __call__(self, fr):
        if fr.type != core.G3FrameType.Calibration:
            if self.buffer is None:
                return None
            self.bp_frame['BolometerProperties'] = self.bp
            buf = self.buffer
            self.buffer = None
            return buf
        if 'BolometerProperties' in fr:
            self.bp = fr.pop('BolometerProperties', None)
            self.bp_frame = fr
        if 'PointingOffsetX' in fr:
            for bolo in self.bp.keys():
                try:
                    self.bp[bolo].x_offset = fr['PointingOffsetX'][bolo]
                    self.bp[bolo].y_offset = fr['PointingOffsetY'][bolo]
                except KeyError:
                    pass
        else:
            self.buffer.append(fr)
        return []

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
def flip_xy(fr):
    if 'BolometerProperties' not in fr:
        return
    bp = fr.pop('BolometerProperties', None)
    for b, prop in bp:
        prop.y_offset *= -1
        prop.x_offset *= -1
    fr['BolometerProperties'] = bp
    return fr
# pipe.Add(flip_xy)

pipe.Add(fix_bp)

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

# Next to source-relative units (XXX: hardcode list of sources here), without the opacity
# correction this is likely meant to find (XXX: what about other uses?)
# Stop CalibrateRawTimestreams at Watts to (optionally) avoid the opacity correction.
pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Power,
    output='TimestreamsWatts')
pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
    OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')

# Basic timestream filtering with dynamic source filter to handle bright point sources with
# unreliable pointing
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreams',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = True,
    poly_order = 4)

# Standard bolometer weighting
pipe.Add(std_processing.weighting.AddSigmaClippedWeight)

# do some very nice flagging
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
pipe.Add(timestreamflagging.FlagBadHousekeeping, 
         ts_key = 'PolyFilteredTimestreams')
pipe.Add(timestreamflagging.flaggingutils.SigmaclipFlagGroupG3MapValue,
         m_key = 'TodWeights', low = 3.0, high = 3.0, flag_reason = 'BadWeight')
# remove flagged detectors
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')

# Split by direction
pipe.Add(std_processing.pointing.split_left_right_scans)

# Split by band
for direction in ['Left', 'Right']:
    pipe.Add(calibration.SplitByBand, input='DeflaggedTimestreams' + direction,
             output_root='DeflaggedTimestreams' + direction)

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams', 'PolyFilteredTimestreams', 'DeflaggedTimestreamsLeft', 'DeflaggedTimestreamsRight'])

# pipe.Add(core.Dump)

# Kick off maps
for band in ['90', '150', '220']: # XXX should be automatic
    for direction in ['Left', 'Right']:
        pipe.Add(mapmaker.MapInjector, 
                 map_id=args.source + '_%sGHz' % band + direction,
                 maps_lst=[smstub], is_stub=True, make_polarized=False, 
                 do_weight=True)

for direction in ['Left', 'Right']:
    pipe.Add(mapmaker.CalculatePointing, 
             map_id=args.source + '_150GHz' + direction,
             pointing_store_key = 'PixelPointing' + direction,
             ts_map_key = 'DeflaggedTimestreams', trans_key='OnlineRaDecRotation')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Scan)
def dropshit(fr):
    for key in fr:
        if 'PixelPointing' not in key:
            fr.pop(key, None)
pipe.Add(dropshit)
pipe.Add(core.G3WRiter, filename = args.output)

'''
# Don't need warnings about missing scan frames
core.set_log_level(core.G3LogLevel.LOG_ERROR, unit = 'MapBinner')
for band in ['90', '150', '220']:
    for direction in ['Left', 'Right']:
        pipe.Add(mapmaker.BinMap, 
                 map_id=args.source + '_%sGHz' % band + direction,
                 ts_map_key='DeflaggedTimestreams{}{}GHz'.format(direction, band),
                 pointing_store_key='PixelPointing' + direction, 
                 timestream_weight_key = 'TodWeights',
                 trans_key='OnlineRaDecRotation')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
'''
pipe.Run()
