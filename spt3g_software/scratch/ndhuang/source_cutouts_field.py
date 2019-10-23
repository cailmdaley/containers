'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import os
import argparse as ap
import numpy as np
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging
import spt3g.mapmaker.mapmakerutils as maputils
import nhutils as nhu
import scipy.stats


# Usage: makecoadd.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
           help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
               required = True, help='Output filename')
P.add_argument('-s', '--source', action='store', 
           default=None, help='name of source')
P.add_argument('-x', '--xlen', action='store', type=float,
           default=1, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
           default=1, help='map height [deg]')
args = P.parse_args()
x_len = args.xlen
y_len = args.ylen
res = .4 * core.G3Units.arcmin

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

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

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
pipe.Add(calibration.ApplyTCalibration, InKCMB=not False,
    OpacityCorrection=not False, Input='TimestreamsWatts', Output='CalTimestreams')

# Basic timestream filtering with dynamic source filter to handle bright point sources with
# unreliable pointing
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreams',
         ts_out_key = 'PolyFilteredTimestreams', 
         use_dynamic_source_filter = True,
         poly_order = 19)

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

# Which maps to make?
cat = nhu.at20_sources_for_field([0, float(args.source.split('dec')[-1])],
                                 [100, 7.5], 2)
if len(cat) < 1:
    raise RuntimeError('No sources')
inds = np.argsort(cat.flux)[::-1]
# The maps should be roughly sorted by increasing flux
for src in cat[inds]:
    # Generate map stub
    sourcename = '{:.1f},{:.1f}'.format(src.ra, src.dec)
    smstub = maputils.FlatSkyMap(x_len = int(x_len * core.G3Units.deg / res), 
                                 y_len = int(y_len * core.G3Units.deg / res), 
                                 res = res, proj = coordinateutils.MapProjection.Proj0,
                                 alpha_center = src.ra * core.G3Units.deg,
                                 delta_center = src.dec * core.G3Units.deg,
                                 pol_type = core.MapPolType.T,
                                 coord_ref = core.MapCoordReference.Equatorial)
    # Kick off maps
    for band in ['90', '150', '220']: # XXX should be automatic
        for direction in ['Left', 'Right']:
            pipe.Add(mapmaker.MapInjector, 
                     map_id = '{}_{}GHz{}'.format(sourcename, band, direction),
                     maps_lst = [smstub], is_stub = True, make_polarized = False, 
                     do_weight = True)

    for direction in ['Left', 'Right']:
        pipe.Add(mapmaker.CalculatePointing, 
                 map_id = '{}_{}GHz{}'.format(sourcename, band, direction),
                 pointing_store_key = 'PixelPointing' + direction + sourcename,
                 ts_map_key = 'DeflaggedTimestreams',
                 trans_key = 'OnlineRaDecRotation')
    
        # Don't need warnings about missing scan frames
        core.set_log_level(core.G3LogLevel.LOG_ERROR, unit = 'MapBinner')
    for band in ['90', '150', '220']:
        for direction in ['Left', 'Right']:
            pipe.Add(mapmaker.BinMap, 
                     map_id = '{}_{}GHz{}'.format(sourcename, band, direction),
                     ts_map_key = 'DeflaggedTimestreams{}{}GHz'.format(direction, 
                                                                       band),
                     pointing_store_key = 'PixelPointing' + direction + sourcename, 
                     timestream_weight_key = 'TodWeights',
                     trans_key = 'OnlineRaDecRotation')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename = args.output)
pipe.Run()
