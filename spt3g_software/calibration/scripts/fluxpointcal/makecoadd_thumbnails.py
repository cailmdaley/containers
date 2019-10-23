'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging, pointing,sources
import scipy.stats
import os, numpy


@core.cache_frame_data(type=core.G3FrameType.Map,pointing_model = 'OnlinePointingModel')
def inject_pointing_model(frame, pointing_model = None):
    if frame.type != core.G3FrameType.Map:
        return
    if pointing_model is not None:
        frame['OnlinePointingModel'] = pointing_model
    return



## Usage: makecoadd_thumbnails.py <input files.g3> -o outputmaps.g3 -s subfield 
##                                                 --point-source-file /path/to/pointsourcefile.txt
##                                                 -N numberofsources -b 90 150 220 --ra-lim -35 35
##                                                 --dec-lim -41 71

P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[], help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3', help='Output filename')
P.add_argument('-s', '--source', action='store', default=None, help='name of source')
P.add_argument('-k', '--source-relative', action='store_true', default=False, 
               help='calibrate in source-relative units rather than K_cmb')
P.add_argument('-p', '--polarized', action='store_true', default=True, help='make polarized maps')
P.add_argument('-r', '--res', action='store', type=float, default=0.25, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float, default=.25, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float, default=.25, help='map height [deg]')
P.add_argument('--point-source-file', action='store',
               default='./code/spt3g_software/sources/1500d_ptsrc_and_decrement_list.txt', 
               help='Location of point source list')
P.add_argument('-N', '--numsources', action='store', type=int, default=10, 
               help='Number of brightest sources in the subfield')
P.add_argument('-b', '--bands', action='store', nargs='+', default=['90'], help='Bands from which to make coadds')
P.add_argument('--ra-lim',action='store',type=float, nargs='+', default=[-35,35], 
               help='RA boundary for ptsrc searching [deg]')
P.add_argument('--dec-lim',action='store',type=float, nargs='+', default=[-41,-70], 
               help='Dec boundary for ptsrc searching [deg]')


args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin



bright_sources = sources.source_utils.get_brightest_sources(args.point_source_file, subfield=args.source, 
                                                            n=args.numsources, ra_lim=args.ra_lim, dec_lim=args.dec_lim)


# Generate map stub
smstubs   = {}

for source in bright_sources:
    RA =bright_sources[source]['RA']
    Dec=bright_sources[source]['Dec']
    smstubs[source] = coordinateutils.FlatSkyMap(x_len=int(float(args.xlen)*core.G3Units.deg/res),
        y_len=int(float(args.ylen)*core.G3Units.deg/res), res=res,
        alpha_center=RA*core.G3Units.deg, delta_center=Dec*core.G3Units.deg,
        proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        pol_type=coordinateutils.MapPolType.T, coord_ref=coordinateutils.MapCoordReference.Equatorial)
 
## hardcoded for now... should change that...
if args.source is not None:
    maskDec = float(args.source.split('ra0hdec')[-1])* core.G3Units.deg
    maskRA  = 0.0 * core.G3Units.deg
    mask_x_len = int(75*core.G3Units.deg/res)
    mask_y_len = int(7*core.G3Units.deg/res)
else:
    maskDec = numpy.mean(args.dec_lim)* core.G3Units.deg
    maskRA  = numpy.mean(args.ra_lim)* core.G3Units.deg
    mask_x_len = abs(int(numpy.diff(args.ra_lim)[0]*core.G3Units.deg/res))
    mask_y_len = abs(int(numpy.diff(args.dec_lim)[0]*core.G3Units.deg/res))
    
# Set up an empty map for point source filtering
ps_params = coordinateutils.FlatSkyMap(x_len=mask_x_len, y_len=mask_y_len, res=res,
        proj=coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
        alpha_center=maskRA, delta_center=maskDec, pol_type=coordinateutils.MapPolType.T,
        coord_ref=coordinateutils.MapCoordReference.Equatorial)

# Fill map with the point source mask
sources.source_utils.make_point_source_map(ps_params, args.point_source_file)
    
ps_mask_id = 'PointSourceMask'


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
pipe.Add(pointing.offline_pointing.ApplyPointingCorrection)

pipe.Add(
    pointing.CalculateCoordTransRotations,
    raw_az_key='RawBoresightAz',
    raw_el_key='RawBoresightEl',
    output='OfflineBoresight',
    transform_store_key='OfflineRaDecRotation',
    model='OfflinePointingModel',
    flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
    #, 'thermolin'] # Thermoline broken as of 4/13/17
)

# do some very nice early flagging
pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
         ts_key='RawTimestreams_I')

pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Power,
    output='TimestreamsWatts')
pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
    OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')

pipe.Add(mapmaker.MapInjector, map_id = ps_mask_id,
         maps_lst = [ps_params,], is_stub=False)

pipe.Add(mapmaker.CalculatePointing, map_id=ps_mask_id,
             pointing_store_key = 'PixelPointing',
             ts_map_key = 'CalTimestreams', trans_key='OfflineRaDecRotation')



# Basic timestream filtering with dynamic source filter to handle bright point sources with
# unreliable pointing
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreams',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = False,
    poly_order = 4, point_source_mask_id = ps_mask_id,point_source_pointing_store_key = 'PixelPointing', )

# Standard bolometer weighting
pipe.Add(std_processing.weighting.AddSigmaClippedWeight)

# do some very nice flagging
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
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
pipe.Add(core.Delete, keys=['PixelPointing','RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams', 'PolyFilteredTimestreams'])

pipe.Add(std_processing.AddMetaData, stats = stats, args = vars(args))


# Kick off maps
for source in smstubs:
    for band in args.bands: # XXX should be automatic
        pipe.Add(mapmaker.MapInjector, map_id=str(source)+ '-%sGHz' % band,
          maps_lst=[smstubs[source]], is_stub=True, make_polarized=args.polarized, do_weight=True)
    ## Calculate pointing for each source's thumbnail.
    pipe.Add(mapmaker.CalculatePointing, map_id=str(source) + '-%sGHz'%band,
             pointing_store_key = 'PixelPointing_%s'%source,
             ts_map_key = 'DeflaggedTimestreams', trans_key='OfflineRaDecRotation')
    for band in args.bands: # XXX should be automatic
        pipe.Add(mapmaker.BinMap, map_id=str(source) + '-%sGHz' % band,
                 ts_map_key='DeflaggedTimestreams%sGHz' % band,
                 pointing_store_key='PixelPointing_%s'%source, timestream_weight_key = 'TodWeights',
                 trans_key='OfflineRaDecRotation')
    pipe.Add(core.Delete, keys = ['PixelPointing_%s'%source])
    
pipe.Add(inject_pointing_model, pointing_model='OnlinePointingModel')

pipe.Add(lambda fr: fr.type is not core.G3FrameType.Scan)
#== core.G3FrameType.Map or fr.type == core.G3FrameType.PipelineInfo) # Drop TOD
#pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()
