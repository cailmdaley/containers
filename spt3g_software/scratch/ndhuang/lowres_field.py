'''
Script to make CMB field maps.
'''
import argparse
import sys
from spt3g import core, std_processing, mapmaker, calibration
from spt3g import timestreamflagging, todfilter, coordinateutils
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations


# Usage: makemaps.py <input files.g3> -o outputmaps.g3

@core.pipesegment
def field_mapmaker(pipe, split_leftright = False, map_params = None, 
                   ps_map = None, output_file = None, verbose = False):
    assert(map_params is not None and ps_map is not None)

    # Cut turnarounds
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
    # Flag junk
    pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
             flag_key = 'Flags', ts_key = 'RawTimestreams_I')
    # Apply calibrations
    pipe.Add(std_processing.CalibrateRawTimestreams, 
             units = core.G3TimestreamUnits.Power,
             output = 'TimestreamsWatts')
    # Turn opacity correction back on at some point
    pipe.Add(calibration.ApplyTCalibration, Source='RCW38', InKCMB = True, 
             OpacityCorrection = False, Input = 'TimestreamsWatts', 
             Output = 'CalTimestreams')
    # In case we're using pre-Feb 2018 data
    pipe.Add(FillCoordTransRotations, transform_store_key = "RaDecTransform", 
             do_bad_transform = True)
    # Calculate detector pointing early and add point source mask
    pipe.Add(mapmaker.MapInjector, map_id = 'PointSourceMask',
             maps_lst = [ps_map,], is_stub=False)
    pipe.Add(mapmaker.mapmakerutils.CalculatePointing, 
             map_id = 'PointSourceMask', 
             pointing_store_key = 'PixelPointing', trans_key='RaDecTransform',
             ts_map_key = 'CalTimestreams')
    # Basic timestream filtering
    pipe.Add(mapmaker.TodFiltering,
             # filtering options
             poly_order = 4, 
             filters_are_ell_based = True, 
             mhpf_cutoff = 50, lpf_filter_frequency = 6600,
             point_source_mask_id = 'PointSourceMask',
             # boiler plate
             ts_in_key='CalTimestreams',
             ts_out_key = 'PolyFilteredTimestreams', 
             point_source_pointing_store_key = 'PixelPointing',
             use_dynamic_source_filter = False,
             boresight_az_key='OnlineBoresightAz',
             boresight_el_key='OnlineBoresightEl')
    # Clean up detritus
    pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 
                                'TimestreamsWatts', 'CalTimestreams'])
    # Calculate Weights
    pipe.Add(std_processing.weighting.AddPSDWeights,
             low_f = .5 * core.G3Units.Hz,
             high_f = 8 * core.G3Units.Hz)
    # Flagging
    pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,
             flag_key = 'Flags', ts_key = 'PolyFilteredTimestreams')
    pipe.Add(timestreamflagging.flaggingutils.SigmaclipFlagGroupG3MapValue,
             m_key = 'TodWeights', low = 2.5, high = 2.5, 
             flag_reason = 'BadWeight')
    pipe.Add(timestreamflagging.GenerateFlagStats, flag_key = 'Flags')
    pipe.Add(timestreamflagging.RemoveFlagged, 
             input_ts_key = 'PolyFilteredTimestreams',
             input_flag_key = 'Flags',
             output_ts_key = 'DeflaggedTimestreams')
    if split_leftright:
        # split left and right scans if requested
        pipe.Add(std_processing.pointing.split_lef_right_scans)
        # Split by band (could do other things: wafer, bolos, etc.)
        for direction in ['Left', 'Right']:
            pipe.Add(calibration.SplitTimestreamsByBand, 
                     input='DeflaggedTimestreams' + direction,
                     output_root='DeflaggedTimestreams')
        if verbose:
            pipe.Add(core.Dump)
        for direction in ['Left', 'Right']:
            for band in ['90', '150', '220']: # XXX should be automatic
                mapid = '%s-%sGHz' %(band, direction),
                pipe.Add(mapmaker.MapInjector, map_id = mapid,
                         maps_lst=[map_params], is_stub=True, 
                         make_polarized=True, do_weight=True)
                pipe.Add(mapmaker.mapmakerutils.BinMap, map_id = mapid,
                         ts_map_key='DeflaggedTimestreams%s_%sGHz' %(direction, 
                                                                     band),
                         trans_key='RaDecTransform',
                         pointing_store_key='PixelPointing', 
                         timestream_weight_key = 'TodWeights')
    else:
        # Split by band (could do other things: wafer, bolos, etc.)
        pipe.Add(calibration.SplitTimestreamsByBand, 
                 input='DeflaggedTimestreams',
                 output_root='DeflaggedTimestreams')
        if verbose:
            pipe.Add(core.Dump)
        # Kick off maps
        for band in ['90', '150', '220']: # XXX should be automatic
            mapid = '%sGHz' % band
            pipe.Add(mapmaker.MapInjector, map_id = mapid,
                     maps_lst=[ps_map], is_stub=True, make_polarized=True, 
                     do_weight=True)
            pipe.Add(mapmaker.mapmakerutils.BinMap, map_id = mapid,
                     ts_map_key='DeflaggedTimestreams%sGHz' % band,
                     trans_key='RaDecTransform',
                     pointing_store_key='PixelPointing', 
                     timestream_weight_key = 'TodWeights')

    # Drop TOD    
    pipe.Add(lambda fr: fr.type == core.G3FrameType.Map)
    # Write to file
    pipe.Add(core.G3Writer, filename = output_file)
    pipe.Run(profile=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Maps for a CMB field')
                          # formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_files', nargs = '+')
    parser.add_argument('-o', '--output', default = 'output.g3',
                   help='Output filename')
    parser.add_argument('-r', '--res', action='store', default = 4.0, 
                   help='resolution [arcmin]')
    parser.add_argument('-x', '--xlen', default = 50, 
                   help='map width [deg]')
    parser.add_argument('-y', '--ylen', default = 30, 
                   help='map height [deg]')
    parser.add_argument('--psfile', default = None,
        help = 'Point source cofiguration file for making a point source mask')
    parser.add_argument('-v', '--verbose', action = 'store_true', 
                        default = False, help = 'Print every frame')
    args = parser.parse_args()

    args.res *= core.G3Units.arcmin
    args.xlen *= core.G3Units.deg
    args.ylen *= core.G3Units.deg
    args.psfile = '/home/ndhuang/spt_code/sptpol_software/config_files/ptsrc_config_ra0hdec-57p5_both_50mJy.txt'
    x_len = int(args.xlen / args.res)
    y_len = int(args.ylen / args.res)

    # Generate map stub
    map_params = coordinateutils.FlatSkyMap(
        x_len = x_len, y_len = y_len, res = args.res,
        alpha_center = 0., delta_center = -57.5 * core.G3Units.deg,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea)

    # Set up a map to feed parameters to the mapmaker
    ps_params = coordinateutils.FlatSkyMap(
        x_len = x_len, y_len = y_len, res = args.res,
        alpha_center = 0., delta_center = -57.5 * core.G3Units.deg,
        pol_type = core.MapPolType.T,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea)
    # Now, fill a map with the point source mask
    mapmaker.pointsourceutils.make_point_source_map(ps_params, args.psfile)

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=args.input_files)
    pipe.Add(field_mapmaker, map_params = map_params, ps_map = ps_params,
             output_file = args.output)
    pipe.Run()
