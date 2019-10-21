'''
Script to make CMB field maps.
'''
import argparse
import sys
from spt3g import core, std_processing, mapmaker, calibration
from spt3g import timestreamflagging, todfilter, coordinateutils
from spt3g import std_processing
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations, AddAzElTrans


# Usage: makemaps.py <input files.g3> -o outputmaps.g3

@core.pipesegment
def field_mapmaker(pipe, split_leftright = False, map_params = None, 
                   output_file = None, verbose = False):
    assert(map_params is not None)

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
    #pipe.Add(FillCoordTransRotations, transform_store_key = "RaDecTransform", 
    #         do_bad_transform = True)

    # Calculate detector pointing early and add point source mask
    pipe.Add(mapmaker.MapInjector, map_id = 'UselessMapInfo',
             maps_lst = [map_params,], is_stub=True)



    pipe.Add(AddAzElTrans, az_key = 'OnlineBoresightAz', el_key = 'OnlineBoresightEl', 
             out_key = 'AzElTransform')

    pipe.Add(mapmaker.mapmakerutils.CalculatePointing, 
             map_id = 'UselessMapInfo', 
             pointing_store_key = 'PixelPointing', trans_key='AzElTransform',
             ts_map_key = 'CalTimestreams')

    # Basic timestream filtering
    pipe.Add(mapmaker.TodFiltering,
             # filtering options
             poly_order = 9, 
             filters_are_ell_based = True, 
             mhpf_cutoff = 100, lpf_filter_frequency = 6600,
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
    pipe.Add(calibration.SplitTimestreamsByBand, 
             input='DeflaggedTimestreams',
             output_root='DeflaggedTimestreams')
    if verbose:
        pipe.Add(core.Dump)
    # Kick off maps
    for band in ['90', '150', '220']: # XXX should be automatic
        mapid = '%sGHz' % band
        pipe.Add(mapmaker.MapInjector, map_id = mapid,
                 maps_lst=[map_params], is_stub=True, make_polarized=False, 
                 do_weight=True)
        pipe.Add(mapmaker.mapmakerutils.BinMap, map_id = mapid,
                 ts_map_key='DeflaggedTimestreams%sGHz' % band,
                 trans_key='AzElTransform',
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
    parser.add_argument('-i', '--input_files', nargs = '+', 
                        default = ['/spt/data/bolodata/downsampled/ra0hdec-57.5/35886130/offline_calibration.g3',
                                   '/spt/data/bolodata/downsampled/ra0hdec-57.5/35886130/0000.g3'])
    parser.add_argument('-o', '--output', default = 'output.g3',
                   help='Output filename')
    parser.add_argument('-v', '--verbose', action = 'store_true', 
                        default = False, help = 'Print every frame')
    args = parser.parse_args()

    map_params = std_processing.sourcemaps.CreateGroundMap(
        res = 0.25 * core.G3Units.deg, delta_center = 50.0 * core.G3Units.deg,
        delta_angular_extent = 30 * core.G3Units.deg)

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=args.input_files)
    pipe.Add(field_mapmaker, map_params = map_params, output_file = args.output, 
             verbose=args.verbose)
    pipe.Run()
