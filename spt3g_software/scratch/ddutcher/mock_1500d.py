'''
mock_1500d.py

Script to mock-observe simulated skies
with same processing as 1500d

Utilizes ddutcher's SimStubs that live with the maps
'''
import argparse
from spt3g import core, std_processing, mapmaker, calibration
from spt3g import timestreamflagging, todfilter, coordinateutils
from spt3g.std_processing.pointing import CalculateCoordTransRotations

# Suppress all warnings about timestream missing from scan frame
core.set_log_level(core.G3LogLevel.LOG_ERROR, 'MapBinner')

parser = argparse.ArgumentParser(description='Mock-observe CMB fields')
parser.add_argument('input_files',
                   help = 'The SimStub containg real scan pointing/flags/weights')
parser.add_argument('-m', '--mapsim', default ='None',
                    help = 'The simulated sky to mock obseve')
parser.add_argument('-o', '--output', default = 'output.g3',
                   help='Output filename')
parser.add_argument('-s', '--source', default = 'ra0hdec-57.5',
                   help='Name of source, to set field center')
parser.add_argument('-r', '--res', default = 2.0, type=float,
                   help='resolution [arcmin]')
parser.add_argument('-x', '--xlen', default = 75, type=float,
                   help='map width [deg]')
parser.add_argument('-y', '--ylen', default = 50, type=float,
                   help='map height [deg]')
parser.add_argument('--psfile', default = None,
                   help = 'Point source configuration file for making a point source mask')
parser.add_argument('-v', '--verbose', action = 'store_true', default=False,
                    help = 'Print every frame')
parser.add_argument('--lr', action = 'store_true',default=False,
                    help = 'Split left-right')
parser.add_argument('--tonly', default = True, action = 'store_false',
                    help = 'Include this flag to make T-only maps')
args = parser.parse_args()

args.res *= core.G3Units.arcmin
args.xlen *= core.G3Units.deg
args.ylen *= core.G3Units.deg
x_len = int(args.xlen / args.res)
y_len = int(args.ylen / args.res)

# Generate map stub
map_params = std_processing.CreateSourceMapStub(
    args.source, x_len = x_len, y_len = y_len, res = args.res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea)

if args.psfile is not None:
    # Set up an empty map for point source filtering
    ps_params = std_processing.CreateSourceMapStub(
        args.source, x_len = x_len, y_len = y_len, res = args.res,
        proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea)
    # Now, fill map with the point source mask
    mapmaker.pointsourceutils.make_point_source_map(ps_params, args.psfile)
    
    ps_map_id = 'PointSourceMask'
else:
    ps_map_id = None

# Begin pipeline
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = args.input_files)

###
# Make fake timestreams
###

# Load the simulated map 
m = mapmaker.mapmakerutils.load_spt3g_map(args.mapsim)
pipe.Add(mapmaker.MapInjector, map_id='sim_map',maps_lst=[m['T'],m['Q'],m['U']], is_stub=False)

# Calculate the sim map pointing to give the sim TODs
# something to look at
def foo(fr, bolo_id_key = 'valid_ids', output = 'valid_ids_tsm'):
    # CalculatePointing wants a TimestreamMap. Fine.
    if bolo_id_key in fr:
        tmp = core.G3TimestreamMap()
        for b in fr[bolo_id_key]:
            tmp[b]=core.G3Timestream()
        fr[output] = tmp
pipe.Add(foo)

# interp_sim is slooooooow
interp_sim = False
if not interp_sim:
    pipe.Add(mapmaker.CalculatePointing,map_id = 'sim_map',
            pointing_store_key = 'SimPixelPointing',
            ts_map_key = 'valid_ids_tsm',
            trans_key='OnlineRaDecRotation')
    
pipe.Add(mapmaker.mapmakerutils.FillSimTodSegment,
        out_ts_key = 'SimTimestreams',
        ts_to_get_sample_rate_key = 'OnlineBoresightRa',
        interp_sim = interp_sim,
        map_is_healpix = True,
        trans_key='OnlineRaDecRotation',
        sim_map_id = 'sim_map',
        valid_ids = 'valid_ids',
        sim_pointing_key = 'SimPixelPointing' if not interp_sim else '')

###
# Repeat processing on simulated timestreams
###

# Add point source mask
if args.psfile is not None:
    pipe.Add(mapmaker.MapInjector, map_id =  ps_map_id,
             maps_lst = [ps_params,], is_stub=False)  

# Calculate the output map pointing so that our
# sim timestreams get binned into something useful 

pipe.Add(mapmaker.MapInjector, map_id = 'bsmap',
         maps_lst = [map_params,], is_stub=False)
pipe.Add(mapmaker.CalculatePointing,map_id = 'bsmap',
        pointing_store_key = 'PixelPointing',
        ts_map_key = 'SimTimestreams',
        trans_key='OnlineRaDecRotation')

pipe.Add(todfilter.polyutils.CommonModeFilter(
    in_ts_map_key = 'SimTimestreams', 
    out_ts_map_key = 'CMFilteredSimTimestreams',
    per_band = False, per_wafer = True, per_squid = False)) 
# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering,
         # filtering options
         poly_order = 19, 
         filters_are_ell_based = True, 
         mhpf_cutoff = 300, lpf_filter_frequency = 6600,
         point_source_mask_id = ps_map_id,
         # boiler plate
         ts_in_key='CMFilteredSimTimestreams',
         ts_out_key = 'PolyFilteredSimTimestreams', 
         point_source_pointing_store_key = 'PixelPointing',
         use_dynamic_source_filter = False,
         boresight_az_key='OnlineBoresightAz',
         boresight_el_key='OnlineBoresightEl')
# Remove leftover flagged bolos
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredSimTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')

# Clean up junk
pipe.Add(core.Delete,keys = ['SimTimestreams', 'CMFilteredSimTimestreams',
                             'PolyFilteredSimTimestreams'])
    
###    
# Make maps with sim timestreams
###

# Kick off maps
if args.lr:
    # split left and right scans if requested
    pipe.Add(std_processing.pointing.split_left_right_scans)
    # Split by band (could do other things: wafer, bolos, etc.)
    for direction in ['Left', 'Right']:
        pipe.Add(calibration.SplitByBand,
                 input='DeflaggedTimestreams' + direction,
                 output_root='DeflaggedTimestreams' + direction)
    if args.verbose:
        pipe.Add(core.Dump)
    for direction in ['Left', 'Right']:
        for band in ['90', '150', '220']: # XXX should be automatic
            mapid = '%s-%sGHz' %(direction, band)
            pipe.Add(mapmaker.MapInjector, map_id = mapid,
                     maps_lst=[map_params], is_stub=True, 
                     make_polarized=args.tonly, do_weight=True)
            pipe.Add(mapmaker.mapmakerutils.BinMap, map_id = mapid,
                     ts_map_key='DeflaggedTimestreams%s%sGHz' %(direction, 
                                                                 band),
                     trans_key='OnlineRaDecRotation',
                     pointing_store_key='PixelPointing', 
                     timestream_weight_key = 'TodWeights')
else:
    # Split by band (could do other things: wafer, bolos, etc.)
    pipe.Add(calibration.SplitByBand,
             input='DeflaggedTimestreams',
             output_root='DeflaggedTimestreams')
    if args.verbose:
        pipe.Add(core.Dump)
    # Kick off maps
    for band in ['90', '150', '220']: # XXX should be automatic
        mapid = '%sGHz' % band
        pipe.Add(mapmaker.MapInjector, map_id = mapid,
                 maps_lst=[map_params], is_stub=True, 
                 make_polarized=args.tonly, do_weight=True)
        pipe.Add(mapmaker.mapmakerutils.BinMap, map_id = mapid,
                 ts_map_key='DeflaggedTimestreams%sGHz' % band,
                 trans_key='OnlineRaDecRotation',
                 pointing_store_key='PixelPointing', 
                 timestream_weight_key = 'TodWeights')

# Drop TOD    
pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan)

# Write to file
pipe.Add(core.G3Writer, filename = args.output)
pipe.Run(profile=True)