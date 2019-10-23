from spt3g import core, todfilter, timestreamflagging, mapmaker, std_processing, calibration, dfmux, gcp, coordinateutils
from glob import glob
import numpy as np
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations

root_folder = '/spt/data/bolodata/downsampled/saturn/34344067/'
cal_file = '/spt/user/production/calibration/calframe/saturn/34344067.g3'
input_file_lst = [cal_file] + sorted(glob(root_folder + '0*.g3'))

input_ts = 'RawTimestreams_I'

map_id = 'SaturnMap'

start_time = None
bolos = None
if not bolos:
    for fname in input_file_lst:
        if (bolos is None) or (start_time is None):
            for frame in core.G3File(fname):
                if 'CalibratorResponseSN' in frame:
                    bolos = [k for k in frame['CalibratorResponseSN'].keys() if frame['CalibratorResponseSN'][k] > 20]
                if input_ts in frame:
                    start_time = frame[input_ts].start
                    break
        
assert(not start_time is None)

# creates the flat sky map we are using to fill in information
# this function returns a coordinateutils.FlatSkyMap with the center coordinates at a source
# you can also call the constructor directly
map_in = std_processing.CreateSourceMapStub(
    'saturn', 600, 600,  #map dimensions
    0.5 * core.G3Units.arcmin, #resolution
    #the map projection is stored in an enum
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea, #map projection
    at_time = start_time #because planets move
    )

# start our mapmaking script
pipe = core.G3Pipeline()

#read our files, the calibration file first, then the TOD files
pipe.Add(core.G3Reader, filename = input_file_lst)

#prevents multiple identitical meta data frames going through
pipe.Add(core.DeduplicateMetadata)

#converts our timestream units to watts
pipe.Add(std_processing.CalibrateRawTimestreams,
         units=core.G3TimestreamUnits.Power,
         output = 'BoloMapTimestreams')

#rejects scans that are the turn arounds of the telescope
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

#remove timestreams that are not in the bolometer map
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, 
         input='BoloMapTimestreams', output='BoloMapTimestreamsFilt')

# Use the timestream derivative value estimated above to flag bolometers
pipe.Add(timestreamflagging.flaggingutils.FlagBadG3MapValueCal,
         flag_key = 'Flags', #where we store the flags
         m_key = 'CalibratorResponseSN', # the key that has the G3MapDouble with values we are flagging on
         flag_reason = 'BadCalResponse', #the string we label the flagged detectors.
         min_val = 20, #values less than this are flagged 
         max_val = 1e30 #values greater than this are flagged
         )

#actually remove the timestreams
pipe.Add(timestreamflagging.flaggingutils.RemoveFlaggedTimestreams,
         input_ts_key = 'BoloMapTimestreamsFilt', 
         input_flag_key = 'Flags', 
         output_ts_key = 'FlaggedTs' 
         )

#this module is used to add a map frame to pipeline.
# it takes the G3SkyMap and puts it in the appropriate frame format
# and then pushes it on the pipe.
# All subsequent modules that need map information get it from the frame in the pipeline
# The map_id is unique string id that is assigned to the map frame that's used to identify
# it to the modules down the pipe
pipe.Add(mapmaker.MapInjector,
         map_id = map_id,
         maps_lst = [map_in],
         is_stub = True,
         make_polarized = False,
         do_weight = True)

pipe.Add(core.Dump)

#figures out the individual detector pointing relative to boresight
#and which pixel in the map they are attached to
pipe.Add(FillCoordTransRotations, transform_store_key = "RaDecTransform", 
         do_bad_transform = True)

pipe.Add(mapmaker.CalculatePointing,
         map_id = map_id,
         trans_key = 'RaDecTransform',
         pointing_store_key = 'SaturnPixelPointing',
         ts_map_key = 'FlaggedTs' )

#applies a 4th order polynomial subtraction to the data
#while excluding really large excursions when using the filter
pipe.Add(mapmaker.TodFiltering,
         ts_in_key = 'FlaggedTs',
         ts_out_key = 'FilteredTs',
         use_dynamic_source_filter = True,
         poly_order = 4)

#calculate weights
pipe.Add(std_processing.weighting.AddPSDWeights,
        input = 'FilteredTs',
        output = 'TimestreamWeights')

#creates our weighted histogram
pipe.Add(mapmaker.BinMap,
         map_id = map_id, 
         ts_map_key = 'FilteredTs',
         pointing_store_key = 'SaturnPixelPointing',
         timestream_weight_key = 'TimestreamWeights'
         #individual_bolos_to_map = bolos,
         )

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename = 'saturn_output_individual_test.g3')

pipe.Run(profile = True)
