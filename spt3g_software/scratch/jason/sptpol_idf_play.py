import spt3g
import time, numpy

from spt3g import core, mapmaker, sptpol
from spt3g.mapmaker import mapmakerutils as MMU

idf_fn = '/data55/sptdat/idf/ra0hdec-57p5_V3/data/ra0hdec-57.5_idf_20130803_012652_150ghz.h5'

reader = sptpol.DirectIdfReader(filename=idf_fn)

# First calculate the min and max of ra and dec to have an idea of how to make a map
ramin,  ramax = (numpy.inf, -numpy.inf)
decmin,decmax = (numpy.inf, -numpy.inf)
frames = True
while frames:
    frames = reader(None)
    print(len(frames))
    for frame in frames:
        print frame
        if frame.type == core.G3FrameType.Scan:
            ramax = max(ramax, frame['BoresightRa'].max() )
            ramin = min(ramin, frame['BoresightRa'].min() )
            decmax = max(decmax, frame['BoresightDec'].max() )
            decmin = min(decmin, frame['BoresightDec'].min() )
            print(ramin/core.G3Units.deg, ramax/core.G3Units.deg, decmin/core.G3Units.deg, decmax/core.G3Units.deg)
            # (-31.401818827049265, 2.040692565518107, -65.439367951462472, -49.293652342628491)


# Now set up map-making parameters
resarcmin = 5.0
x_len=int(numpy.ceil(62*60/resarcmin*numpy.cos(numpy.deg2rad(-57.5)))) # 400
y_len=int(numpy.ceil(20*60/resarcmin)) # 240

out_map = coordinateutils.FlatSkyMap(x_len=x_len, y_len=y_len,
                              res=resarcmin * core.G3Units.arcmin,
                              delta_center=-57.5 * core.G3Units.deg,
                              alpha_center=0.0 * core.G3Units.deg,
                              proj=mapmaker.MapProjection.Proj5)
                              
# Now setup the pipeline to actually make the map.
pipe = core.G3Pipeline()
pipe.Add(sptpol.DirectIdfReader, filename = idf_fn)
pipe.Add(mapmaker.FlaggedRemover,
         in_ts_map_key="CalTimestreams",  out_ts_map_key = "FlaggedRemover",
         flag_key = "Flags", scan_flag_key = "ScanIsBad" )
pipe.Add(MMU.MakeMap, map_in=out_map, 
         map_id="test_out", ts_in_key="FlaggedRemover", 
         poly_order = 1,
         do_weight = True, keep_all_data = True)
pipe.Add(core.Dump)
pipe.Run(graph = True)
core.plot_frame_processing_info(pipe)

