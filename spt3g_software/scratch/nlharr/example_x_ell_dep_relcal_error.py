import numpy as np
import cPickle as pickle
from spt3g import core, mapmaker
from directidfreader import DirectIdfReader
from spt3g.mapmaker.mapmakerutils import MakeMap, ExtractTheMaps, FillSimTimestreams, remove_weight
import copy

idf_fn = mapmaker.get_test_files_path() + 'ra23h30dec-55_idf_20130430_003228_150ghz.h5'
sim_map_np = np.load(mapmaker.get_test_files_path() + 'silly_sim.npy')

@core.indexmod
def ElDepChange(frame, 
                input_ts_map_key, 
                output_ts_map_key,
                boresight_el_key = 'BoresightEl'):
    if frame.type != core.G3FrameType.Scan:
        return
    ts_map = copy.copy(frame[input_ts_map_key])
    el = np.mean(frame[boresight_el_key])/core.G3Units.deg
    resp_adjust = 1.0 + (el-55)/10.0
    for k in ts_map.keys():
        ts_map[k] *= resp_adjust
    frame[output_ts_map_key] = ts_map

out_map_parameters = coordinateutils.FlatSkyMap(x_len = 780, 
                                         y_len = 780,
                                         res = 1 * core.G3Units.arcmin, #units are important
                                         phi_center = 352.5 * core.G3Units.deg,
                                         theta_center = -55.0* core.G3Units.deg,
                                         proj= mapmaker.MapProjection.Proj5) #we use enums now
sim_maps = []
for i, ptype in enumerate( [core.MapPolType.T, core.MapPolType.Q, core.MapPolType.U]):
    sim_maps.append(coordinateutils.FlatSkyMap( sim_map_np[i,:,:],
                                         res = 1 * core.G3Units.arcmin, #units are important
                                         phi_center = 352.5 * core.G3Units.deg,
                                         theta_center = -55.0* core.G3Units.deg,
                                         proj= mapmaker.MapProjection.Proj5, 
                                         pol_type = ptype))

pipe = core.G3Pipeline()
pipe.Add(DirectIdfReader, filename = idf_fn, load_bolo_data = False)
pipe.Add(FillSimTimestreams,
         out_ts_key = 'SimulatedTimestreams',
         T = sim_maps[0], Q = sim_maps[1], U = sim_maps[2],  
         sim_map_id = 'test_sim_in')

pipe.Add(ElDepChange, input_ts_map_key = 'SimulatedTimestreams',
         output_ts_map_key = 'ElDepTimstreams')

pipe.Add(MakeMap, map_in = out_map_parameters,
         map_id="test_sim_out", ts_in_key='ElDepTimstreams',
         do_weight = True)
map_extractor = ExtractTheMaps()
pipe.Add(map_extractor)
pipe.Add(core.Dump)
pipe.Run()

T = map_extractor.maps["test_sim_out"]["T"]
Q = map_extractor.maps["test_sim_out"]["Q"]
U = map_extractor.maps["test_sim_out"]["U"]
W = map_extractor.maps["test_sim_out"]["Wpol"]
Tnw, Qnw, Unw = remove_weight(T, Q, U, W)
