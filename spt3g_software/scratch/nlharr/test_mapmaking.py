from spt3g import mapmaker, todfilter
from spt3g import core
from spt3g.mapmaker import mapmakerutils as MMU

from directidfreader import DirectIdfReader

import numpy as np


#######################
#make the input sim map
#######################
#core.G3Logger.global_logger.set_level(core.G3LogLevel.LOG_DEBUG)

#################################
#make the shell o the output map #
#################################
nx = 472
ny = 200
proj = mapmaker.MapProjection.Proj1
phi_center = 0 * core.G3Units.deg
theta_center = -57.5 * core.G3Units.deg
res = 7.5 * core.G3Units.arcmin
is_polarized = True

out_map = mapmaker.FlatMapMat1x3(n_x = 472, n_y = 200,
                                 is_weight = False, is_polarized = True, initialize_map = False, res = res,
                                 theta_center = -57.5 * core.G3Units.deg, 
                                 phi_center = 7.5 * core.G3Units.arcmin,
                                 proj=mapmaker.MapProjection.Proj1)

 #############################
# Make the point source maps  #
 #############################
from glob import glob
#idf_fns = ['/home/nlharr/spt_code/spt3g_software_AdventuresInForcePushing/scratch/sptpol_idfs/ra23h30dec-55_idf_20130430_003228_150ghz.h5']
import usefulmodules
outm = None
outw = None
map_extractor = MMU.ExtractTheMaps()

#idf_fn = '/home/nlharr/spt_code/spt3g_software/scratch/test_detritus/sptpol_idfs/ra23h30dec-55_idf_20130430_003228_150ghz.h5'
idf_fn = '/data/sptdat/idf/ra23h30dec-55_bb/data/ra23h30dec-55_idf_20130430_003228_150ghz.h5'
pipe = core.G3Pipeline()
pipe.Add(DirectIdfReader, filename = idf_fn, preload_data = True)
pipe.Add(usefulmodules.print_frame)
pipe.Add(mapmaker.FlaggedRemover,
         in_ts_map_key="CalTimestreams",  out_ts_map_key = "FlaggedRemover",
         flag_key = "Flags", scan_flag_key = "ScanIsBad" )
pipe.Add(MMU.make_map_pipe_segment, map_in=out_map, 
         poly_order = 4, mhpf_cutoff = 0.1, lpf_filter_frequency = 30,
         map_id="test_out", ts_in_key="FlaggedRemover")
pipe.Run(profile=True)
    
