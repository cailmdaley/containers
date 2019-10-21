#!/usr/bin/env python
import unittest, h5py
import numpy as np
import cPickle as pickle
from spt3g.sptpol.directidfreader import DirectIdfReader
from spt3g.mapmaker import mapmakerutils as MMU
from spt3g import core, mapmaker, coordinateutils

def plot_real():
    test_data_path = mapmaker.get_test_files_path()
    idf_fn = test_data_path+'ra23h30dec-55_idf_20130430_003228_150ghz.h5'
    out_map = coordinateutils.FlatSkyMap(x_len=780, y_len=780,
                                  res=1.0 * core.G3Units.arcmin,
                                  delta_center=-55.0 * core.G3Units.deg,
                                  alpha_center=352.5 * core.G3Units.deg,
                                  proj=coordinateutils.MapProjection.Proj5)
    pipe = core.G3Pipeline()
    pipe.Add(DirectIdfReader, filename = idf_fn,  number_of_scans_to_read = 3)
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

def plot_sim():
    res = 0.25 * core.G3Units.arcmin
    theta_center = -55.0* core.G3Units.deg
    phi_center = 352.5* core.G3Units.deg
    proj =  coordinateutils.MapProjection.Proj5
    idf_fn = 'test_data_files/ra23h30dec-55_idf_20130430_003228_150ghz.h5'
    print('loading')
    out_map = coordinateutils.FlatSkyMap(x_len = 3400, y_len = 3400,
                                  res = res, theta_center = theta_center,
                                  phi_center = phi_center, proj = proj)
    
    sim_m_orig = pickle.load(open('test_data_files/OneTEBSimMap.pkl', 'rb'))
    sim_proj = coordinateutils.MapProjection.Proj5
    print 'into flatsky'
    tmp_map = np.ascontiguousarray(sim_m_orig[:,:,0])
    Tsim = coordinateutils.FlatSkyMap( tmp_map,
                                res = res, theta_center = theta_center,
                                phi_center = phi_center, proj = proj,
                                pol_type = coordinateutils.MapPolType.T
                            )
    print 'Q'
    tmp_map = np.ascontiguousarray(sim_m_orig[:,:,1])
    Qsim = coordinateutils.FlatSkyMap( tmp_map,
                                res = res, theta_center = theta_center,
                                phi_center = phi_center, proj = proj,
                                pol_type = coordinateutils.MapPolType.Q)
    print 'U'
    tmp_map = np.ascontiguousarray(sim_m_orig[:,:,2])
    Usim = coordinateutils.FlatSkyMap( tmp_map,
                                res = res, theta_center = theta_center,
                                phi_center = phi_center, proj = proj,
                                pol_type = coordinateutils.MapPolType.U)
    pipe = core.G3Pipeline()
    pipe.Add(DirectIdfReader, filename = idf_fn, load_bolo_data = False, 
             number_of_scans_to_read = 3)
    pipe.Add(MMU.FillSimTimestreams,
             T = Tsim, Q= Qsim, U = Usim,
             sim_map_id = 'test_sim_in')
    pipe.Add(MMU.MakeMap, map_in=out_map,
             map_id="test_sim_out", ts_in_key="TestOut",
             poly_order = -1, do_weight = True)
    pipe.Run(graph = True)
    core.plot_frame_processing_info(pipe)

if __name__ == '__main__':
    plot_real()
