from spt3g import core, mapmaker, sptpol, todfilter, dfmux, util
import spt3g.mapmaker.mapmakerutils as mmu
import spt3g.std_processing as std_processing
from spt3g.calibration.apply_t_cal import ApplyTCalibration
from spt3g import coordinateutils, calibration
from copy import copy
import numpy as np

import pickle

def call_the_func(cal_frame_fn, tod_fn, out_fn):
    pntsrc_mask = pickle.load(open('point_source_mask.pkl'))
    input_ts = 'RawTimestreams_I'
    kcmb_ts = 'BoloMapTimestreams'
    
    individual_bolos_to_map = std_processing.get_bolos_in_obs(tod_fn)
    
    pipe = core.G3Pipeline()

    precache = True
    if precache:
        pipe.Add(core.G3Reader, filename = [cal_frame_fn, tod_fn])
        
        #pipe.Add(core.InjectDebug, type = core.G3FrameType.Calibration)
        
        pipe.Add(mapmaker.RejectTurnArounds)
        
        pipe.Add(dfmux.ConvertTimestreamUnits, Input=input_ts, 
                 Output='BoloMapTimestreamsWatts', 
                 Units=core.G3TimestreamUnits.Power)
        
        pipe.Add(ApplyTCalibration, 
                 Input='BoloMapTimestreamsWatts', 
                 Output = kcmb_ts , SkyCal='RCW38FluxCalibration')
        
        pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
                 az_timestream='RawBoresightAz', el_timestream='RawBoresightEl',
                 ra_timestream='RawBoresightRa', dec_timestream='RawBoresightDec')         
        
        pipe.Add(todfilter.dftutils.AddNoiseInfo,
                 ts_map_key = kcmb_ts,
                 plot_psds = False,
                 save_folder = '/home/nlharr/spt3g_software/scratch/nlharr/rcw38map/plots/',
                 save_tag = 'psdvals'
             )
        
        pipe.Add(todfilter.flaggingutils.AddNumGlitches,
                 thresholds = [10, 15, 20],
                 input_ts_key = kcmb_ts)
        

        pipe.Add(todfilter.flaggingutils.GroupAverageG3MapValue,
                 per_band = True,
                 per_wafer = True,
                 per_squid = True,
                 input_g3_map_key = 'MaxTimestreamDerivative',
                 output_g3_map_key = 'BandMedianAvMaxTimestreamDerivative',
                 average_func = np.median)

        #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
        

    if precache:
        pipe.Add(core.G3Writer, filename = 'cachedtod.g3')
    else:
        pipe.Add(core.G3Reader, filename = 'cachedtod.g3')

    #pipe.Add(core.Dump)
    #pipe.Add(mapmaker.AddTodNoiseWeights, psd_key = 'TimestreamBinnedPsd_20.00_30.00Hz')
        #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
        
    pipe.Add(todfilter.flaggingutils.FlagGlitches,
             thresholds = [10], num_glitches = [1],
             input_ts_key = kcmb_ts,
             flag_key = 'Flags')


    pipe.Add(todfilter.flaggingutils.FlagGroupAverageG3MapValue,
             g3_map_key = 'BandMedianAvMaxTimestreamDerivative',
             min_val = 0,
             max_val = 3e-8,
             per_band = True,
             per_wafer = True,
             per_squid = True,
             flag_reason = 'BadMedianDeriv')


    pipe.Add(todfilter.flaggingutils.FlagBadG3MapValue,
             flag_key = 'Flags',
             m_key = 'TimestreamVariance',
             tag = 'BadRmsNoise',
             min_val = 1e-7,
             max_val = 2e-3)


    pipe.Add(todfilter.flaggingutils.FlagBadG3MapValue,
             flag_key = 'Flags',
             m_key = 'MaxTimestreamDerivative',
             tag = 'BadDeriv',
             min_val = 1e-11,
             max_val = 1e-7)



    pipe.Add(sptpol.LoadPolAngsFromIDF, 
             fn = '/home/nlharr/tmp/ra0hdec-57.5_idf_20150520_103653_150ghz.h5')


    pipe.Add(todfilter.flaggingutils.FlagToBand, band = 1)

    pipe.Add(todfilter.flaggingutils.RemoveFlaggedTimestreams,
             input_ts_key = kcmb_ts,
             input_flag_key = 'Flags',
             output_ts_key = 'FlaggedTs')

    #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)


    map_in = std_processing.CreateSourceMapStub('cena', 100, 100, 0.25 * core.G3Units.arcmin, 
                                 proj = mapmaker.MapProjection.Proj5 )

    

    #pipe.Add(core.InjectDebug, type = core.G3FrameType.Calibration)
    pipe.Add(mapmaker.MakeMap,
             map_in = map_in, 
             map_id = 'dummy',
             lpf_filter_frequency = 30 * core.G3Units.Hz,
             ts_in_key = 'FlaggedTs',

             do_weight = True,
            
             poly_order = 4,
             use_dynamic_source_filter = False,
             point_source_mask = pntsrc_mask,
             point_source_mask_id = 'pntsrc_mask',

             make_polarized = True,
             #individual_bolos_to_map = individual_bolos_to_map,
           
             use_boresight_pointing = False,
             fill_in_unity_weights = True
    )


    pipe.Add(lambda frame: frame.type != core.G3FrameType.Scan)

    pipe.Add(core.G3Writer, filename = out_fn)

    pipe.Run(profile = True)
    

#fn = 'cena-20150404_031558.g3'
fn = 'cena-20150405_160706.g3'
cal_frame_fn = '/home/nlharr/spt3g_software/scratch/nlharr/rcw38map/tmpdata/cal/' + fn
tod_fn = '/home/nlharr/spt3g_software/scratch/nlharr/rcw38map/tmpdata/tod/' + fn

out_fn= 'cena_bolo_maps_3.g3'
#out_fn= 'cena_sum_map_2.g3'

#out_fn= 'prof.g3'
call_the_func(cal_frame_fn, tod_fn, out_fn)
