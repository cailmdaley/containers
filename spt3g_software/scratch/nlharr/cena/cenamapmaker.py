from spt3g import core, mapmaker, sptpol, todfilter, dfmux, util, timestreamflagging
import spt3g.mapmaker.mapmakerutils as mmu
import spt3g.std_processing as std_processing
from copy import copy
import numpy as np
import pickle
#np.seterr(all='raise')

def get_point_source_mask(first_pass_map_fn, 
                          point_source_fn = '/home/nlharr/spt_code/spt3g_software/scratch/nlharr/cena/biggercenamask.pkl'):

    mp_frame = [fr for fr in core.G3File(first_pass_map_fn)][2]
    mp = mp_frame['T']

    alpha_center = mp_frame['BrightestPixelRaOrig']
    delta_center = mp_frame['BrightestPixelDecOrig']

    data = pickle.load(open(point_source_fn))
    return coordinateutils.FlatSkyMap(data, res = mp.res, proj = mp.proj, 
                               alpha_center = alpha_center, delta_center = delta_center,
                               pol_type = core.MapPolType.T)


def make_cena_map(cal_frame_fn, tod_fn, out_fn, is_first_pass = True, first_pass_map_fn = None):
    input_ts = 'RawTimestreams_I'
    kcmb_ts = 'KcmbTimestreams'
    
    if is_first_pass:
        individual_bolos_to_map = core.G3VectorString()
        pntsrc_mask = None
        use_dynamic_source_filter = True
    else:
        individual_bolos_to_map = std_processing.get_bolos_in_obs(tod_fn)
        pntsrc_mask = get_point_source_mask(first_pass_map_fn = first_pass_map_fn)
        use_dynamic_source_filter = False
    
    pipe = core.G3Pipeline()
    
    pipe.Add(core.G3Reader, filename = [cal_frame_fn, tod_fn])
    pipe.Add(mapmaker.RejectTurnArounds)
    pipe.Add(std_processing.PreprocessScanData, include_noise_info = True)
    

    pipe.Add(timestreamflagging.flaggingutils.GroupAverageG3MapValue,
         per_band = True,
         per_wafer = True,
         per_squid = True,
         input_g3_map_key = 'TimestreamMaxDerivative',
         output_g3_map_key = 'AvMaxTimestreamDerivative',
         average_func = np.median)
    
    pipe.Add(timestreamflagging.glitchfinding.FlagGlitches,
             min_num_glitch_map = {float(10):1},
         flag_key = 'Flags')
    
    pipe.Add(timestreamflagging.flaggingutils.FlagGroupAverageG3MapValue,
             g3_map_key = 'AvMaxTimestreamDerivative',
             min_val = 0,
             max_val = 3e-8,
             per_band = True,
             per_wafer = True,
         per_squid = True,
             flag_reason = 'BadMedianDeriv')
    
    pipe.Add(timestreamflagging.flaggingutils.FlagBadG3MapValue,
             flag_key = 'Flags',
             m_key = 'TimestreamVariance',
             flag_reason = 'BadRmsNoise',
             min_val = 1e-7,
             max_val = 2e-3)

    pipe.Add(timestreamflagging.flaggingutils.FlagBadG3MapValue,
             flag_key = 'Flags',
             m_key = 'TimestreamMaxDerivative',
             flag_reason = 'BadDeriv',
             min_val = 1e-11,
             max_val = 1e-7)
    
    pipe.Add(timestreamflagging.flaggingutils.RemoveFlaggedTimestreams,
             input_ts_key = kcmb_ts,
             input_flag_key = 'Flags',
             output_ts_key = 'FlaggedTs')

    map_in = std_processing.CreateSourceMapStub('cena', 200, 200, 0.25 * core.G3Units.arcmin, 
                                            proj = mapmaker.MapProjection.Proj5 )


    pipe.Add(mapmaker.MakeMap,
             map_in = map_in, 
             map_id = 'CenaMap',
             lpf_filter_frequency = 30 * core.G3Units.Hz,
             ts_in_key = 'FlaggedTs',
             
             do_weight = True,
             
             poly_order = 4,
             use_dynamic_source_filter = use_dynamic_source_filter,
             point_source_mask = pntsrc_mask,
             
             make_polarized = False,
             individual_bolos_to_map = individual_bolos_to_map,
             force_make_sum_map = True,
             
             use_boresight_pointing = False,
             fill_in_unity_weights = True,

             #point_source_mask_debug = True

         )
    
    pipe.Add(core.Dump)
    pipe.Add(mapmaker.mapmakerutils.RemoveWeightModule)

    pipe.Add(mapmaker.summingmaps.OffsetTMapFromPreviousOffset,
                 offset_map_id = 'CenaMap', desired_center_pixel = np.array([100,100]))
    
    
    pipe.Add(lambda frame: frame.type != core.G3FrameType.Scan)
    
    pipe.Add(core.G3Writer, filename = out_fn)

    pipe.Run(profile = False)


    

cena_line_up = {
    'cena-20150312_131911.g3':'rcw38-20150312_123711.g3',
    'cena-20150325_151459.g3':'rcw38-20150325_143300.g3',
    'cena-20150301_230758.g3':'rcw38-20150301_222600.g3',
    'cena-20150314_025927.g3':'rcw38-20150314_021727.g3',
    'cena-20150303_053730.g3':'rcw38-20150303_045532.g3',
    'cena-20150315_161543.g3':'rcw38-20150315_153345.g3',
    'cena-20150303_122107.g3':'rcw38-20150303_113907.g3',
    'cena-20150317_065622.g3':'rcw38-20150317_061426.g3',
    'cena-20150304_174541.g3':'rcw38-20150304_170344.g3',
    'cena-20150318_194938.g3':'rcw38-20150318_190454.g3',
    'cena-20150401_014312.g3':'rcw38-20150401_103312.g3',
    'cena-20150306_071532.g3':'rcw38-20150306_063336.g3',
    'cena-20150320_091251.g3':'rcw38-20150320_083055.g3',
    'cena-20150307_204739.g3':'rcw38-20150307_200253.g3',
    'cena-20150321_225618.g3':'rcw38-20150321_221418.g3',
    'cena-20150309_094036.g3':'rcw38-20150309_085837.g3',
    'cena-20150322_114432.g3':'rcw38-20150322_110232.g3',
    'cena-20150310_235454.g3':'rcw38-20150310_231255.g3',
    'cena-20150324_013009.g3':'rcw38-20150324_004811.g3',
}

tod_dir = '/data52/nwhitehorn/3g-calib/cena/tod/'
cal_dir = '/data/nlharr/cena_maps/cal_frames/'
out_dir_1 = '/data/nlharr/cena_maps/pass_1/'
out_dir_2 = '/data/nlharr/cena_maps/pass_2_bigger/'

is_first = False

if is_first:
    for k, v in cena_line_up.items():
        out_fn = out_dir_1 + '/map_'+k
        tod_fn = tod_dir+'/'+k
        cal_frame_fn = cal_dir +'/' + k
        make_cena_map(cal_frame_fn, tod_fn, out_fn)
else:
    for k, v in cena_line_up.items():
        out_fn = out_dir_2 + '/map_'+k
        first_pass_map_fn = out_dir_1 + '/map_'+k
        tod_fn = tod_dir+'/'+k
        cal_frame_fn = cal_dir +'/' + k
        make_cena_map(cal_frame_fn, tod_fn, out_fn, is_first_pass = False, 
                      first_pass_map_fn = first_pass_map_fn)
