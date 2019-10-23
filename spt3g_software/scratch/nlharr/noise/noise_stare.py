from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk
from spt3g import todfilter, coordinateutils, timestreamflagging
from spt3g.util.extractdata import ExtractDataFromPipeline, ExtractCombinedTimestreams
from spt3g.calibration.template_groups import get_template_groups
from spt3g.util.genericutils import sort_two_lists, str_cmp_with_numbers_sorted

import pickle, argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--dir', help='noise stare directory', default='/spt/data/bolodata/downsampled/RCW38-pixelraster/32009821/')
parser.add_argument('--label', help='string that gets prepended to the saved pngs', default = 'tacostacostacos')
args = parser.parse_args()


save_label = args.label
noise_dir = args.dir 
input_files = [noise_dir + 'offline_calibration.g3',
               noise_dir + '0000.g3']

def print_stat(frame):
    if frame.type != core.G3FrameType.Scan:
        return
    ts = frame['CalTimestreams_I']
    ts = todfilter.polyutils.poly_filter_g3_timestream_map(ts, 8)
    for v in ts.values():
        print( (np.var(v)/(core.G3Units.uK**2.0) * v.sample_rate/core.G3Units.hz))**0.5
               

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = input_files)
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)
pipe.Add(core.Dump)


         

#verbose, but can easily be expanded to Q if you want
for tag in ['I']:
    pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
             ts_key = 'RawTimestreams_%s'%tag)

    pipe.Add(timestreamflagging.RemoveFlaggedTimestreams,
             input_ts_key = 'RawTimestreams_%s'%tag,
             output_ts_key = 'FlaggedRawTimestreams_%s'%tag,
             input_flag_key = 'Flags')

    pipe.Add(dfmux.ConvertTimestreamUnits, Input='FlaggedRawTimestreams_%s'%tag, 
             Output='TimestreamsWatts_%s'%tag, Units=core.G3TimestreamUnits.Power)

    pipe.Add(dfmux.ConvertTimestreamUnits, Input='FlaggedRawTimestreams_%s'%tag, 
             Output='TimestreamsAmps_%s'%tag, Units=core.G3TimestreamUnits.Current)

    pipe.Add(calibration.ApplyTCalibration, Source='RCW38', InKCMB=True,
             OpacityCorrection=True, Input='TimestreamsWatts_%s'%tag, 
             Output='CalTimestreams_%s' % tag)

    pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,
             variance_prefilter_poly_scale = 8,
             ts_key = 'CalTimestreams_%s'%tag)
    pipe.Add(timestreamflagging.RemoveFlaggedTimestreams,
             input_ts_key = 'CalTimestreams_%s'%tag,
             output_ts_key = 'FlaggedCalTimestreams_%s'%tag,
             input_flag_key='Flags')

    pipe.Add(timestreamflagging.GenerateFlagStats, flag_key='Flags')


#pipe.Add(print_stat)
pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
meta_data = ExtractDataFromPipeline(keys=['BolometerProperties', 'WiringMap', 
                                          'DfMuxHousekeeping'])
pipe.Add(meta_data)
ts_data = ExtractCombinedTimestreams(keys=[
        'FlaggedCalTimestreams_I',
        'CalTimestreams_I', #'CalTimestreams_Q', 
        'RawTimestreams_I', #'RawTimestreams_Q',
        'TimestreamsWatts_I',#, 'TimestreamsWatts_Q'
        'TimestreamsAmps_I'
                                           ])
#pipe.Add(ts_data)

pipe.Run()


def sort_frequency(hk_map, wiring_map, ts_map):
    ts_keys = ts_map.keys()
    freqs = []
    for k in ts_keys:
        bhk = dfmux.HousekeepingForBolo(hk_map, wiring_map, k)
        freqs.append( bhk.carrier_frequency )
    freqs, ts_keys = sort_two_lists(freqs, ts_keys)
    return ts_keys, freqs


def sort_board(wiring_map, bolo_props, ts_map):
    tgs = get_template_groups(bolo_props, wiring_map, per_band = False, 
                              per_board = True, per_squid = True, per_pixel = True,
                              include_keys = True)
    sorted_keys = []
    boards = []

    for k in sorted(tgs.keys()):
        board = k.split('_')[0]
        sq = k.split('_')[1]
        for ch_id in tgs[k]:
            sorted_keys.append(ch_id)
            boards.append(board)
    return sorted_keys, boards

def sort_band_board(wiring_map, bolo_props, ts_map):
    tgs = get_template_groups(bolo_props, wiring_map, per_band = True, 
              
                              per_board = True, per_squid = True, per_pixel = True,
                              include_keys = True)
    sorted_keys = []
    labels = []

    for k in sorted(tgs.keys(), str_cmp_with_numbers_sorted):
        band = k.split('_')[0]
        board = k.split('_')[1]
        sq = k.split('_')[2]
        for ch_id in tgs[k]:
            sorted_keys.append(ch_id)
            labels.append(band+'_'+board)
    return sorted_keys, labels


def convert_ts_map_to_grid(hk_map, wiring_map, bolo_props, ts_map, q_ts_map=None,
                           sort_rule = 'board'):
    #find the keys in some meaningful order
    if q_ts_map is None:
        sf = 1
    else:
        sf = 2
    ts_grid = np.zeros((len(ts_map.keys()) * sf,
                        len(ts_map[ts_map.keys()[0]]) ))


    if sort_rule == 'board':
        sorted_keys, labels = sort_board(wiring_map, bolo_props, ts_map)
    elif sort_rule == 'freq':
        sorted_keys, labels = sort_frequency(hk_map, wiring_map, ts_map)
    elif sort_rule == 'band_board':
        sorted_keys, labels = sort_band_board(wiring_map, bolo_props, ts_map)
    valid_keys = []
    valid_labels = []

    ind = 0
    for i, k in enumerate(sorted_keys):
        if k in ts_map:
            ts_grid[ind*sf,:] = ts_map[k]
            if not q_ts_map is None:
                ts_grid[ind*sf+1,:] = q_ts_map[k]
            ind += 1
            valid_keys.append(k)
            valid_labels.append(labels[i])
    return ts_grid, valid_labels


def normalize_ts_map_variance(ts_map):
    import copy
    out_tsm = copy.copy(ts_map)
    for k in out_tsm.keys():
        out_tsm[k] /= np.std(out_tsm[k])
    return out_tsm

i_data = todfilter.polyutils.poly_filter_g3_timestream_map(ts_data.CalTimestreams_I, 27)
#q_data = todfilter.polyutils.poly_filter_g3_timestream_map(ts_data.CalTimestreams_Q, 14)

#psd not used, but hey, it' here I guess
if 0:
    psd = todfilter.dftutils.get_psd_of_ts_map(i_data)
    for k in psd[0].keys():
        plt.plot((psd[1]/core.G3Units.hz)[2:], psd[0][k][2:], 'k.', alpha = 0.02)
    plt.show()

#correlation
if 0:
    #for sort_rule in ['board', 'freq', 'band_board']:
    for sort_rule in ['board']:
        i_data_norm = normalize_ts_map_variance(i_data)
        print("Gridding")
        grid, grid_labels  = convert_ts_map_to_grid(
            meta_data.DfMuxHousekeeping,
            meta_data.WiringMap, 
            meta_data.BolometerProperties, 
            i_data_norm, 
            sort_rule = sort_rule)

        plt.clf()
        if sort_rule == 'board' or sort_rule == 'band_board':
            print("Setting Labels")
            board_boundaries = [0]
            for i in range(1, len(grid_labels)):
                if grid_labels[i] != grid_labels[i-1]:
                    board_boundaries.append(i)
            board_boundaries.append(np.shape(grid)[0])


            tick_locations = []
            tick_labels = []
            for i in range(1, len(board_boundaries)):
                tl = (board_boundaries[i]+board_boundaries[i-1])/2.0
                if (board_boundaries[i] - board_boundaries[i-1] == 0):
                    continue
                tick_locations.append(tl)
                tick_labels.append( grid_labels[int(tl)] )
            plt.xticks(tick_locations, tick_labels,rotation=90, size=4)
            plt.yticks(tick_locations, tick_labels, size = 4)
        print("Correlation")
        corr = np.cov(grid)
        print("Plotting")


        #plt.imshow(np.log(np.abs(corr)), interpolation='none',cmap='Greys_r')
        plt.imshow(np.abs(corr), interpolation='none',cmap='Greys_r')
        plt.savefig("%s_noise_stare_%s.png"%(save_label, sort_rule), dpi=400)



