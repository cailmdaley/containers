import os, time
import argparse as ap
import numpy as np
from glob import glob

from spt3g import core, dfmux, calibration, coordinateutils, gcp, mapmaker, std_processing, todfilter, timestreamflagging

# Usage: get_obs_psd.py <input files.g3> -o outputmaps.g3
P = ap.ArgumentParser(description='Timestream PSDS',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
           help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
P.add_argument('--polyfilt',action='store_true')
args = P.parse_args()

'''
Calculate the by-scan median psds of timestream key 'ts_key' for an observation.
Saves the result.

Update July 2018: Exclude w201 and its nasty lines.
'''

@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def FlagWafer(frame, wafer, flag_key = 'Flags', flag_reason = 'BadWafer',
               bolo_props=None):
    '''
    Flags any detectors on a specified wafer (e.g. w201)
    '''

    assert(not bolo_props is None)
    bad_bolos = []
    for k in bolo_props.keys():
        if bolo_props[k].wafer_id == wafer:
            bad_bolos.append(k)
    timestreamflagging.add_flag(frame, flag_key, flag_reason, bad_bolos)

def get_median_psd(frame, input):
    if input in frame:
        psd, freqs = todfilter.dftutils.get_psd_of_ts_map(frame[input])
        med_psd = np.zeros(len(freqs))
        num_bolos = len(psd.keys())
        resp_matrix = np.zeros((num_bolos,len(freqs)))
        for i, bolo in enumerate(psd.keys()):
            resp_matrix[i, :] = psd[bolo]
        for j in range(np.shape(resp_matrix)[1]):
            med_psd[j] = np.median(resp_matrix[:,j])
        
        frame['PSD'] = core.G3VectorDouble(med_psd)
        frame['PSD_Freqs'] = core.G3VectorDouble(freqs)
    else:
        return

def get_info(fr):
    if 'TrackerStatus' in fr:
        if np.gradient(fr['TrackerStatus'].az_pos)[0]>0:
            fr['Direction'] = 'right'
        else:
            fr['Direction'] = 'left'
        fr['Time'] = fr['TrackerStatus'].time[0]
    
def clean_frame(frame):
    extra = []
    for k in frame.keys():
        if k not in extra and \
        k not in ['PSD','PSD_Freqs','Direction','Time']:
            extra += [k]
    core.Delete(frame, keys=extra)
    
processing_start = time.time()

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = args.input_files)
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
pipe.Add(get_info)
pipe.Add(core.Delete,keys='RawTimestreams_Q')

# Most basic flags
# These make up std_processing.flagsegments.FlagInvalidData
# pipe.Add(timestreamflagging.FlagNaNs,ts_key='RawTimestreams_I',flag_key='Flags')
# pipe.Add(std_processing.flagsegments.FlagBadHousekeeping, ts_key='RawTimestreams_I', 
#          flag_key='Flags')

pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion,
         flag_key='Flags', ts_key = 'RawTimestreams_I')

# Convert to watts
pipe.Add(std_processing.CalibrateRawTimestreams, output='TimestreamsWatts',
        units = core.G3TimestreamUnits.Power)

pipe.Add(std_processing.flagsegments.FieldFlaggingPostKcmbConversion,
         flag_key='Flags', ts_key = 'TimestreamsWatts')

if args.polyfilt == True:
    ts = 'PolyFilteredTimestreams'
    pipe.Add(mapmaker.TodFiltering, ts_in_key='TimestreamsWatts',
            ts_out_key=ts, poly_order=9)
else:
    ts = 'TimestreamsWatts'

pipe.Add(std_processing.weighting.AddPSDWeights,
         input = ts, output = 'TodWeights',
         low_f = .5 * core.G3Units.Hz,
         high_f = 8 * core.G3Units.Hz)
pipe.Add(calibration.SplitByBand, 
         input = 'TodWeights', output_root = 'TodWeights')
for band in ['90', '150', '220']:
    pipe.Add(timestreamflagging.flaggingutils.SigmaclipFlagGroupG3MapValue,
             m_key = 'TodWeights%sGHz'%band, low = 2.5, high = 2.5, 
             flag_reason = 'BadWeight', flag_key = 'Flags')
pipe.Add(FlagWafer, wafer = 'w201', flag_key = 'Flags')
pipe.Add(timestreamflagging.GenerateFlagStats, flag_key = 'Flags')
# Drop flagged timestreams
pipe.Add(timestreamflagging.RemoveFlaggedTimestreams,
         input_ts_key = ts, input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')

pipe.Add(get_median_psd, input = 'DeflaggedTimestreams')
pipe.Add(lambda fr: fr.type == core.G3FrameType.Scan)
pipe.Add(clean_frame)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename = args.output)
pipe.Run(profile=True)

print('Took %f minutes'%((time.time()-processing_start)/60.))
