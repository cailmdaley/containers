'''
Script to make individual bolometer maps with boresight pointing for use in
pointing calibration (e.g. for RCW38).
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk
from spt3g.mapmaker.mapmakerutils import MakeMap
import os

# Usage: makebolomaps.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
		      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
	       help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
	       help='Output filename')
P.add_argument('--input-ts', action='store',
	       choices=['RawTimestreams_I', 'RawTimestreams_Q'],
	       default='RawTimestreams_I', help='Input timestream type')
P.add_argument('-s', '--source', action='store', 
	       default='rcw38', help='name of source')
P.add_argument('-c', '--calpath', action='store', default=None,
               help='path to directory of calibrator results')
P.add_argument('-x', '--xtalkpath', action='store', default=None,
               help='path to directory of crosstalk results')
P.add_argument('-r', '--res', action='store', 
	       default=0.5, help='resolution [arcmin]')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# get start time of first file
fname = args.input_files[1]
for frame in core.G3File(fname):
    if args.input_ts in frame:
        starttime = frame[args.input_ts].start
        break

print('starttime is: ')
print starttime
print('\n')

# Guess the list of bolos to use
bolos = None
if not bolos:
    bolos = []
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if frame.type == core.G3FrameType.Calibration:
                keys = frame['CalibratorResponseSN'].keys()
                for key in keys:
                    if frame['CalibratorResponseSN'][key] > 5. and frame['ElnodSNSlopes'][key] < -10.:
                        bolos.append(key)

print len(bolos)
print bolos[0:10]

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = 3.*core.G3Units.deg/res,
    y_len = 3.*core.G3Units.deg/res, res = res,
    proj = mapmaker.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

input = []

if args.calpath is not None:
    calib_data = calibration.build_cal_frames.CalResultsForTime(starttime, {'calib': os.path.join(args.calpath, 'calib-%(fileformatstr)s.g3')})
    input.append(calib_data['calib'])

if args.xtalkpath is not None:
    xtalk_data = calibration.build_cal_frames.CalResultsForTime(starttime, {'xtalk': os.path.join(args.xtalkpath, 'xtalk-%(fileformatstr)s.g3')})
    input.append(xtalk_data['xtalk'])

input += args.input_files
pipe.Add(core.G3Reader, filename=input)
pipe.Add(core.DeduplicateMetadata)

# Combine our various in-progress cal data
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

# Invert crosstalk and convert to Watts
if args.xtalkpath is not None:
    pipe.Add(xtalk.CrosstalkCleanedTimestreams, input=args.input_ts, output='BoloMapTimestreams', units=core.G3TimestreamUnits.Watts)
else:
    pipe.Add(dfmux.ConvertTimestreamUnits, Input=args.input_ts, Output='BoloMapTimestreams', Units=core.G3TimestreamUnits.Watts)

## Reset all bolometer pointing that we think we already know to 0 to get boresight-pointing maps
#pipe.Add(std_processing.MakeBoresightBolometerProperties)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(core.Dump)

# Kick off maps
pipe.Add(MakeMap,
         # Timestream filtering -- these really shouldn't be part of MakeMap
         poly_order = 4,
         use_dynamic_source_filter = True, # Enable dynamic PS filtering

         # Actual map-making options
         map_in = smstub,
         map_id = 'bsmap',
         ts_in_key = 'BoloMapTimestreams',

         make_polarized = False,
         do_weight = True,
         fill_in_unity_weights = True,
         use_boresight_pointing = True,
         individual_bolos_to_map = bolos,
         boresight_ra_key = 'BoresightRa',
         boresight_dec_key = 'BoresightDec')

pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

