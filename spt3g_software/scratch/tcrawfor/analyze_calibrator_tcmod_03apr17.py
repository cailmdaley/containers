#!/usr/bin/env python

import numpy, sys, os
from spt3g import core, dfmux, calibration, xtalk
import argparse as ap

# Usage: analyze_calibrator.py -o output.g3 <input files>

P = ap.ArgumentParser(description='Analyze calibrator data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-x', '--xtalkpath', action='store', default=None,
               help='path to directory of crosstalk results')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
P.add_argument('-p', '--phase_type', default='single', action='store',
               help='Do by-bolo phase analysis? (Default is false, set arg to bybolo for True.)')
args = P.parse_args()

pipe = core.G3Pipeline()

# Get obs ID for observation, so we can get the right matrix to do crosstalk inversion
obs = None
for fname in args.input:
    for frame in core.G3File(fname):
        if 'ObservationID' in frame:
            obs = frame['ObservationID']
            break
    if obs is not None:
        break

input = []
if args.xtalkpath is not None:
    xtalk_data = calibration.build_cal_frames.MostRecentCalResultsForObs(obs,
        {'xtalk': os.path.join(args.xtalkpath, 'xtalk-%(fileformatstr)s.g3')})
    input.append(xtalk_data['xtalk'])

input += args.input

pipe.Add(core.G3Reader, filename=input)

# Invert crosstalk and convert to Watts
if args.xtalkpath is not None:
    pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='RawTimestreams_I',
        output='CalTimestreams_I', units=core.G3TimestreamUnits.Watts)
else:
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='RawTimestreams_I',
        Output='CalTimestreams_I', Units=core.G3TimestreamUnits.Watts)
if args.phase_type == 'bybolo':
    pipe.Add(calibration.calibrator.AnalyzeCalibratorData, Input='CalTimestreams_I', CalcPerBoloPhase=True)
else:
    pipe.Add(calibration.calibrator.AnalyzeCalibratorData, Input='CalTimestreams_I')
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(calibration.calibrator.MakeCalibratorFrame)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()
