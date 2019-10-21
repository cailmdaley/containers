#!/usr/bin/env python

import numpy, sys, os
from spt3g import core, dfmux, calibration, std_processing, dfmux
import argparse as ap

# Usage: analyze_calibrator.py -o output.g3 <input files>

P = ap.ArgumentParser(description='Analyze calibrator data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
P.add_argument('-p', '--phaseIandQ', action='store_true', default=False,
               help='Analyze I and Q separately (for phase rotation purposes)')
args = P.parse_args()

pipe = core.G3Pipeline()

def drop_start_stop_glitches(frame, cut_time = .25 * core.G3Units.s):
    '''
    The docstring below is false, but preserved for posterity.  The
    actual problem was with the downsampler causing huge ringing,
    which was fixed later.

    Due to the time resolution of feature bits in GCP, we get a little
    bit of time before and after the shutter opens.  This leads to
    large DC glitches that the calibrator analysis locks into.  We fix
    this by cutting a little bit of time on either end (.25 s by
    default).
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    new_tsm = core.G3TimestreamMap()
    cal = frame.pop('CalibratorOn', None)
    cut = int(cal.sample_rate * cut_time)
    frame['CalibratorOn'] = cal[cut:-cut]
    for name, bolo in frame['RawTimestreams_I'].iteritems():
        new_tsm[name] = bolo[cut:-cut]
    
    frame.pop('RawTimestreams_I', None)
    frame['RawTimestreams_I'] = new_tsm
    return frame

pipe.Add(core.G3Reader, filename=args.input)
# pipe.Add(drop_start_stop_glitches)

# Interpolate over NaNs
pipe.Add(std_processing.InterpOverNans)

if args.phaseIandQ:
    # SG Jan 19: only convert to Watts. IQ Rotation happens in AnalyzeCalibratorData. Can crosstalk inversion be done pre-rotation?
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='RawTimestreams_I',Output='CalTimestreamsI',Units=core.G3TimestreamUnits.Power)
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='RawTimestreams_Q',Output='CalTimestreamsQ',Units=core.G3TimestreamUnits.Power)
    Input_Q = 'CalTimestreamsQ'

else:
    # Invert crosstalk, apply Elnod IQ rotation, and convert to Watts
    pipe.Add(std_processing.CalibrateRawTimestreams,
             units=core.G3TimestreamUnits.Power,
             keep_original=False, calibrator_iq_rotation=False)
    Input_Q = None

# Drop calibration frames so that we don't pollute the output file
#pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(calibration.calibrator.AnalyzeCalibratorData, Input='CalTimestreams', Input_Q = Input_Q, Do_IQ_Rotation = args.phaseIandQ)
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(calibration.calibrator.MakeCalibratorFrame)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run(profile=True)
