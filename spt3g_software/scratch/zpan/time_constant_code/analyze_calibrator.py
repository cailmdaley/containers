#!/usr/bin/env python

import numpy, sys, os
from spt3g import core, dfmux, calibration, std_processing
import argparse as ap
# use the customized calibrator_analysis in the same folder
import calibrator_analysis

# Usage: analyze_calibrator.py -o output.g3 <input files>

P = ap.ArgumentParser(description='Analyze calibrator data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
args = P.parse_args()

pipe = core.G3Pipeline()

def drop_start_stop_glitches(frame, cut_time = .25 * core.G3Units.s):
    '''
    The docstring below is false, but preserver for posterity.  The
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

# Invert crosstalk, apply IQ rotation, and convert to Watts
pipe.Add(std_processing.CalibrateRawTimestreams,
         units=core.G3TimestreamUnits.Power)
# Drop calibration frames so that we don't pollute the output file
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q'])

pipe.Add(calibrator_analysis.AnalyzeCalibratorData, Input='CalTimestreams')
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(calibrator_analysis.MakeCalibratorFrame)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()
