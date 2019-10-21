#!/usr/bin/env python

import numpy as np
from spt3g import core, calibration, std_processing, mapmaker, timestreamflagging
import argparse as ap

# Usage: analyze_elnod.py -o output.g3 <input files>

P = ap.ArgumentParser(description='Analyze elnod data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
P.add_argument('--i-data-key', default='RawTimestreams_I', action='store',
               help='I Timestream key to read')
P.add_argument('--q-data-key', default=None, action='store',
               help='Q Timestream key to read')
args = P.parse_args()

pipe = core.G3Pipeline()

if args.q_data_key is None:
    if args.i_data_key.endswith('_I'):
        args.q_data_key = args.i_data_key[:-2] + '_Q'
    else:
        raise ValueError('Unable to automatically determine Q data key')

pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(std_processing.InterpOverNans)
pipe.Add(mapmaker.TodFiltering, ts_in_key = args.i_data_key, 
        ts_out_key = args.i_data_key + '_LPF', poly_order = 0, 
        lpf_filter_frequency = 5.*core.G3Units.Hz)
pipe.Add(mapmaker.TodFiltering, ts_in_key = args.q_data_key, 
        ts_out_key = args.q_data_key + '_LPF', poly_order = 0, 
        lpf_filter_frequency = 5.*core.G3Units.Hz)
pipe.Add(calibration.elnod.AnalyzeElnod,
         i_data_key=args.i_data_key + '_LPF',
         q_data_key=args.q_data_key + '_LPF')
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration or fr.type == core.G3FrameType.EndProcessing)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()
