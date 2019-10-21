#!/usr/bin/env python

import numpy as np
from spt3g import core, calibration, std_processing,timestreamflagging
from spt3g.std_processing import flagsegments
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

pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

pipe.Add(std_processing.CalibrateRawTimestreams, q_output='CalTimestreamsQ')

pipe.Add(lambda fr: 'CalTimestreams' in fr)

pipe.Add(calibration.elnod.AnalyzeElnod,
         i_data_key='CalTimestreams',
         q_data_key='CalTimestreamsQ')

pipe.Add(lambda fr: 'ElnodSlopes' in fr)
pipe.Add(core.Rename, keys={'ElnodSlopes': 'KcmbElnodSlopes'})
pipe.Add(core.Delete, keys=['ElnodDrifts','ElnodEigenvalueDominantVectorI','ElnodEigenvalueDominantVectorQ','ElnodEigenvalueRatio','ElnodRSquared','ElnodSNSlopes','ElnodSigmaSlopes','ElnodSlopes','ElnodVariances'])
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

