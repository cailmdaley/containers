#!/usr/bin/env python

import numpy as np
from spt3g import core, calibration, std_processing
import argparse as ap

# Compute elnod slope in K_cmb*arb/airmass

# Usage: elnod-kcmb.py -o output.g3 <input files>

P = ap.ArgumentParser(description='Analyze elnod data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
args = P.parse_args()

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(std_processing.CalibrateRawTimestreams)
pipe.Add(calibration.SplitTimestreamsByBand, input='CalTimestreams',
    output_root='CalTimestreams')
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)
pipe.Add(calibration.elnod.AnalyzeElnod, i_data_key='CalTimestreams150GHz',
         q_data_key='')
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
pipe.Add(core.Rename, keys={'ElnodSNSlopes': 'ElnodKcmbSNSlopes', 'ElnodSlopes': 'ElnodKcmbSlopes', 'ElnodRSquared': 'ElnodKcmbRSquared'})
pipe.Add(core.Delete, keys=['ElnodDrifts', 'ElnodEigenvalueDominantVectorI', 'ElnodEigenvalueDominantVectorQ', 'ElnodEigenvalueRatio', 'ElnodSigmaSlopes', 'ElnodVariances'])
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()
