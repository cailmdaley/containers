#!/usr/bin/env python

from spt3g import core, calibration
from spt3g.pointing import offline_pointing
import os, glob
import argparse as ap

# Usage: build_cal_frame.py -o output.g3 input.g3
# Meta-script to process a bunch of results from other cal observations
# and come up with the final cal frame for analysis.

P = ap.ArgumentParser(description='Combine calibration frames using'
                      'calibration.build_cal_frames.MergeCalibrationFrames',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
               help='Output filename')
P.add_argument('-c', '--config_dir', action='store', default=None,
               help='Directory of config files for text-based data (offline '
                    'pointing')
P.add_argument('-t', '--time', action='store', default=None,
               help='Observation time for interpolation from config-table '
                    'based modules')
args = P.parse_args()

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
if args.config_dir is not None and os.path.exists(os.path.join(args.config_dir, 'pointing')):
    def grabpointing(fr):
        fr['OfflinePointingModel'] = offline_pointing.OfflinePointingParamsAtTime(args.time, glob.glob(os.path.join(args.config_dir, 'pointing', '*')))
    pipe.Add(grabpointing)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)

pipe.Run()

