#!/usr/bin/env python

from spt3g import core, calibration
import sys

# Usage: build_bolo_props.py <inputs> output.g3
# Meta-script to process a bunch of results from other cal observations
# and come up with the quasi-time-independent data (pointing offsets,
# polarization angles, etc.) that live in Cal['BolometerProperties']

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=sys.argv[1:-1])
pipe.Add(calibration.build_cal_frames.BuildBoloPropertiesMap,
         drop_original_frames=True, pointing_cut_threshold=2.*0.00054) # <= 2x pixel-spacing
pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
def AddNotes(fr):
    fr['BolometerPropertiesBasedOn'] = core.G3VectorString(sys.argv[1:-1])
pipe.Add(AddNotes)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=sys.argv[-1])

pipe.Run()

