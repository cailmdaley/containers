#!/usr/bin/env python
import sys
from spt3g import core

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=sys.argv[1])


def calcmean(frame):
    detmeans = core.G3MapDouble()
    if 'RawTimestreams_I' in frame:
        detlist = frame['RawTimestreams_I'].keys()
    for det in detlist:
        bolodata = frame['RawTimestreams_I'][det]
        thismean = bolodata.mean()
        detmeans[det] = thismean
        print(thismean)
    frame['means'] = detmeans


pipe.Add(calcmean)

pipe.Add(core.G3Writer, filename=sys.argv[2])

pipe.Run()
