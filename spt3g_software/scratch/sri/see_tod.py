import numpy as np, os, sys, glob
import spt3g
from spt3g import core

fname = '/spt/data/bolodata/downsampled/ra0hdec-57.5/11961968/0000.g3'

for frame in core.G3File(fname):
	if not 'RawTimestreams_I' in frame: continue
	print (frame, frame.type)
	all_tod = frame['RawTimestreams_I'] #dictinoary with detector id as keys
	example_tod = all_tod.values()[0]
	sys.exit()
