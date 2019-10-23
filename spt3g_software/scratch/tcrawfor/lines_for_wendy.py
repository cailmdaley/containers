import numpy as np
import scipy
import pickle
import glob
import pylab
from spt3g import core, std_processing, gcp
from spt3g.scratch.tcrawfor import tctools

dir1 = '/poleanalysis/sptdaq/20170129_rawdata/'
fname0 = dir1 + '20170129_194237.g3'
f1 = core.G3File(fname0)
wframe = f1.next()
wmap = wframe['WiringMap']
names = np.asarray(wmap.keys())
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

plate_scale_correction = 42./15.
xoffs = np.asarray([bp[name].x_offset*42./15./core.G3Units.arcmin for name in names])
yoffs = np.asarray([bp[name].y_offset*42./15./core.G3Units.arcmin for name in names])

pylab.figure()
pylab.plot(xoffs,yoffs,'o')
pylab.xlabel('X offset from boresight [arcmin]')
pylab.ylabel('Y offset from boresight [arcmin]')
tctools.annotate_focal_plane()

