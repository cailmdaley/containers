import numpy as np
import scipy
import pickle
import glob
import pylab
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

#f1=core.G3File('/poleanalysis/sptdaq/calresult/calibration/boloproperties/32393#000.g3')
#bpframe = f1.next()
#bp = bpframe['BolometerProperties']

fname0 = '/spt_data/bolodata/fullrate/calibrator/34198196/nominal_online_cal.g3'
f1 = core.G3File(fname0)
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

names = bp.keys()
xoffs = np.asarray([bp[name].x_offset/core.G3Units.arcmin for name in names])
yoffs = np.asarray([bp[name].y_offset/core.G3Units.arcmin for name in names])
wafnames = np.asarray([bp[name].physical_name[0:4].upper() for name in names])

pylab.figure()
pylab.plot(xoffs,yoffs,'o')
pylab.xlim(-80,80)
pylab.xlabel('X offset from boresight [arcmin]')
pylab.ylabel('Y offset from boresight [arcmin]')
tctools.annotate_focal_plane()

