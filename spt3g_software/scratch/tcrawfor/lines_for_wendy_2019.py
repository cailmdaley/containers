import numpy as np
import scipy
import pickle
import glob
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

#bp = core.G3File('/poleanalysis/sptdaq/calresult/calibration/boloproperties/62679603.g3').next()['BolometerProperties']
fname0 = '/spt_data/bolodata/fullrate/calibrator/66579096/nominal_online_cal.g3'
f1 = core.G3File(fname0)
cframe = f1.next()
bp = cframe['NominalBolometerProperties']

names = bp.keys()
xoffs = np.asarray([bp[name].x_offset/core.G3Units.arcmin for name in names])
yoffs = np.asarray([bp[name].y_offset/core.G3Units.arcmin for name in names])
wafnames = np.asarray([bp[name].physical_name[0:4].upper() for name in names])

plot(xoffs,-yoffs,'o')
xlim(-82,82)
ylim(-62,62)
xlabel('X offset from boresight [arcmin]')
ylabel('Y offset from boresight [arcmin]')
tctools.annotate_focal_plane()

