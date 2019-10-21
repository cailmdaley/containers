from spt3g import core, dfmux, std_processing, todfilter
from spt3g.util import tctools
import numpy as np
import pickle
from scipy import ndimage

nframe = core.G3File('/spt/user/production/calibration/noise/64037271.g3').next()
#nframe = core.G3File('/spt/user/production/calibration/noise/66705025.g3').next()
netdict = nframe['NET_10.0Hz_to_15.0Hz']
names = netdict.keys()
bpframe = core.G3File('/spt/user/production/calibration/boloproperties/60000000.g3').next()
bp = bpframe['BolometerProperties']
wafnames = np.asarray([bp[name].wafer_id for name in names])

bands = np.zeros(len(names))
xoffs = np.zeros(len(names))
yoffs = np.zeros(len(names))
nets = np.zeros(len(names))
for i in np.arange(len(names)):
    name = names[i]
    bands[i] = bp[name].band/10.
    xoffs[i] = bp[name].x_offset/core.G3Units.arcmin
    yoffs[i] = bp[name].y_offset/core.G3Units.arcmin
    if 'y' in bp[name].physical_name:
        xoffs[i] += 1.
        yoffs[i] += 1.
    nets[i] = netdict[name]

wh90 = np.where(bands == 90.)
wh150 = np.where(bands == 150.)
wh220 = np.where(bands == 220.)

xoffs90 = xoffs[wh90]
yoffs90 = yoffs[wh90]
nets90 = nets[wh90]
xoffs150 = xoffs[wh150]
yoffs150 = yoffs[wh150]
nets150 = nets[wh150]
xoffs220 = xoffs[wh220]
yoffs220 = yoffs[wh220]
nets220 = nets[wh220]

