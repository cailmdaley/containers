from spt3g import core, dfmux, mapmaker
import numpy as np
import pickle
import scipy
from scipy import ndimage

rfiles=glob.glob('/spt/user/production/calibration/RCW38-pixelraster/maps/4*.g3')
rfiles.sort()
nrfiles = len(rfiles)

id1 = "005.10.1.1.1714"

mjds = np.zeros(nrfiles)
xps = np.zeros(nrfiles)
yps = np.zeros(nrfiles)

for i in np.arange(nrfiles):
    f1 = core.G3File(rfiles[i])
    for j in np.arange(10):
        frame = f1.next()
        if frame.type is core.G3FrameType.Observation:
            mjds[i] = frame['ObservationStart'].mjd
        if frame.type is core.G3FrameType.Map:
            if frame['Id'] == id1:
                mapshape = np.shape(frame['T'])
                amap_sm = ndimage.gaussian_filter(frame['T'],2.)
                ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
                xps[i] = xcenter
                yps[i] = ycenter
