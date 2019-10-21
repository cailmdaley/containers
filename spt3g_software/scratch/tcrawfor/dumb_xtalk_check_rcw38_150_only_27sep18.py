from spt3g import core, dfmux, mapmaker
from spt3g.util import tctools
import numpy as np
import pickle
import glob

obsid = '54036687'

f1 = core.G3File('/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsid+'.g3')

mapdict = {}

for frame in f1:
    if frame.type is core.G3FrameType.Calibration:
        bp = frame['BolometerProperties']
        names = []
        calsndict = frame['CalibratorResponseSN']
        for name in calsndict.keys():
            if calsndict[name] > 20.:
                names.append(name)
        nbolo = len(names)

        if frame.type is core.G3FrameType.Map:
#            if frame['Id'] != 'bsmap' and frame['Id'] != 'azmap' and frame['Id'] != 'elmap':
            if 'Wunpol' in frame:
                try:
                    twt = np.asarray(frame['Wunpol'].TT)
                    twt2 = twt.copy()
                    twt2[np.where(twt == 0.)] = 1.
                except:
                    pass
            if frame['Id'] in names:
                if bp[frame['Id']].band/10. == 150.:
                    thismap = np.asarray(frame['T'])
                    thismap /= twt2
                    mapdict[frame['Id']] = thismap
            print(notavariable)

mapnames = mapdict.keys()
nmap = len(mapnames)
xmin = 50
xmax = 310
ymin = 50
ymax = 310
xtalk_matrix = np.zeros([nmap,nmap])
for i in np.arange(nmap):
    mapi = mapdict[mapnames[i]]
    for j in np.arange(nmap-i) + i:
        mapj = mapdict[mapnames[j]]
        xtalk_matrix[i,j] = np.sum(mapi[ymin:ymax,xmin:xmax]*mapj[ymin:ymax,xmin:xmax])
        xtalk_matrix[j,i] = xtalk_matrix[i,j]

xtalk_matrix_unnorm = xtalk_matrix.copy()
for i in np.arange(nmap):
    for j in np.arange(nmap):
        xtalk[i,j] /= np.sqrt(xtalk[i,i]*xtalk[j,j])
