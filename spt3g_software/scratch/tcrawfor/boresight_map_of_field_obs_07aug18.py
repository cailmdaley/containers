import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp, dfmux
import os
import glob

obsid = '50192225'

# get obs files
ofiles = glob.glob('/spt/data/bolodata/downsampled/ra0hdec-44.75/'+obsid+'/0*.g3')
ofiles.sort()

# other stuff
cframe = core.G3File('/spt/user/production/calibration/calframe/ra0hdec-44.75/'+obsid+'.g3').next()
bp2 = cframe['BolometerProperties']

# get decent bolos
gnames = []
for name in bp2.keys():
    try:
        if cframe['CalibratorResponseSN'][name] > 50. and cframe['ElnodSNSlopes'][name] > 200. and bp2[name].band == 1500.:
            gnames.append(name)
    except:
        pass
gnames = np.asarray(gnames)

#nx = 1200
#ny = 200
nx = 600
ny = 100
ra0 = 0.
dec0 = -45.5
cosdec = np.cos(dec0*np.pi/180.)
#map_weighted = np.zeros(nx*ny)
#map_weights = np.zeros(nx*ny)
map_weighted = np.zeros([ny,nx])
map_weights = np.zeros([ny,nx])
reso_arcmin = 8.0

poly_order = 9

converter = dfmux.unittransforms.ConvertTimestreamUnits(Input='RawTimestreams_I')

for file1 in ofiles:
    f1 = core.G3File(file1)
    for frame in f1:
        if frame.type is core.G3FrameType.Wiring:
            converter(frame)
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
                converter(frame)
                npts = len(frame['OnlineBoresightRa'])
                dtemp = np.zeros([len(gnames),npts])
                nbolo = 0
                for name in gnames:
                    try:
                        dtemp[nbolo,:] = frame['CalTimestreams'][name] - np.mean(frame['CalTimestreams'][name])
                        nbolo += 1
                    except:
                        pass
                dtemp = dtemp[0:nbolo,:]
                commonmode = np.sum(dtemp,0)/np.float(nbolo)
                commonmode -= np.mean(commonmode)
                tcm2 = np.sum(commonmode**2)
                wtemp = np.zeros([nbolo,npts])
                xtemp = np.arange(npts,dtype='float')
                xtemp -= np.mean(xtemp)
                xtemp /= np.max(xtemp)
                for i in np.arange(nbolo):
                    dt2 = dtemp[i,:]
                    coeff1 = np.sum(dt2*commonmode)/tcm2
                    dt2 -= coeff1*commonmode
                    wtemp[i] = 1./(np.std(dt2))**2
                sumts_wtd = np.sum(dtemp*wtemp,0)
                sumts_wts = np.sum(wtemp,0)
                ra_wrap = frame['OnlineBoresightRa']/core.G3Units.deg
                ra_wrap[np.where(ra_wrap > 180.)] -= 360.
                xptemp = (np.round((ra_wrap-ra0)*cosdec*60./reso_arcmin)).astype(int) + nx/2
                yptemp = (np.round((frame['OnlineBoresightDec']/core.G3Units.deg-dec0)*60./reso_arcmin)).astype(int) + ny/2
#                iptemp = yptemp*nx + xptemp
                for i in np.arange(len(xptemp)):
                    map_weighted[yptemp[i],xptemp[i]] += sumts_wtd[i]
                    map_weights[yptemp[i],xptemp[i]] += sumts_wts[i]
#                for i in np.arange(len(iptemp)):
#                    map_weighted[iptemp[i]] += sumts_wtd[i]
#                    map_weights[iptemp[i]] += sumts_wts[i]
#
#    print(notavariable)

                



#                    coeffs2 = np.polyfit(xtemp,dt2,poly_order)
#                    print(notavariable)
#                    for j in np.arange(poly_order):
#                        dt2 -= xtemp**(poly_order-j-1)

wt2 = map_weights.copy()
wt2[np.where(map_weights == 0.)] = 1.
map_unw = map_weighted.copy()
map_unw /= wt2
