import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
import glob
from spt3g.scratch.tcrawfor import tctools, gauss_fit

mapfiles = glob.glob('/home/nwhitehorn/focus-quasar*.g3')
mapfiles.sort()
mapfiles = mapfiles[1:]
#mapfiles = [mapfiles[0]]
nmf = len(mapfiles)
params = np.zeros([nmf,7])
mapfiles_fineres = glob.glob('/home/nwhitehorn/focus-quasar*fine*.g3')
mapfiles_fineres.sort()

npts_model = 20
xmodel = np.zeros([npts_model,npts_model])
ymodel = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model,dtype=int)
for i in np.arange(npts_model):
    xmodel[i,:] = xtemp
    ymodel[i,:] = i
    
for i in np.arange(len(mapfiles)):
    mfile = mapfiles[i]
    f1 = core.G3File(mfile)
    for frame in f1:
        if frame.type == core.G3FrameType.Map:
            mframe = frame
    map_weighted = mframe['T']
    wts = np.asarray(mframe['Wunpol'].TT)
    whn0 = np.where(wts > 0.)
    map_unw = (np.asarray(map_weighted)).copy()
    map_unw[whn0] /= wts[whn0]
    map_unw = -map_unw
    mapshape = np.shape(map_unw)
    ycenter, xcenter = np.unravel_index(np.argmax(map_unw[150:280,120:250]),[130,130])
    ycenter += 150
    xcenter += 120
# kludge
    if mfile == '/home/nwhitehorn/focus-quasar.g3':
        xcenter = 191
        ycenter = 216
    map_small = map_unw[ycenter-npts_model/2:ycenter+npts_model/2,xcenter-npts_model/2:xcenter+npts_model/2]
    ptemp = np.asarray(gauss_fit.fit2Dgaussian(map_small))
    modelfunc = gauss_fit.twoDgaussian(ptemp[0],ptemp[1],ptemp[2],ptemp[3],ptemp[4],ptemp[5],ptemp[6])
    model = modelfunc(xmodel,ymodel)
    params[i,:] = ptemp

reso_arcmin = 0.5
ellip = params[:,4]/params[:,5]
whlt1 = np.where(ellip < 1.)
ellip[whlt1] = 1./ellip[whlt1]
fwhms = np.sqrt(params[:,4]*params[:,5])*reso_arcmin*np.sqrt(8.*np.log(2.))

mfile = mapfiles_fineres[0]
f1 = core.G3File(mfile)
for frame in f1:
    if frame.type == core.G3FrameType.Map:
        mframe = frame
map_weighted = mframe['T']
reso_arcmin = map_weighted.res/core.G3Units.arcmin
wts = np.asarray(mframe['Wunpol'].TT)
whn0 = np.where(wts > 0.)
map_unw = (np.asarray(map_weighted)).copy()
map_unw[whn0] /= wts[whn0]
map_unw = -map_unw
mapshape = np.shape(map_unw)
map_sm = ndimage.gaussian_filter(map_unw[680:760,590:670],8)
ycenter, xcenter = np.unravel_index(np.argmax(map_sm),[80,80])
ycenter += 680
xcenter += 590
npts_model = 40
xmodel2 = np.zeros([npts_model,npts_model])
ymodel2 = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model,dtype=int)
map_small = map_unw[ycenter-npts_model/2:ycenter+npts_model/2,xcenter-npts_model/2:xcenter+npts_model/2]
#map_small = ndimage.gaussian_filter(map_small,2)
ptemp2 = np.asarray(gauss_fit.fit2Dgaussian(map_small))
modelfunc2 = gauss_fit.twoDgaussian(ptemp2[0],ptemp2[1],ptemp2[2],ptemp2[3],ptemp2[4],ptemp2[5],ptemp2[6])
model2 = modelfunc2(xmodel2,ymodel2)
params_fineres = ptemp2
ellip_fineres = params_fineres[4]/params_fineres[5]
fwhm_fineres = np.sqrt(params_fineres[4]*params_fineres[5])*reso_arcmin*np.sqrt(8.*np.log(2.))

