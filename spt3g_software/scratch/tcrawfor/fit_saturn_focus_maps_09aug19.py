import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
import glob
from spt3g.util import fitting, maths

mapfiles = glob.glob('/spt/user/production/calibration/saturn/maps/8*.g3')
mapfiles = np.asarray(mapfiles)
mapfiles.sort()
obsids = np.asarray([(mf.split('/')[-1]).split('.g3')[0] for mf in mapfiles])
sobsids = np.argsort(obsids)
obsids = obsids[sobsids]
mjds = np.asarray([std_processing.obsid_to_g3time(obsid).mjd for obsid in obsids])
mapfiles = mapfiles[sobsids]
#mapfiles = mapfiles[1:]
#mapfiles = [mapfiles[0]]
nmf = len(mapfiles)
params = np.zeros([3,nmf,7])

npts_model = 300
maps_small = np.zeros([3,nmf,npts_model,npts_model])
xmodel = np.zeros([npts_model,npts_model])
ymodel = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model)
for i in np.arange(npts_model):
    xmodel[i,:] = xtemp
    ymodel[i,:] = i
    
for i in np.arange(len(mapfiles)):
    mfile = mapfiles[i]
    print(mfile)
    f1 = core.G3File(mfile)
    j = 0
    for frame in f1:
        if frame.type is core.G3FrameType.Map:
            mapshape = np.shape(frame['T'])
            map_weighted = frame['T']
            wts = np.asarray(frame['Wunpol'].TT)
            whn0 = np.where(wts > 0.)
            map_unw = (np.asarray(map_weighted)).copy()
            map_unw[whn0] /= wts[whn0]
            xcenter = mapshape[0]//2
            ycenter = mapshape[0]//2
            map_small = map_unw[ycenter-npts_model//2:ycenter+npts_model//2,xcenter-npts_model//2:xcenter+npts_model//2]
            maps_small[j,i,:,:] = map_small
            ptemp = np.asarray(fitting.fit_gaussian2d(map_small))
            modelfunc = maths.gaussian2d(ptemp[0],ptemp[2],ptemp[1],ptemp[4],ptemp[3],ptemp[5],height=ptemp[6])
            model = modelfunc(xmodel,ymodel)
            params[j,i,:] = ptemp
#            print(notavariable)
            j += 1

reso_arcmin = frame['T'].res/core.G3Units.arcmin
ellip = params[:,:,3]/params[:,:,4]
whlt1 = np.where(ellip < 1.)
ellip[whlt1] = 1./ellip[whlt1]
fwhms = np.sqrt(params[:,:,3]*params[:,:,4])*reso_arcmin*np.sqrt(8.*np.log(2.))
x_off_arcsec = (params[:,:,0] - np.float(npts_model/2.))*reso_arcmin*60.
y_off_arcsec = (params[:,:,3] - np.float(npts_model/2.))*reso_arcmin*60.

outd = {}
outd['mjds'] = mjds
outd['ellip'] = ellip
outd['fwhms'] = fwhms
outd['x_off_arcsec'] = x_off_arcsec
outd['y_off_arcsec'] = y_off_arcsec
outfile = '/spt/user/tcrawfor/fit_saturn_focus_maps_09aug19.pkl'
pickle.dump(outd,open(outfile,'wb'))