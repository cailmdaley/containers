import numpy as np
import scipy
from scipy import ndimage
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
import glob
from spt3g.util import tctools, gauss_fit

mapfiles = glob.glob('/poleanalysis/sptdaq/calresult/calibration/RCW38/maps/6*.g3')
mapfiles.sort()
mapfiles = mapfiles[8:]
#mapfiles = [mapfiles[0]]
obsids = np.asarray([(mf.split('/')[-1]).split('.g3')[0] for mf in mapfiles])
mjds = np.asarray([std_processing.obsid_to_g3time(obsid).mjd for obsid in obsids])
nmf = len(mapfiles)
params = np.zeros([3,nmf,7])

#print(notavariable)

npts_model = 40
maps_small = np.zeros([3,nmf,npts_model,npts_model])
xmodel = np.zeros([npts_model,npts_model])
ymodel = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model,dtype=int)
for i in np.arange(npts_model):
    xmodel[i,:] = xtemp
    ymodel[i,:] = i
    
for i in np.arange(len(mapfiles)):
    mfile = mapfiles[i]
#    print mfile
    f1 = core.G3File(mfile)
    for j in np.arange(3):
        frame = f1.next()
        map_weighted = frame['T']
        wts = np.asarray(frame['Wunpol'].TT)
        whn0 = np.where(wts > 0.)
        map_unw = (np.asarray(map_weighted)).copy()
        map_unw[whn0] /= wts[whn0]
        if np.max(map_unw[160:200,160:200]) < np.abs(np.min(map_unw[160:200,160:200])):
            map_unw = -map_unw
        mapshape = np.shape(map_unw)
#        print mapshape
#        ycenter, xcenter = np.unravel_index(np.argmax(map_unw[160:200,160:200]),[40,40])
#        ycenter += 160
#        xcenter += 160
        xcenter = mapshape[0]//2
        ycenter = mapshape[0]//2
        map_small = map_unw[ycenter-npts_model//2:ycenter+npts_model//2,xcenter-npts_model//2:xcenter+npts_model//2]
        maps_small[j,i,:,:] = map_small
        ptemp = np.asarray(gauss_fit.fit2Dgaussian(map_small))
        modelfunc = gauss_fit.twoDgaussian(ptemp[0],ptemp[1],ptemp[2],ptemp[3],ptemp[4],ptemp[5],ptemp[6])
        model = modelfunc(xmodel,ymodel)
        params[j,i,:] = ptemp
#        if j == 2:
#            print(notavariable)


reso_arcmin = frame['T'].res/core.G3Units.arcmin
ellip = params[:,:,4]/params[:,:,5]
whlt1 = np.where(ellip < 1.)
ellip[whlt1] = 1./ellip[whlt1]
fwhms = np.sqrt(params[:,:,4]*params[:,:,5])*reso_arcmin*np.sqrt(8.*np.log(2.))
x_off_arcsec = (params[:,:,2] - np.float(npts_model/2.))*reso_arcmin*60.
y_off_arcsec = (params[:,:,3] - np.float(npts_model/2.))*reso_arcmin*60.

outd = {}
outd['mjds'] = mjds
outd['ellip'] = ellip
outd['fwhms'] = fwhms
outd['x_off_arcsec'] = x_off_arcsec
outd['y_off_arcsec'] = y_off_arcsec
#outfile = '/home/tcrawfor/fit_rcw38_maps_15jan19.pkl'
#outfile = '/home/tcrawfor/fit_rcw38_maps_24jan19.pkl'
outfile = '/home/tcrawfor/fit_rcw38_maps_28jan19.pkl'
pickle.dump(outd,open(outfile,'wb'))