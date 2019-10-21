import numpy as np
import scipy
from scipy import ndimage, interpolate
import pickle
import argparse as ap
from spt3g import core, std_processing, gcp
import os
import glob
from spt3g.util import tctools, filtered_gauss_fit

#mapdir = '/poleanalysis/tcrawfor/radio_day_2019/'
mapdir = '/poleanalysis/tcrawfor/radio_day_2019_nopm_taken_out/'
mapfiles = glob.glob(mapdir + '*/*.g3')
srcnames = np.asarray([mf.split('/')[-2] for mf in mapfiles])
srcdecs = np.asarray([np.float((ntemp.split('-')[-1])[0:2]) for ntemp in srcnames])
asd = np.argsort(srcdecs)
mapfiles = (np.asarray(mapfiles)[asd])[::-1]
srcnames = (srcnames[asd])[::-1]
srcdecs = (srcdecs[asd])[::-1]

nmf = len(mapfiles)
params = np.zeros([3,nmf,7])
wnoise = np.zeros([3,nmf])
crudesn = np.zeros([3,nmf])
fwhms = np.zeros([3,nmf])
pos_unc_arcsec = np.zeros([3,nmf])
ra_src_meas = np.zeros([3,nmf])
dec_src_meas = np.zeros([3,nmf])

npts_model = 80
srange_low = 100
srange_high = 260
maps_small = np.zeros([3,nmf,npts_model,npts_model])
models_small = np.zeros([3,nmf,npts_model,npts_model])
xmodel = np.zeros([npts_model,npts_model])
ymodel = np.zeros([npts_model,npts_model])
xtemp = np.arange(npts_model,dtype=int)
for i in np.arange(npts_model):
    xmodel[i,:] = xtemp
    ymodel[i,:] = i
    
# get associated calfiles
f2 = open('/home/tcrawfor/spt_code/spt3g_software/calibration/scripts/fluxpointcal/radio_day_60000000.scr','r')
lines = f2.readlines()
calfiles = []
caldir = '/poleanalysis/sptdaq/calresult/calibration/calibrator/'
for mfile in mapfiles:
    thisobsid_str = mfile.split('/')[-1].split('.g3')[0]
    for line in lines:
        if thisobsid_str in line:
            thiscalobsid_str = line.split(caldir)[-1].split('.g3')[0]
            calfiles.append(caldir+thiscalobsid_str+'.g3')
medcals = np.zeros(nmf)

obsids = []
for i in np.arange(nmf):
    mfile = mapfiles[i]
    print(mfile)
    thisobsid_str = mfile.split('/')[-1].split('.g3')[0]
    obsids.append(np.int(thisobsid_str))
    fstemp = os.stat(mfile)
    if fstemp.st_size < 100:
        continue
    f1 = core.G3File(mfile)
    cframe = core.G3File(calfiles[i]).next()
    thesecals = []
    thesecsn = []
    for name in cframe['CalibratorResponse'].keys():
        if np.isfinite(cframe['CalibratorResponse'][name]) and np.isfinite(cframe['CalibratorResponseSN'][name]):
            thesecals.append(cframe['CalibratorResponse'][name])
            thesecsn.append(cframe['CalibratorResponseSN'][name])
    thesecals_good = np.asarray(thesecals)[np.where(np.asarray(thesecsn) > 20.)]
    medcal = np.median(thesecals_good)
    medcals[i] = medcal
    do_invert_T = medcal < 0.
    for j in np.arange(3):
        frame = f1.next()
        map_weighted = np.asarray(frame['T'])
        if do_invert_T:
            map_weighted *= -1.
        mapshape = np.shape(map_weighted)
        xcenter = mapshape[0]//2
        ycenter = mapshape[0]//2
        if i == 0 and j == 0:
            maps_weighted = np.zeros([3,nmf,mapshape[0],mapshape[1]])          
            maps_unw = np.zeros([3,nmf,mapshape[0],mapshape[1]])
        maps_weighted[j,i,:,:] = map_weighted
        wts = np.asarray(frame['Wunpol'].TT)
        whn0 = np.where(wts > 0.)
        map_unw = (np.asarray(map_weighted)).copy()
        map_unw[whn0] /= wts[whn0]
        map_small = map_unw[ycenter-npts_model//2:ycenter+npts_model//2,xcenter-npts_model//2:xcenter+npts_model//2]
        maps_small[j,i,:,:] = map_small
        maps_unw[j,i,:,:] = map_unw
        ptemp = np.asarray(filtered_gauss_fit.fitFiltered2Dgaussian(map_small))
        modelfunc = filtered_gauss_fit.twoDgaussian(ptemp[0],ptemp[1],ptemp[2],ptemp[3],ptemp[4],ptemp[5],ptemp[6])
        model = filtered_gauss_fit.hpfilt(modelfunc(*np.indices(map_small.shape)))
        models_small[j,i,:,:] = model
        wnoise_offsource = np.nanstd((map_small-model)[0:10,0:10])
        wnoise[j,i] = wnoise_offsource
        crudesn[j,i] = np.sqrt(np.nansum((model/wnoise_offsource)**2))
        reso_arcmin = frame['T'].res/core.G3Units.arcmin
        fwhms[j,i] = np.sqrt(np.abs(ptemp[4])*np.abs(ptemp[5]))*reso_arcmin*np.sqrt(8.*np.log(2.))
        pos_unc_arcsec[j,i] = fwhms[j,i]*60./crudesn[j,i]
        params[j,i,:] = ptemp
        if j == 0:
            ra_grid = np.zeros([npts_model,npts_model])
            dec_grid = np.zeros([npts_model,npts_model])
            for k in np.arange(npts_model):
                for q in np.arange(npts_model):
                    angstemp = frame['T'].pixel_to_angle(np.int(k+xcenter-npts_model//2),np.int(q+ycenter-npts_model//2))
                    ra_grid[q,k] = angstemp[0]/core.G3Units.deg
                    dec_grid[q,k] = angstemp[1]/core.G3Units.deg
        try:
            ra_src_meas[j,i] = scipy.interpolate.interpn((xtemp,xtemp),ra_grid,np.asarray([ptemp[3],ptemp[2]]))
            dec_src_meas[j,i] = scipy.interpolate.interpn((xtemp,xtemp),dec_grid,np.asarray([ptemp[3],ptemp[2]]))
        except:
            pass
#        print(notavariable)

ellip = params[:,:,4]/params[:,:,5]
whlt1 = np.where(ellip < 1.)
ellip[whlt1] = 1./ellip[whlt1]
x_off_arcsec = (params[:,:,2] - (np.float(npts_model)/2.-0.5))*reso_arcmin*60.
y_off_arcsec = (params[:,:,3] - (np.float(npts_model)/2.-0.5))*reso_arcmin*60.

outd = {}
outd['ampls'] = params[:,:,1]
outd['ellip'] = ellip
outd['fwhms'] = fwhms
outd['x_off_arcsec'] = x_off_arcsec
outd['y_off_arcsec'] = y_off_arcsec
outd['srcnames'] = srcnames
outd['obsids'] = obsids
outd['crudesn'] = crudesn
outd['pos_unc_arcsec'] = pos_unc_arcsec
outd['ra_src_meas'] = ra_src_meas
outd['dec_src_meas'] = dec_src_meas
userdir = os.path.expanduser('~/')
#outfile = userdir + '/fit_radio_day_maps_01feb19.pkl
#outfile = userdir + '/fit_radio_day_maps_02feb19.pkl
#outfile = userdir + '/fit_radio_day_maps_02feb19_newboloprops_nopm.pkl'
outfile = userdir + '/fit_radio_day_maps_04feb19_newboloprops_nopm.pkl'
pickle.dump(outd,open(outfile,'wb'))
