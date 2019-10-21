import numpy as np
import scipy
from scipy import ndimage
import pickle, glob
from spt3g import core, std_processing, gcp
from spt3g.util import tctools

ny = 360
nx = 360

reso_arcmin = 0.5

sources = np.asarray([
'J1617-5848', 
'J0904-5735', 
'J1329-5608', 
'J1326-5256', 
'J0210-5101', 
'J2235-4835', 
'J1107-4449', 
'J1427-4206', 
'J0522-3627', 
'J1626-2951', 
'J1037-2934', 
'J0453-2807', 
'J2258-2758', 
'J1337-1257', 
'J1512-0905', 
'J0006-0623'])

nobs = len(sources)
mdir = '/poleanalysis/tcrawfor/calresult/calibration/'

maps = np.zeros([nobs,3,ny,nx])
wm90 = np.zeros([nobs,ny,nx])
wm150 = np.zeros([nobs,ny,nx])
wm220 = np.zeros([nobs,ny,nx])
obsids = []
mjds = []
for i in np.arange(nobs):
    thisfile = glob.glob(mdir + sources[i] + '/*.g3')
    if len(thisfile) > 1:
        thisfile.sort()
#        thisfile = thisfile[-1]
        thisfile = thisfile[-2]
    else:
        thisfile = thisfile[0]
    obsid = (thisfile.split('/'))[-1].split('.')[0]
    obsids.append(obsid)
    mjds.append(std_processing.obsid_to_g3time(obsid).mjd)
    f1 = core.G3File(thisfile)
    for j in np.arange(3):
        frame = f1.next()
        mt1=np.asarray(frame['T'])
        if j == 0: wm90[i,:,:] = mt1
        if j == 1: wm150[i,:,:] = mt1
        if j == 2: wm220[i,:,:] = mt1
        wt1=np.asarray(frame['Wunpol'].TT)
        wh1=np.where(wt1 > 0)
        mt1[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt1

print(notavariable)

isgood90 = [0, 1, 2, 3, 4, 5, 6, 8]#, 12, 13]
isgood150 = [0, 1, 2, 3, 4, 5, 6]

xmin = 130
xmax = 230
npts = 200
npts_ampl = 30
xvec = np.arange(xmax-xmin)
xarr = np.zeros([xmax-xmin,xmax-xmin])
yarr = np.zeros([xmax-xmin,xmax-xmin])
for i in np.arange(xmax-xmin):
    xarr[i,:] = xvec
    yarr[i,:] = xvec[i]
centervec = (np.arange(npts)/float(npts)-0.5)*4. + (xmax-xmin)/2.
amplvec = (np.arange(npts_ampl)/float(npts_ampl)-0.5)*0.4
goodmaps = np.zeros([len(isgood90),xmax-xmin,xmax-xmin])
models = np.zeros([len(isgood90),xmax-xmin,xmax-xmin])
ybest = np.zeros(len(isgood90))
xbest = np.zeros(len(isgood90))
for i in np.arange(len(isgood90)):
    j = isgood90[i]
#    thismap = ndimage.gaussian_filter(maps[j,0,xmin:xmax,xmin:xmax],2)
    thismap = maps[j,0,xmin:xmax,xmin:xmax]
    if i in isgood150:
#        thismap += ndimage.gaussian_filter(maps[j,1,xmin:xmax,xmin:xmax],3)*6.
        thismap += ndimage.gaussian_filter(maps[j,1,xmin:xmax,xmin:xmax],1)*6.
    chitemp = np.zeros([npts_ampl,npts,npts])
    for q in np.arange(npts):
        cenx = centervec[q]
        xrt2 = (xarr - cenx)**2
        for qq in np.arange(npts):
            ceny = centervec[qq]
            yrt2 = (yarr - ceny)**2
#            model_norm = np.exp(-(xrt2+yrt2)/2./1.7**2)
            model_norm = np.exp(-(xrt2+yrt2)/2./1.4**2)
            for qqq in np.arange(npts_ampl):
                model = amplvec[qqq]*model_norm
                resid = thismap - model
                chitemp[qqq,qq,q] = np.sum(resid**2)
    qqqmin,qqmin,qmin = np.unravel_index(np.argmin(chitemp),(npts_ampl,npts,npts))
    ybest[i] = centervec[qqmin]
    xbest[i] = centervec[qmin]
    goodmaps[i,:,:] = thismap
#    models[i,:,:] = amplvec[qqqmin]*np.exp(-((yarr-ybest[i])**2+(xarr-xbest[i])**2)/2./1.7**2)
    models[i,:,:] = amplvec[qqqmin]*np.exp(-((yarr-ybest[i])**2+(xarr-xbest[i])**2)/2./1.4**2)
