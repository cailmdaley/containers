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

decstemp = np.asarray([int(name[-5:-2]) for name in sources])
sdt = np.argsort(decstemp)
sources = sources[sdt]

nobs = len(sources)
mdir = '/poleanalysis/tcrawfor/calresult/calibration/'

maps = np.zeros([nobs,3,ny,nx])
wm90 = np.zeros([nobs,ny,nx])
wm150 = np.zeros([nobs,ny,nx])
wt150 = np.zeros([nobs,ny,nx])
wm220 = np.zeros([nobs,ny,nx])
obsids = []
mjds = []
obsid_min = 33480000
obsid_max = 33530000
for i in np.arange(nobs):
    thesefiles = glob.glob(mdir + sources[i] + '/calhack/*.g3')
    if len(thesefiles) > 1:
        for file1 in thesefiles:
            obsid1 = np.int((file1.split('/'))[-1].split('.')[0])
            if obsid1 > obsid_min and obsid1 < obsid_max:
                thisfile = file1
    else:
        thisfile = thesefiles[0]
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
        if j == 1: wt150[i,:,:] = np.asarray(frame['Wunpol'].TT)
        wh1=np.where(wt1 > 0)
        mt1[wh1]/=wt1[wh1]
        maps[i,j,:,:] = mt1

print(notavariable)

isgood90 = np.arange(11)
isgood150 = np.arange(6)

xmin = 130
xmax = 230
#npts = 200
#npts_ampl = 40
npts = 20
npts_ampl = 10
xvec = np.arange(xmax-xmin)
xarr = np.zeros([xmax-xmin,xmax-xmin])
yarr = np.zeros([xmax-xmin,xmax-xmin])
for i in np.arange(xmax-xmin):
    xarr[i,:] = xvec
    yarr[i,:] = xvec[i]
centervec = (np.arange(npts)/float(npts)-0.5)*10. + (xmax-xmin)/2. - 2.
amplvec = (np.arange(npts_ampl)/float(npts_ampl)-0.5)*0.6 + 1.0
goodmaps = np.zeros([len(isgood90),xmax-xmin,xmax-xmin])
models = np.zeros([len(isgood90),xmax-xmin,xmax-xmin])
ybest = np.zeros(len(isgood90))
xbest = np.zeros(len(isgood90))
amplbest = np.zeros(len(isgood90))
for i in np.arange(len(isgood90)):
#for i in np.arange(2):
    j = isgood90[i]
#    thismap = ndimage.gaussian_filter(maps[j,0,xmin:xmax,xmin:xmax],2)
    thismap = maps[j,0,xmin:xmax,xmin:xmax]
    if i in isgood150:
#        thismap += ndimage.gaussian_filter(maps[j,1,xmin:xmax,xmin:xmax],3)*6.
        thismap += ndimage.gaussian_filter(maps[j,1,xmin:xmax,xmin:xmax],1)*6.
    ampl_guess = np.max(thismap[25:75,25:75])*1.2
    thismap /= ampl_guess
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
    amplbest[i] = amplvec[qqqmin]
    ybest[i] = centervec[qqmin]
    xbest[i] = centervec[qmin]
    goodmaps[i,:,:] = thismap
#    models[i,:,:] = amplvec[qqqmin]*np.exp(-((yarr-ybest[i])**2+(xarr-xbest[i])**2)/2./1.7**2)
    models[i,:,:] = amplbest[i]*np.exp(-((yarr-ybest[i])**2+(xarr-xbest[i])**2)/2./1.4**2)
