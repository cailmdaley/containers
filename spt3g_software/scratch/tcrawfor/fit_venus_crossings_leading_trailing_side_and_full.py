import numpy as np
import scipy
from scipy import ndimage, optimize
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

# Create a function which returns a Gaussian (normal) distribution.                                 
def gauss(x, *p):
    a, b, c, d = p
    y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d

    return y

obsidstr = '67936775'
g3files = glob.glob('/spt/data/bolodata/fullrate/venus-pixelraster/'+obsidstr+'/0*.g3')
g3files.sort()

## !!!
#g3files = g3files[12:14]
## !!!

calfile = '/spt/user/production/calibration/calframe/venus-pixelraster/'+obsidstr+'.g3'
cframe = core.G3File(calfile).next()
bp = cframe['BolometerProperties']
csn = cframe['CalibratorResponseSN']
names1 = np.asarray(bp.keys())
names2 = np.asarray(csn.keys())
names3 = np.intersect1d(names1,names2)
goodnames = [name for name in names2 if csn[name] > 20.]

# output dictionary
outdict = {}

# obs and fit settings
el_venus = 20.*np.pi/180.
scanspeed_onsky = 0.3*np.cos(el_venus)
arcmin_per_sample = 60.*scanspeed_onsky/152.6
nx = 300
xtemp = (np.arange(nx)-nx/2)*arcmin_per_sample

# step through files and scans, eyeball if there's a venus crossing; if so, try to fit the leading and trailing halves to separate gaussians

## !!!
#goodnames = ['2019.000','2019.05v']
## !!!

for g3file in g3files:
    print(g3file)
    f1 = core.G3File(g3file)
    for frame in f1:
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
#                print frame
#                print(notavariable)
                for name in goodnames:
#                    print(name)
#                    try:
                    dtemp = -np.asarray(frame['RawTimestreams_I'][name])

                    y0 = np.mean(dtemp[0:50])
                    y1 = np.mean(dtemp[-50:])
                    slope = np.arange(len(dtemp))/np.float(len(dtemp))*(y1-y0) + y0
                    dtemp -= slope

                    imax = np.argmax(dtemp)
                    sdtemp = np.std(dtemp[0:50])
#                    print dtemp[imax]
#                    print sdtemp
                    if np.abs(dtemp[imax])/sdtemp > 10. and imax > nx/2 and imax < len(dtemp)-nx/2: # Venus detected!
                        dtemp2 = dtemp[imax-nx/2:imax+nx/2]
                        p_initial = [dtemp[imax], 0.0, 0.6, 0.0]
                        try:
                            errs1 = np.zeros(nx) + sdtemp
                            errs1[nx/2+1:nx/2+nx/6] *= 1000.
                            popt1, pcov1 = optimize.curve_fit(gauss, xtemp, dtemp2, p0=p_initial, sigma=errs1)
                            npars = len(popt1)
                            y_fit_1 = gauss(xtemp, *popt1)
                            errs2 = np.zeros(nx) + sdtemp
                            errs2[nx/2-nx/6:nx/2] *= 1000.
                            popt2, pcov2 = optimize.curve_fit(gauss, xtemp, dtemp2, p0=p_initial, sigma=errs2)
                            y_fit_2 = gauss(xtemp, *popt2)
                            errs3 = np.zeros(nx) + sdtemp
                            popt3, pcov3 = optimize.curve_fit(gauss, xtemp, dtemp2, p0=p_initial, sigma=errs3)
                            y_fit_3 = gauss(xtemp, *popt3)
                            ra = frame['OnlineBoresightRa']
                            if name not in outdict.keys():
                                outdict[name] = {}
                                outdict[name]['band'] = np.int(bp[name].band/10.)
                                outdict[name]['left'] = {}
                                outdict[name]['left']['data'] = np.zeros(nx)
                                outdict[name]['left']['fit_lead'] = np.zeros(nx)
                                outdict[name]['left']['fit_trail'] = np.zeros(nx)
                                outdict[name]['left']['fit_full'] = np.zeros(nx)
                                outdict[name]['left']['pars_lead'] = np.zeros(npars)
                                outdict[name]['left']['pars_trail'] = np.zeros(npars)
                                outdict[name]['left']['pars_full'] = np.zeros(npars)
                                outdict[name]['right'] = {}
                                outdict[name]['right']['data'] = np.zeros(nx)
                                outdict[name]['right']['fit_lead'] = np.zeros(nx)
                                outdict[name]['right']['fit_trail'] = np.zeros(nx)
                                outdict[name]['right']['fit_full'] = np.zeros(nx)
                                outdict[name]['right']['pars_lead'] = np.zeros(npars)
                                outdict[name]['right']['pars_trail'] = np.zeros(npars)
                                outdict[name]['right']['pars_full'] = np.zeros(npars)
                            if ra[-1] > ra[0]:
                                if popt1[0] > outdict[name]['left']['pars_lead'][0]:
                                    outdict[name]['left']['data'] = dtemp2
                                    outdict[name]['left']['fit_lead'] = y_fit_1
                                    outdict[name]['left']['pars_lead'] = popt1
                                    outdict[name]['left']['fit_trail'] = y_fit_2
                                    outdict[name]['left']['pars_trail'] = popt2
                                    outdict[name]['left']['fit_full'] = y_fit_3
                                    outdict[name]['left']['pars_full'] = popt3
                            else:
                                if popt1[0] > outdict[name]['right']['pars_lead'][0]:
                                    outdict[name]['right']['data'] = dtemp2
                                    outdict[name]['right']['fit_lead'] = y_fit_1
                                    outdict[name]['right']['pars_lead'] = popt1
                                    outdict[name]['right']['fit_trail'] = y_fit_2
                                    outdict[name]['right']['pars_trail'] = popt2
                                    outdict[name]['right']['fit_full'] = y_fit_3
                                    outdict[name]['right']['pars_full'] = popt3
                        except(RuntimeError):
                            print("Fit failed for detector "+name+"\n")
#                    print(notavariable)
#                    except:
#                        pass
#                    if np.abs(dtemp[imax])/sdtemp > 100.:
#                        print(notavariable)

pickle.dump(outdict,open('/spt/user/tcrawfor/public/venus_fit_out_28feb19.pkl','w'))
