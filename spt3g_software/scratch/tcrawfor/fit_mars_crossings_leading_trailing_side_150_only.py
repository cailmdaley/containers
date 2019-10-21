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

obsidstr = '54868035'
g3files = glob.glob('/spt/data/bolodata/fullrate/mars-pixelraster/'+obsidstr+'/0*.g3')
g3files.sort()

# !!!
g3files = g3files[12:14]
# !!!

calfile = '/spt/user/production/calibration/calframe/mars-pixelraster/'+obsidstr+'.g3'
cframe = core.G3File(calfile).next()
bp = cframe['BolometerProperties']
csn = cframe['CalibratorResponseSN']
names1 = np.asarray(bp.keys())
names2 = np.asarray(csn.keys())
names3 = np.intersect1d(names1,names2)
gn150 = [name for name in names2 if bp[name].band==1500. and csn[name] > 20.]

# output dictionaries
pars_lead_left = {}
pars_lead_right = {}
pars_trail_left = {}
pars_trail_right = {}
datafit_left = {}
datafit_right = {}

# obs and fit settings
el_mars = 22.75*np.pi/180.
scanspeed_onsky = 0.3*np.cos(el_mars)
arcmin_per_sample = 60.*scanspeed_onsky/152.6
nx = 300
xtemp = (np.arange(nx)-nx/2)*arcmin_per_sample

# step through files and scans, eyeball if there's a mars crossing; if so, try to fit the leading and trailing halves to separate gaussians

# !!!
gn150 = ['005.12.1.2.2548','005.12.1.2.2548']
# !!!

for g3file in g3files:
    print(g3file)
    f1 = core.G3File(g3file)
    for frame in f1:
        if frame.type is core.G3FrameType.Scan:
            if 'Turnaround' not in frame:
#                print frame
                for name in gn150:
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
                    if np.abs(dtemp[imax])/sdtemp > 100. and imax > nx/2 and imax < len(dtemp)-nx/2: # Mars detected!
                        dtemp2 = dtemp[imax-nx/2:imax+nx/2]
                        p_initial = [dtemp[imax], 0.0, 0.6, 0.0]
                        try:
                            errs1 = np.zeros(nx) + sdtemp
                            errs1[nx/2+1:nx/2+nx/6] *= 1000.
                            popt1, pcov1 = optimize.curve_fit(gauss, xtemp, dtemp2, p0=p_initial, sigma=errs1)
                            y_fit_1 = gauss(xtemp, *popt1)
                            errs2 = np.zeros(nx) + sdtemp
                            errs2[nx/2-nx/6:nx/2] *= 1000.
                            popt2, pcov2 = optimize.curve_fit(gauss, xtemp, dtemp2, p0=p_initial, sigma=errs2)
                            y_fit_2 = gauss(xtemp, *popt2)
                            ra = frame['OnlineBoresightRa']
                            if ra[-1] > ra[0]:
                                if name in pars_lead_left.keys():
                                    if popt1[0] > pars_lead_left[name][0]:
                                        pars_lead_left[name] = popt1
                                        pars_trail_left[name] = popt2
                                        datafit_left[name] = np.zeros([3,nx])
                                        datafit_left[name][0,:] = dtemp2
                                        datafit_left[name][1,:] = y_fit_1
                                        datafit_left[name][2,:] = y_fit_2
                                else:
                                    pars_lead_left[name] = popt1
                                    pars_trail_left[name] = popt2
                                    datafit_left[name] = np.zeros([3,nx])
                                    datafit_left[name][0,:] = dtemp2
                                    datafit_left[name][1,:] = y_fit_1
                                    datafit_left[name][2,:] = y_fit_2
                            else:
                                if name in pars_lead_right.keys():
                                    if popt1[0] > pars_lead_right[name][0]:
                                        pars_lead_right[name] = popt1
                                        pars_trail_right[name] = popt2
                                        datafit_right[name] = np.zeros([3,nx])
                                        datafit_right[name][0,:] = dtemp2
                                        datafit_right[name][1,:] = y_fit_1
                                        datafit_right[name][2,:] = y_fit_2
                                else:
                                    pars_lead_right[name] = popt1
                                    pars_trail_right[name] = popt2
                                    datafit_right[name] = np.zeros([3,nx])
                                    datafit_right[name][0,:] = dtemp2
                                    datafit_right[name][1,:] = y_fit_1
                                    datafit_right[name][2,:] = y_fit_2
                        except(RuntimeError):
                            print("Fit failed for detector "+name+"\n")
#                    except:
#                        pass
#                    if np.abs(dtemp[imax])/sdtemp > 100.:
#                        print(notavariable)
