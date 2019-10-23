import numpy as np
import pickle
import time

#bdata = pickle.load(open('bdata_focus_onebolo_04feb17.pkl'))

nbolo = (np.shape(bdata))[0]
fwhms_angle_all = np.zeros([nbolo,9])

azrate = 0.4*core.G3Units.deg
samplerate = 152.6

for j in np.arange(nbolo):
#for j in np.arange(1):

#    if ibands[j] != 150:
    if ibands[j] != 220:
        continue

    bdtemp = bdata[j,:]

    bdt2 = bdtemp[400000:424000]
    bdt2 -= np.median(bdt2)
#    sdtemp = np.std(bdt2)
    sdtemp = tctools.robust_sigma(bdt2)
    whsaturn = np.where(bdt2 < -5.*sdtemp)
    if len(whsaturn[0]) < 50:
        continue

    d1 = (703226-25806)/9
    xmins = []
    for i in np.arange(9):
        dtemp = bdtemp[d1*i:d1*(i+1)]
        xmins.append(d1*i+np.argmin(dtemp))
    fwhms_samples = np.zeros(9)
    for i in np.arange(9):
        thisbd = bdtemp[xmins[i]-100:xmins[i]+100]
        thisbd -= np.mean(thisbd[0:20])
        thisbd_interp = np.interp(np.arange(2000)/10.,np.arange(200),thisbd)
        xmin = np.argmin(thisbd_interp)
        dhm = np.abs(thisbd_interp - thisbd_interp[xmin]/2.)
        asdhm = np.argsort(dhm)
        fwhms_samples[i] = np.abs(asdhm[0] - asdhm[1])/10.
#        pylab.plot(thisbd)
#        pylab.xlim(80,120)
        
    fwhms_sec = fwhms_samples / samplerate
    fwhms_angle = fwhms_sec*azrate

    fwhms_angle_all[j,:] = fwhms_angle
