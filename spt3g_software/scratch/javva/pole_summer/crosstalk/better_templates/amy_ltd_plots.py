#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import argparse as ap
import numpy as np 
import pickle
import matplotlib.pyplot as plt
import glob
from matplotlib.colors import LogNorm

#Import some backround info for freuqency cuts, etc.

carrier_freqs = pickle.load(open('/spt/user/javva/crosstalk/templates/63643116.pkl','rb'))['carrier_freqs']                                                                    
wiring_info = pickle.load(open('/spt/user/javva/crosstalk/templates/63643116.pkl','rb'))['wiring_info']
params = pickle.load(open('/spt/user/javva/crosstalk/templates/63643116.pkl','rb'))['params']
params_amps = params['official_params'] 
pix_info = params_amps
#import the crosstalk fit information for 1000 bolometers

a = pickle.load(open('thousand_bolos_6_and_7_withtransition.pkl','rb'))
b = pickle.load(open('thousand_bolos_6_and_7_withtransition_only_news.pkl','rb'))

#import the main bolometer rcw38 value for all the bolometers in the focal plane (for normalization)
mba = pickle.load(open('only_main_fit_7sand6s_plustransition.pkl','rb'))


#break the indivitual observations out into dictionaries
obs_63180236 = a['63180236']
obs_63643116 = a['63643116']
obs_64146996 = a['64146996']
obs_70792281 = a['70792281']
obs_63653698 = a['63653698']
obs_64371421 = a['64371421']
obs_65256447 = a['65256447']
obs_71387365 = a['71387365']
obs_72053645 = a['72053645']
obs_72623667 = a['72623667']
obs_73196729= a['73196729']
obs_73730640= a['73730640']
obs_68940727 = a['68940727']
obs_68994321 = a['68994321']
obs_70228757 = a['70228757']

obs_68187179 = b['68187179']
obs_65632251 = b['65632251']
obs_66450627 = b['66450627']
obs_66071855 = b['66071855']
obs_65752840 = b['65752840']
obs_65337318 = b['65337318']
obs_65884978 = b['65884978']
obs_66451604 = b['66451604']
obs_66539587 = b['66539587']


#some options for plotting other subsets of the data
'''

#middle ones
other_obs = [obs_68187179,obs_65632251,obs_66450627,obs_66071855,obs_65752840,obs_65337318,obs_65884978,obs_66451604,obs_66539587]

oosmo = ['68187179','66450627','66071855','65752840','65337318','65884978','66451604','66539587']


#6s

other_obs = [obs_64146996,obs_63180236,obs_63653698,obs_64371421,obs_65256447, obs_63643116, obs_68940727, obs_68994321]
oos = ['64146996','63180236','63653698','64371421','65256447','63643116', '68940727', '68994321']


#7s

other_obs = [obs_71387365,obs_72053645,obs_72623667,obs_73196729,obs_73730640, obs_70792281, obs_70228757]
oos = ['71387365', '72053645','72623667','73196729','73730640','70792281','70228757']


#both
'''

other_obs = [obs_64146996,obs_63180236,obs_63653698,obs_64371421,obs_65256447, obs_63643116, obs_70792281,obs_71387365,obs_72053645,obs_72623667,obs_73196729,obs_73730640, obs_68940727, obs_68994321, obs_70228757]
oos = ['64146996','63180236','63653698','64371421','65256447','63643116','70792281','71387365', '72053645','72623667','73196729','73730640','68940727','68994321','70228757']

#put all of the main bolometers into one mba2 

oosmo = ['68187179','65632251','66450627','66071855','65752840','65337318','65884978','66451604','66539587']
mba2 = pickle.load(open('only_main_fit_7sand6s_plustransition_other_news.pkl','rb'))

file = pickle.load(open('only_main_fit_7sand6s_plustransition_other_news.pkl','rb'))

for i in oosmo:
    if i not in mba.keys():
        mba[i] = mba2[i]

oos = sorted(list(mba.keys()))

#remove bad points

oos = ['63180236',
 '63643116',
 '63653698',
 '64146996',
 '64371421',
 '65256447',
 '65337318',
 '70228757',
 '70792281',
 '71387365',
 '72053645',
 '72623667',
 '73196729',
 '73730640']


other_obs = [obs_63180236,obs_63643116, obs_63653698, obs_64146996, obs_64371421,obs_65256447, obs_65337318, obs_70228757, obs_70792281, obs_71387365, obs_72053645, obs_72623667, obs_73196729,obs_73730640]


#Distance Cut Plot

def find_distance(boloa,bolob):
    mbx1 = (360./2)+params_amps[boloa]['x_offset']
    mby1 = (360./2)+params_amps[boloa]['y_offset']
    mbx2 = (360./2)+params_amps[bolob]['x_offset']
    mby2 = (360./2)+params_amps[bolob]['y_offset']
    return np.sqrt((mbx1-mbx2)**2 +(mby1-mby2)**2)

for o in other_obs:
    1+1

#set the distance cut
lim = 4

plt.figure()

#make the list of means 

means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            continue
        cst = []
        chb = None
        doplt = 0
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['on_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if np.abs(np.mean(cst)) > 18 and np.std(cst) > 10:
                continue
            means.append((np.mean(cst)))

print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Signal (with cut)', histtype ='step')

scut = np.nan_to_num(means)

#do the same for noise,and S, N w/o cuts

means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        #if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            #continue
        cst = []
        chb = None
        doplt = 0
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['on_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if True:
                means.append((np.mean(cst)))

print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Signal (with NO cut)', histtype ='step')

means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            continue
        cst = []
        chb = None
        doplt = 0
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['off_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if np.abs(np.mean(cst)) > 18 and np.std(cst) > 10:
                continue
            means.append((np.mean(cst)))


print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Noise (with cut)', histtype ='step')


means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
#        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
#            continue
        cst = []
        chb = None
        doplt = 0
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['off_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if True:
                means.append((np.mean(cst)))


print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Noise (with NO cut)', histtype ='step')

plt.legend()
plt.xlim(-20,20)
plt.title('Singal and Noise Histogram for 1000 bolometers')
plt.xlabel('Percentage Crosstalk [%]')
plt.ylabel('Number of Bolometer Coefficients')
plt.ylim(0,200)

#make xtalk vs freq sep plot

plt.figure()
means = []
varis = []
cs = []
freqs = []
ct = 0
for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            continue
        cst = []
        chb = None
        doplt = 0
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['on_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if np.abs(np.mean(cst)) > 18 and np.std(cst) > 10:
                ct += 1
                continue
            means.append((np.mean(cst)))
            freqs.append(np.abs(mcf-(carrier_freqs[o[b]['other_bolos_order'][i]]/core.G3Units.MHz)))
print('cut %s outliers'%ct)
a = plt.hist2d(means,freqs, bins = [np.linspace(-100,100,6200),np.linspace(0,1,200)],norm =LogNorm())
plt.clf()
bins = a[2]
xta = means
fhist = {}
for i in bins:
    fhist[i] = []
for i,f in enumerate(freqs):
    b = min(bins, key = lambda x:abs(x-f))
    fhist[b] = np.append(fhist[b], xta[i])
for i in fhist.keys():
    plt.errorbar(i, np.median(fhist[i]),np.std(fhist[i])/np.sqrt(len(fhist[i])),marker = '.', color = 'k')
plt.plot([0,.3],[0,0],'k')
plt.xlim(0,.2)
plt.ylim(-2,2)
plt.ylabel('Median Percentage Crosstalk')
plt.xlabel('Frequency Separation (MHz)')
plt.plot([0.029,0.029], [-3,3], 'r--',label = '29 kHz Spacing LC Board Design Spec')
plt.title('Frequency Separation vs. Median crosstalk (With Distance and std Cut)')
plt.legend()


plt.figure()
plt.hist(np.abs(np.nan_to_num(scut)), bins = np.linspace(0,100,1000), cumulative = True, normed = True)

plt.legend()
plt.xlim(0,20)
plt.title('Cumumulative Crosstalk Histogram for 1000 bolometers')
plt.xlabel('Percentage Crosstalk [%]')
plt.ylabel('Fraction of Bolometer Coefficients (cumulative)')
#plt.ylim(0,200)

#plot with weird close detectors clut

plt.figure()

means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            continue
        p = b
        doplt = 0
        if np.isfinite(np.sqrt(pix_info[p]['x_offset']**2+pix_info[p]['y_offset']**2)):
            if 2*(np.sqrt(pix_info[p]['x_offset']**2+pix_info[p]['y_offset']**2)) <30:
                doplt = 1
        cst = []
        chb = None
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['on_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if np.abs(np.mean(cst)) > 18 and np.std(cst) > 10:
                continue
            means.append((np.mean(cst)))

print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Signal (with cut)', histtype ='step')

means = []
varis = []
cs = []

for endd, b in enumerate(other_obs[0].keys()):
    print(endd)
    if b == 'summary':
        continue
    nope = 0
    for lz in other_obs:
        if b not in lz:
            nope = 1
    if nope ==1 :
        continue
    for i in range(len(obs_70792281[b]['other_bolos_order'])):
        mcf = carrier_freqs[b]/core.G3Units.MHz 
        if find_distance(b,o[b]['other_bolos_order'][i])< lim:
            continue
        p = b
        doplt = 0
        if np.isfinite(np.sqrt(pix_info[p]['x_offset']**2+pix_info[p]['y_offset']**2)):
            if 2*(np.sqrt(pix_info[p]['x_offset']**2+pix_info[p]['y_offset']**2)) <30:
                doplt = 1
        cst = []
        chb = None
        for k, o in enumerate(other_obs):
            if chb != None:
                if chb != o[b]['other_bolos_order'][i]:
                    doplt = 1
            if chb == None:
                chb = o[b]['other_bolos_order'][i]
            if o[b]['other_bolos_order'][i] in mba[oos[k]].keys():
                cst.append(100*o[b]['off_source'][0][i]/mba[oos[k]][o[b]['other_bolos_order'][i]])
        if doplt == 0:
            if np.abs(np.mean(cst)) > 18 and np.std(cst) > 10:
                continue
            means.append((np.mean(cst)))

print('plotting, minus %s Nans'%len(np.where(np.isfinite(means) == False)[0])
)
plt.hist(np.nan_to_num(means), bins = np.linspace(-500,100,4200), label = 'Noise (with cut)', histtype ='step')
plt.legend()
plt.xlim(-20,20)
plt.title('Singal and Noise Histogram for 1000 bolometers - cut strange overlapping ones')
plt.xlabel('Percentage Crosstalk [%]')
plt.ylabel('Number of Bolometer Coefficients')
plt.ylim(0,200)


plt.show()
