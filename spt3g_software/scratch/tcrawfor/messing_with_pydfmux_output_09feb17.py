import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

files_baseline_09feb = glob.glob('/poleanalysis/pydfmux_output/20170209/20170209_021337_take_rawdump_baseline/data/*.pkl')

files_injection_cold_09feb = glob.glob('/poleanalysis/pydfmux_output/20170209/20170209_021440_take_rawdump_4p2MHz_broadcasted_at_slot_16/data/*.pkl')

files_injection_warm_09feb = glob.glob('/poleanalysis/pydfmux_output/20170209/20170209_025637_take_rawdump/data/*.pkl')

i_injected = 10499
i_buck = 6249
i_volcano = [5300,6100]

boarddict = {}
#for i in np.arange(16):
for i in np.arange(10):
    bind = i+1
    mzdict = {}
    for j in np.arange(2):
        mzind = j+1
        mddict = {}
        for k in np.arange(4):
            mdind = k+1
            tempdict = {}
            fb = '/poleanalysis/pydfmux_output/20170209/20170209_021337_take_rawdump_baseline/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
            d1 = pickle.load(open(fb))
            tempdict['baseline'] = d1['freq_domain']['y']
            tempdict['baseline_inj'] = tempdict['baseline'][i_injected]
            tempdict['baseline_buck'] = tempdict['baseline'][i_buck]
            tempdict['baseline_vol'] = np.sqrt(np.mean(tempdict['baseline'][i_volcano[0]:i_volcano[1]]**2))
            fc = '/poleanalysis/pydfmux_output/20170209/20170209_021440_take_rawdump_4p2MHz_broadcasted_at_slot_16/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
            d2 = pickle.load(open(fc))
            tempdict['cold'] = d2['freq_domain']['y']
            tempdict['cold_inj'] = tempdict['cold'][i_injected]
            tempdict['cold_buck'] = tempdict['cold'][i_buck]
            tempdict['cold_vol'] = np.sqrt(np.mean(tempdict['cold'][i_volcano[0]:i_volcano[1]]**2))
            fw = '/poleanalysis/pydfmux_output/20170209/20170209_025637_take_rawdump/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
            d3 = pickle.load(open(fw))
            tempdict['warm'] = d3['freq_domain']['y']
            tempdict['warm_inj'] = tempdict['warm'][i_injected]
            tempdict['warm_buck'] = tempdict['warm'][i_buck]
            tempdict['warm_vol'] = np.sqrt(np.mean(tempdict['warm'][i_volcano[0]:i_volcano[1]]**2))
            mddict[mdind] = tempdict
        mzdict[mzind] = mddict
    boarddict[bind] = mzdict
    
f_mhz = d1['freq_domain']['x']/1e6

baseline_inj = []
baseline_buck = []
baseline_vol = []
cold_inj = []
cold_buck = []
cold_vol = []
warm_inj = []
warm_buck = []
warm_vol = []

for bkey in boarddict.keys():
    for mzkey in boarddict[bkey].keys():
        for mdkey in boarddict[bkey][mzkey].keys():
            dtemp = boarddict[bkey][mzkey][mdkey]
            baseline_inj.append(dtemp['baseline_inj'])
            baseline_buck.append(dtemp['baseline_buck'])
            baseline_vol.append(dtemp['baseline_vol'])
            cold_inj.append(dtemp['cold_inj'])
            cold_buck.append(dtemp['cold_buck'])
            cold_vol.append(dtemp['cold_vol'])
            warm_inj.append(dtemp['warm_inj'])
            warm_buck.append(dtemp['warm_buck'])
            warm_vol.append(dtemp['warm_vol'])

baseline_inj = np.asarray(baseline_inj)
baseline_buck = np.asarray(baseline_buck)
baseline_vol = np.asarray(baseline_vol)
cold_inj = np.asarray(cold_inj)
cold_buck = np.asarray(cold_buck)
cold_vol = np.asarray(cold_vol)
warm_inj = np.asarray(warm_inj)
warm_buck = np.asarray(warm_buck)
warm_vol = np.asarray(warm_vol)


files_baseline_08feb = glob.glob('/poleanalysis/pydfmux_output/20170208/20170208_100913_take_rawdump/data/*.pkl')

files_injection_cold_08feb = glob.glob('/poleanalysis/pydfmux_output/20170208/20170208_102409_take_rawdump/data/*.pkl')

files_injection_warm_08feb = glob.glob('/poleanalysis/pydfmux_output/20170208/20170208_025637_take_rawdump/data/*.pkl')

boarddict2 = {}
#for i in np.arange(16):
for i in np.arange(10):
#for i in np.arange(1)+15:
    bind = i+1
    mzdict = {}
    for j in np.arange(2):
#    for j in np.arange(1)+1:
        mzind = j+1
        mddict = {}
        for k in np.arange(4):
#        for k in np.arange(1)+2:
            mdind = k+1
            tempdict = {}
            try:
#                fb = '/poleanalysis/pydfmux_output/20170208/20170208_100913_take_rawdump/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
                fb = '/poleanalysis/pydfmux_output/20170208/20170208_102409_take_rawdump/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
                d1 = pickle.load(open(fb))
                tempdict['baseline'] = d1['freq_domain']['y']
                tempdict['baseline_inj'] = tempdict['baseline'][i_injected]
                tempdict['baseline_buck'] = tempdict['baseline'][i_buck]
                tempdict['baseline_vol'] = np.sqrt(np.mean(tempdict['baseline'][i_volcano[0]:i_volcano[1]]**2))
                fc = '/poleanalysis/pydfmux_output/20170208/20170208_115658_take_rawdump/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
                d2 = pickle.load(open(fc))
                tempdict['cold'] = d2['freq_domain']['y']
                tempdict['cold_inj'] = tempdict['cold'][i_injected]
                tempdict['cold_buck'] = tempdict['cold'][i_buck]
                tempdict['cold_vol'] = np.sqrt(np.mean(tempdict['cold'][i_volcano[0]:i_volcano[1]]**2))
                fw = '/poleanalysis/pydfmux_output/20170208/20170208_102409_take_rawdump/data/IceCrate_005.IceBoard_'+str(bind)+'.Mezz_'+str(mzind)+'.ReadoutModule_'+str(mdind)+'_OUTPUT.pkl'
                d3 = pickle.load(open(fw))
                tempdict['warm'] = d3['freq_domain']['y']
                tempdict['warm_inj'] = tempdict['warm'][i_injected]
                tempdict['warm_buck'] = tempdict['warm'][i_buck]
                tempdict['warm_vol'] = np.sqrt(np.mean(tempdict['warm'][i_volcano[0]:i_volcano[1]]**2))
            except:
                pass
            mddict[mdind] = tempdict
        mzdict[mzind] = mddict
    boarddict2[bind] = mzdict
    
baseline_inj1 = []
baseline_buck1 = []
baseline_vol1 = []
cold_inj1 = []
cold_buck1 = []
cold_vol1 = []
warm_inj1 = []
warm_buck1 = []
warm_vol1 = []
baseline_inj2 = []
baseline_buck2 = []
baseline_vol2 = []
cold_inj2 = []
cold_buck2 = []
cold_vol2 = []
warm_inj2 = []
warm_buck2 = []
warm_vol2 = []

for bkey in boarddict2.keys():
    for mzkey in boarddict2[bkey].keys():
        for mdkey in boarddict2[bkey][mzkey].keys():
            dtemp2 = boarddict2[bkey][mzkey][mdkey]
            if len(dtemp2.keys()) > 1:
                baseline_inj2.append(dtemp2['baseline_inj'])
                baseline_buck2.append(dtemp2['baseline_buck'])
                baseline_vol2.append(dtemp2['baseline_vol'])
                cold_inj2.append(dtemp2['cold_inj'])
                cold_buck2.append(dtemp2['cold_buck'])
                cold_vol2.append(dtemp2['cold_vol'])
                warm_inj2.append(dtemp2['warm_inj'])
                warm_buck2.append(dtemp2['warm_buck'])
                warm_vol2.append(dtemp2['warm_vol'])
                dtemp1 = boarddict[bkey][mzkey][mdkey]
                baseline_inj1.append(dtemp1['baseline_inj'])
                baseline_buck1.append(dtemp1['baseline_buck'])
                baseline_vol1.append(dtemp1['baseline_vol'])
                cold_inj1.append(dtemp1['cold_inj'])
                cold_buck1.append(dtemp1['cold_buck'])
                cold_vol1.append(dtemp1['cold_vol'])
                warm_inj1.append(dtemp1['warm_inj'])
                warm_buck1.append(dtemp1['warm_buck'])
                warm_vol1.append(dtemp1['warm_vol'])

baseline_inj1 = np.asarray(baseline_inj1)
baseline_buck1 = np.asarray(baseline_buck1)
baseline_vol1 = np.asarray(baseline_vol1)
cold_inj1 = np.asarray(cold_inj1)
cold_buck1 = np.asarray(cold_buck1)
cold_vol1 = np.asarray(cold_vol1)
warm_inj1 = np.asarray(warm_inj1)
warm_buck1 = np.asarray(warm_buck1)
warm_vol1 = np.asarray(warm_vol1)
baseline_inj2 = np.asarray(baseline_inj2)
baseline_buck2 = np.asarray(baseline_buck2)
baseline_vol2 = np.asarray(baseline_vol2)
cold_inj2 = np.asarray(cold_inj2)
cold_buck2 = np.asarray(cold_buck2)
cold_vol2 = np.asarray(cold_vol2)
warm_inj2 = np.asarray(warm_inj2)
warm_buck2 = np.asarray(warm_buck2)
warm_vol2 = np.asarray(warm_vol2)






