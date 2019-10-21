import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

files = [
    '/poleanalysis/pydfmux_output/20170212/20170212_042311_measure_noise_othercrates_off/data/352_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_044718_measure_noise_additional_spectrum_control_filter_3MHz/data/352_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_101530_measure_noise_dan_on_pi_on_rf/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_102314_measure_noise_dan_off_pi_on_rf/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_104248_measure_noise_dan_on_pi_off/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_104643_measure_noise_dan_off_pi_off/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_110332_measure_noise_dan_on_pi_on_mezz/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170212/20170212_110734_measure_noise_dan_off_pi_on_mezz/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_010516_measure_noise_dan_on_pi_off/data/357_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_011038_measure_noise_dan_off_pi_off/data/357_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_012436_measure_noise_dan_on_pi_mezz/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_012649_measure_noise_dan_off_pi_mezz/data/356_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_014247_measure_noise_dan_on_pi_rf/data/355_BOLOS_INFO_AND_MAPPING.pkl',
    '/poleanalysis/pydfmux_output/20170213/20170213_014444_measure_noise_dan_off_pi_rf/data/355_BOLOS_INFO_AND_MAPPING.pkl']

fkeys = [
    'filter_off_0',
    'filter_at_mezz_0',
    'filter_at_rf_1',
    'filter_at_rf_1_dan_off',
    'filter_off_1',
    'filter_off_1_dan_off',
    'filter_at_mezz_1', 
    'filter_at_mezz_1_dan_off',
    'filter_off_2',
    'filter_off_2_dan_off',
    'filter_at_mezz_2',
    'filter_at_mezz_2_dan_off',
    'filter_at_rf_2',
    'filter_at_rf_2_dan_off']

keys_all = []
for file1 in files:
    dtemp = pickle.load(open(file1))
    if len(keys_all) == 0:
        keys_all = dtemp.keys()
    else:
        for key in keys_all:
            if key not in dtemp.keys():
                keys_all.remove(key)

fbias = np.asarray([dtemp[key]['frequency'] for key in keys_all])

ndict = {}
for file1, fkey1 in zip(files,fkeys):
    dtemp = pickle.load(open(file1))
    ndict[fkey1] = np.asarray([dtemp[key]['noise']['median_noise'] for key in keys_all])

## make some plots
#
###  reproducibility
##keystub1 = 'filter_at_mezz'
##keystub2 = 'filter_off'
##ilist = ['0','1','2']
##figure()
##for istr in ilist:
##    key1 = keystub1+'_'+istr
##    key2 = keystub2+'_'+istr
##    plot(fbias/1e6,ndict[key1]/ndict[key2],'o',label=key1+'/'+key2)
##ylim(0,2)
##xlabel('bias frequency [MHz]')
##ylabel('ratio of filter in to filter out')
##legend()
##savefig('reproducibility_filter_at_mezz.png')
##close('all')
##
##keystub1 = 'filter_at_rf'
##keystub2 = 'filter_off'
##ilist = ['1','2']
##figure()
##for istr in ilist:
##    key1 = keystub1+'_'+istr
##    key2 = keystub2+'_'+istr
##    plot(fbias/1e6,ndict[key1]/ndict[key2],'o',label=key1+'/'+key2)
##ylim(0,2)
##xlabel('bias frequency [MHz]')
##ylabel('ratio of filter in to filter out')
##legend()
##savefig('reproducibility_filter_at_rf.png')
##close('all')
##
###  mean filter in / filter out
##
##keystub1 = 'filter_at_rf'
##keystub2 = 'filter_off'
##ilist = ['1','2']
##figure()
##data1 = np.zeros(len(keys_all))
##data2 = np.zeros(len(keys_all))
##for istr in ilist:
##    key1 = keystub1+'_'+istr
##    key2 = keystub2+'_'+istr
##    data1 += ndict[key1]/np.float(len(ilist))
##    data2 += ndict[key2]/np.float(len(ilist))
##plot(fbias/1e6,data1/data2,'o',label=keystub1+'/'+keystub2)
##ylim(0,2)
##xlabel('bias frequency [MHz]')
##ylabel('ratio of filter in to filter out')
##legend()
##savefig('average_filter_at_rf.png')
##close('all')
##
##keystub1 = 'filter_at_mezz'
##keystub2 = 'filter_off'
##ilist = ['0','1','2']
##figure()
##data1 = np.zeros(len(keys_all))
##data2 = np.zeros(len(keys_all))
##for istr in ilist:
##    key1 = keystub1+'_'+istr
##    key2 = keystub2+'_'+istr
##    data1 += ndict[key1]/np.float(len(ilist))
##    data2 += ndict[key2]/np.float(len(ilist))
##plot(fbias/1e6,data1/data2,'o',label=keystub1+'/'+keystub2)
##ylim(0,2)
##xlabel('bias frequency [MHz]')
##ylabel('ratio of filter in to filter out')
##legend()
##savefig('average_filter_at_mezz.png')
##close('all')
##
##keystub1 = 'filter_off'
##keystub2 = 'filter_at_mezz'
##ilist = ['0','1','2']
##keystub3 = 'filter_at_rf'
##figure()
##data1 = np.zeros(len(keys_all))
##data2 = np.zeros(len(keys_all))
##data3 = np.zeros(len(keys_all))
##for istr in ilist:
##    key1 = keystub1+'_'+istr
##    key2 = keystub2+'_'+istr
##    data1 += ndict[key1]/np.float(len(ilist))
##    data2 += ndict[key2]/np.float(len(ilist))
##    if istr != '0':
##        key3 = keystub3+'_'+istr
##        data3 += ndict[key3]/2.
##plot(fbias/1e6,data2/data1,'o',label=keystub2+'/'+keystub1)
##plot(fbias/1e6,data3/data1,'o',label=keystub3+'/'+keystub1)
##ylim(0,2)
##xlabel('bias frequency [MHz]')
##ylabel('ratio of filter in to filter out')
##legend()
##savefig('average_filter_at_mezz_and_rf.png')
##close('all')
#
##  mean filter in / filter out, DAN off
#
#keystub1 = 'filter_at_rf'
#keystub2 = 'filter_off'
#ilist = ['1','2']
#figure()
#data1 = np.zeros(len(keys_all))
#data2 = np.zeros(len(keys_all))
#for istr in ilist:
#    key1 = keystub1+'_'+istr+'_dan_off'
#    key2 = keystub2+'_'+istr+'_dan_off'
#    data1 += ndict[key1]/np.float(len(ilist))
#    data2 += ndict[key2]/np.float(len(ilist))
#plot(fbias/1e6,data1/data2,'o',label=keystub1+'/'+keystub2+', DAN off')
#ylim(0,2)
#xlabel('bias frequency [MHz]')
#ylabel('ratio of filter in to filter out')
#legend()
#savefig('average_filter_at_rf_dan_off.png')
#close('all')
#
#keystub1 = 'filter_at_mezz'
#keystub2 = 'filter_off'
#ilist = ['1','2']
#figure()
#data1 = np.zeros(len(keys_all))
#data2 = np.zeros(len(keys_all))
#for istr in ilist:
#    key1 = keystub1+'_'+istr+'_dan_off'
#    key2 = keystub2+'_'+istr+'_dan_off'
#    data1 += ndict[key1]/np.float(len(ilist))
#    data2 += ndict[key2]/np.float(len(ilist))
#plot(fbias/1e6,data1/data2,'o',label=keystub1+'/'+keystub2+', DAN off')
#ylim(0,2)
#xlabel('bias frequency [MHz]')
#ylabel('ratio of filter in to filter out')
#legend()
#savefig('average_filter_at_mezz_dan_off.png')
#close('all')
#
#keystub1 = 'filter_off'
#keystub2 = 'filter_at_mezz'
#ilist = ['1','2']
#keystub3 = 'filter_at_rf'
#figure()
#data1 = np.zeros(len(keys_all))
#data2 = np.zeros(len(keys_all))
#data3 = np.zeros(len(keys_all))
#for istr in ilist:
#    key1 = keystub1+'_'+istr+'_dan_off'
#    key2 = keystub2+'_'+istr+'_dan_off'
#    data1 += ndict[key1]/np.float(len(ilist))
#    data2 += ndict[key2]/np.float(len(ilist))
#    key3 = keystub3+'_'+istr+'_dan_off'
#    data3 += ndict[key3]/np.float(len(ilist))
#plot(fbias/1e6,data2/data1,'o',label=keystub2+'/'+keystub1+', DAN off')
#plot(fbias/1e6,data3/data1,'o',label=keystub3+'/'+keystub1+', DAN off')
#ylim(0,2)
#xlabel('bias frequency [MHz]')
#ylabel('ratio of filter in to filter out')
#legend()
#savefig('average_filter_at_mezz_and_rf_dan_off.png')
##close('all')



