#!/usr/bin/env python
import numpy, sys, os
import scipy.ndimage, scipy.interpolate, scipy.optimize
from spt3g import core, mapmaker, calibration
import argparse as ap
import numpy as np 
import pickle
import matplotlib.pyplot as plt
import glob

only_main = glob.glob('/spt/user/javva/crosstalk/fit_outputs/injected_source_70792281_smset_fasttry/*')
'''
ons = []
onsm = []
offs = []
offsom = []
offsi = []
offsomi = []
'''
obs_73730640 = {}
obs_70792281 = {}
obs_71387365 = {}
obs_72053645 = {}
obs_72623667 = {}
obs_68994321 = {}
obs_73196729 = {}
obs_70228757 = {}
obs_68940727 = {}
obs_73730640 = {}
obs_70792281 = {}
obs_71387365 = {}
obs_72053645 = {}
obs_72623667 = {}
obs_73196729 = {}
obs_63180236 = {}
obs_63398418 = {}
obs_63643116 = {}
obs_64146996 = {}
obs_70792281 = {}
obs_63653698 = {}
obs_64371421 = {}
obs_65256447 = {}
obs_71387365 = {}
obs_68187179 = {}
obs_65632251 = {}
obs_66450627 = {}
obs_66071855 = {}
obs_65752840 = {}
obs_65337318 = {}
obs_65884978 = {}
obs_66451604 = {}
obs_66539587 = {}


# add in 63398418 at some point

other_obs = [obs_71387365,obs_72053645,obs_72623667,obs_73196729,obs_73730640]
oos = ['71387365', '72053645','72623667','73196729','73730640']

other_obs = [obs_71387365,obs_72053645,obs_72623667,obs_73196729,obs_73730640,obs_64146996,obs_63180236,obs_63653698,obs_64371421,obs_65256447, obs_63643116, obs_68994321, obs_70228757, obs_68940727]
oos = ['71387365', '72053645','72623667','73196729','73730640','64146996','63180236','63653698','64371421','65256447','63643116', '68994321','70228757','68940727']

other_obs = [obs_68187179,obs_65632251,obs_66450627,obs_66071855,obs_65752840,obs_65337318,obs_65884978,obs_66451604,obs_66539587]

oos = ['68187179','65632251','66450627','66071855','65752840','65337318','65884978','66451604','66539587']
import copy
for k,i in enumerate(sorted(only_main)):
    print(k)
    bolo_num = i.split('_')[-1].split('.pkl')[0]
    om1 = pickle.load(open(i,'rb'))
    obs_70792281[bolo_num] = om1['on_source_only_main'][0][0]
    for en,ood in enumerate((other_obs)):
        try:
            oo = oos[en]
            om = pickle.load(open('/spt/user/javva/crosstalk/fit_outputs/injected_source_%s_smset_fasttry/%s_injected_source_%s_smset_fasttry_%s.pkl'%(oo,oo,oo,bolo_num), 'rb'))
            ood[bolo_num] = om['on_source_only_main'][0][0]
        except:
            print('Couldnt find %s for %s'%(bolo_num, oos[en]))
    print(k)


other_d = {}
for en,o in enumerate(other_obs):
    other_d[oos[en]] = o
with open('only_main_fit_7sand6s_plustransition_other_news.pkl', 'wb') as handle:
    pickle.dump(other_d, handle, protocol=pickle.HIGHEST_PROTOCOL)
