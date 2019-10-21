from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import numpy as np
import os
import pylab
import scipy.ndimage as ndimage
try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError
import glob

#defines the input and output sims
inputs_list = glob.glob('/spt/user/javva/lowell/maps/ra0hdec-57.5/sim_input_maps_sasha_noninterp_lpf/-15997020/k_space/*') 
outputs_dir = '/spt/user/javva/lowell/maps/ra0hdec-57.5/sims_maps_filtered_nonsense_p4_sasha_nointerp_npf/-15997020/k_space/'

ep = {}

for i in inputs_list:
        fle_out =outputs_dir+i.split('/')[-1]
        if not os.path.isfile(fle_out):
                continue
        print(i)
        ep[i] = {}
        fi = pickle.load(open(i,'rb'))
        fo = pickle.load(open(fle_out,'rb'))
        ko = list(fo.keys())[0]
        ki = list(fi.keys())[0]
        #make a stack of (output/input)**2
        ep[i]['q_ft_sn'] = np.abs((fo[ko]['q']/fi[ki]['q'])**2)
        ep[i]['u_ft_sn'] = np.abs((fo[ko]['u']/fi[ki]['u'])**2)
        ep[i]['t_ft_sn'] = np.abs((fo[ko]['t']/fi[ki]['t'])**2)
#copy these just to get the right shape
eq = copy.copy(ep[list(ep.keys())[0]]['q_ft_sn'])                           
eu = copy.copy(ep[list(ep.keys())[0]]['u_ft_sn'])                           
et = copy.copy(ep[list(ep.keys())[0]]['t_ft_sn'])                           

print('doing the averages')
#for each spot in k-space, take the median of all the values. despite being noiseless sims, there can be divide by zero errors that could throw off the mean                                                              
for m in range(0,1500):                                               
    print(m)                                                          
    for i in range(0, 2700):                                          
            q_list = [ep[map]['q_ft_sn'][m,i] for map in list(ep.keys())]       
            t_list = [ep[map]['t_ft_sn'][m,i] for map in list(ep.keys())]       
            u_list = [ep[map]['u_ft_sn'][m,i] for map in list(ep.keys())]       
            eq[m,i]=1./np.median(q_list)                                     
            et[m,i]=1./np.median(t_list)                                     
            eu[m,i]=1./np.median(u_list)                                     

dict_e = {}
dict_e['q'] = eq
dict_e['u'] = eu
dict_e['t'] = et

import pickle                                                               
with open('simkspacetf_100_correctapod_median_pluscor.pkl', 'wb') as handle:
        pickle.dump(dict_e, handle, protocol=pickle.HIGHEST_PROTOCOL)

