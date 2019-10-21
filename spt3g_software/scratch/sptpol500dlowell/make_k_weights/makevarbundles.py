from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import numpy as np
import os
import pylab
import scipy.ndimage as ndimage
import glob

'''
This code makes inverse variance weights maps for bundles of 2-d fourier transform maps
Jessica Avva 2018
'''
#This is doing a silly copy thing to allocate in memory an object of the right shape you can overwrite. 
for bn in range(1,35):
        ep = pickle.load(open("/spt/user/javva/lowell/2dpsd/p4_k_arrs_%s.pkl"%bn,"rb"))
        eq = copy.copy(ep[list(ep.keys())[0]]['q'])                           
        eu = copy.copy(ep[list(ep.keys())[0]]['u'])                           
        et = copy.copy(ep[list(ep.keys())[0]]['t'])                           
        for m in range(0,1500):                                               
                print(m)
                for i in range(0, 2700):                                       
                        q_list = [ep[map]['q'][m,i] for map in list(ep.keys())]       
                        t_list = [ep[map]['t'][m,i] for map in list(ep.keys())]       
                        u_list = [ep[map]['u'][m,i] for map in list(ep.keys())]       
                        eq[m,i]=1./np.var(q_list)                                     
                        et[m,i]=1./np.var(t_list)                                     
                        eu[m,i]=1./np.var(u_list)                                     
        with open('/spt/user/javva/lowell/2dpsd/bun%s.pkl'%bn, 'wb') as handle:
                pickle.dump(dict_e, handle, protocol=pickle.HIGHEST_PROTOCOL)
