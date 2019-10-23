import numpy as np
import copy
import numpy
import pickle
import glob
from spt3g import core,dfmux,std_processing, util,mapspectra,mapmaker,coordinateutils
from spt3g.mapspectra import *

'''
This code makes 2d fourier tranform arrays for a stack of sim maps

'''
files = glob.glob('/spt/user/javva/lowell/maps/ra0hdec-57.5/sims_maps_filtered_nonsense_p4_sasha_nointerp_npf/-15997020/*.g3') #define where your stack of sims is
pk_out = '/spt/user/javva/lowell/maps/ra0hdec-57.5/sims_maps_filtered_nonsense_p4_sasha_noi\
nterp_npf/-15997020/k_space/' #the directory you want the final pkl files to go in 

def read_in(file):
    return list(core.G3File(file))[0]

def apod_mask(fr, r=100):
        edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
        tt = copy.copy(fr['Wpol'].TT)
        return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)


print('making apod mask')
apod= pickle.load(open("/spt/user/javva/lowell/good_apodization_mask.pkl","rb"))
print('done with apod mask')


def return2dfft(t, q, u, apod_mask, ell_bins, ell_weights_2d = None,
                kspace_filt = None,  include_kspace_norm = True,
                qu = False, include_cross_terms = False):
    res = t.res
    if kspace_filt is None:
        kspace_filt = np.ones(np.shape(apod_mask))
    if qu:
        e_ft = mapspectra.basicmaputils.map_to_ft(q,apod_mask,kspace_filt)
        b_ft = mapspectra.basicmaputils.map_to_ft(u,apod_mask,kspace_filt)
    else:
        e_ft,b_ft = qu_to_eb_ft(q,u, apod_mask, kspace_filt)
    t_ft =  mapspectra.basicmaputils.map_to_ft(t,apod_mask,kspace_filt)
    return t_ft, e_ft, b_ft


for bun in [1]:
    print('Doing bundle %s'%bun)
    maps_var = {}
    for fle in files:
        print(fle)
        fr = read_in(fle)
        t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'],fr['Wpol'])
        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)
        cen_ells = numpy.linspace(10, 2500, 200) 
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
        tt_bun,qq_bun,uu_cbun = return2dfft(t,qf,uf, apod, ell_bins, qu=True)
        maps_var[fle] = {}
        maps_var[fle]['t']= tt_bun
        maps_var[fle]['q']= qq_bun
        maps_var[fle]['u'] = uu_cbun

        k = fle
        num = k.split('/')[-1].split('_')[-1].split('.')[0]
        nu_dict = {}
        nu_dict[k] = maps_var[k]
        with open(pk_out+'%s.pkl'%num,'wb') as handle:
            pickle.dump(nu_dict,handle,protocol=pickle.HIGHEST_PROTOCOL)
        print(k,num)
