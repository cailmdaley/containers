import numpy as np
import copy
import numpy
import pickle
import glob
from spt3g import core,dfmux,std_processing, util,mapspectra,mapmaker,coordinateutils
from spt3g.mapspectra import *

'''
This code saves 2d fourier transforms for a stack of maps
'''

#grab the maps that go into each bundle
ep = pickle.load(open("bundle_defs_chron_remove_noisy_maps_stringent_cuts_poly4_correct_pickle.pkl","rb"))

def read_in(file):
    return list(core.G3File(file))[0]

apod= pickle.load(open("/spt/user/javva/lowell/good_apodization_mask.pkl","rb"))
print('done with apod mask')

#quick function that returns a 2dfft. note, if you want to look at these plots I recommend fftshift
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
#bundle by bundle, make dictionaries that are arrays of each maps fft

for bun in range(1,35):
    print('Doing bundle %s'%bun)
    maps_var = {}
    for fle in ep[bun]:
        print(fle)
        fr = read_in('/spt/user/nwhitehorn/sptpol/lowell/maps/'+fle+'.g3')
        t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'],fr['Wpol'])
        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)
        cen_ells = numpy.linspace(10, 2500, 200) 
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
        tt_bun,qq_bun,uu_cbun = return2dfft(t,qf,uf, apod, ell_bins, qu=True)
        maps_var[fle] = {}
        maps_var[fle]['t']= tt_bun
        maps_var[fle]['q']= qq_bun
        maps_var[fle]['u'] = uu_cbun


    import pickle
    with open('/big_scratch/javva/lowell_2dpsd/corrapod_p4_k_arrs_%s.pkl'%bun, 'wb') as handle:
        pickle.dump(maps_var, handle, protocol=pickle.HIGHEST_PROTOCOL)




