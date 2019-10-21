import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters

'''
name = 'chron_weights_remove_noisy_maps_stringent_cuts_poly4_correct'

map_dir = '/big_scratch/javva/filtering_bundles/chron_weights_remove_noisy_maps_stringent_cuts_poly8_correct/'

base_maps = [map_dir+'bundles_3gpipe_'+str(i)+'.g3' for i in range(1,19)]

#base_maps = glob.glob(map_dir+'/*')

dirslist_pre = [d.split('/')[-1] for d in base_maps]


coaddweights = None
coadd = {}
#some bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1
defs_dict = {}
fil_list = []
bads = []

print('we are going to add %s maps'%len(base_maps))

for fle in base_maps:
    fr = list(core.G3File(fle))[0]
    print('adding in '+ fle+'.g3')
    fl_nm = map_dir+fle+'.g3'

    num_files = num_files + 1
    fil_list = np.append(fil_list,fle)
    defs_dict[bun_num] = fil_list

    if coaddweights is None:
        for p in ('T', 'Q', 'U'):
            coadd[p] = fr[p]
        coaddweights = fr['Wpol']
    else:
        for p in ('T', 'Q', 'U'):
            coadd[p] += fr[p]
        coaddweights += fr['Wpol']
    mean_weight = np.mean(np.asarray(coaddweights.QQ))
    print('bundle', bun_num, 'mean weight of coadd = ', mean_weight, '; number of files = ', num_files) 
 
f = core.G3Frame(core.G3FrameType.Map)
for p in ('T', 'Q', 'U'):
    f[p] = coadd[p]
f['Wpol'] = coaddweights
core.G3Writer( '/big_scratch/javva/filtering_bundles/chron_weights_remove_noisy_maps_stringent_cuts_poly8_correct/total_coadd_19maps_p8.g3')(f)
print('wrote out')
num_files = 0 
defs_dict[str(bun_num)+'_wts'] = mean_weight
bun_num = bun_num + 1
coaddweights = None
coadd = {}
fil_list = []

import pickle
with open('current_chron_weights_remove_noisy_maps_stringent_cuts_poly8_correct_kmaskbun_coadd.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




