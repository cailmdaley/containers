import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters

'''
name = 'chron_weights_uncleanedmaps'

map_dir = '/spt/user/nwhitehorn/sptpol/lowell/uncleanedmaps/'

base_maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/uncleanedmaps/*')

dirslist_pre = [d.split('/')[-1] for d in base_maps]

dirslist = [d.split('.')[0] for d in dirslist_pre]

g3_bundles = {}

a = pickle.load(open('bundle_defs_chron_even_weights.pkl','rb'))

g3_bundles['sorted'] = sorted(dirslist)

coaddweights = None
coadd = {}
#some bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []


for l in np.arange(1,30):
    fl = a[l]
    for fle in fl:
        print(fle)
        fr = list(core.G3File(map_dir+fle+'.g3'))[0]
        if 'Id' not in fr or fr['Id'] != 'ra0hdec-57.5-150GHz':
            continue
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
    core.G3Writer( '/big_scratch/javva/jackknives/%s/bundles_3gpipe_%s.g3'%(name,bun_num))(f)
    num_files = 0 
    defs_dict[str(bun_num)+'_wts'] = mean_weight
    bun_num = bun_num + 1
    coaddweights = None
    coadd = {}
    fil_list = []

import pickle
with open('bundle_defs_uncleaned.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




