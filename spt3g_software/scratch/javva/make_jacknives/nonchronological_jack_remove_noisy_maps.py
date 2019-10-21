import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters

'''
name = 'non_chron_weights_remove_noisy_maps'

map_dir = '/spt/user/nwhitehorn/sptpol/lowell/maps/'

base_maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*')

dirslist_pre = [d.split('/')[-1] for d in base_maps]

dirslist = [d.split('.')[0] for d in dirslist_pre]

g3_bundles = {}

np.random.shuffle(dirslist)

g3_bundles['sorted'] = dirslist



#get rid of bad maps

tf = open('/spt/user/nwhitehorn/sptpol/lowell/scripts/cleanedmapnoise.txt','r')
lines = tf.read().split('\n') 
lines2 = [d.split('\t') for d in lines] 

mn_dict = {}
for i in lines2[1:-1]:
    mn_dict[i[0]] = {}
    mn_dict[i[0]]['Q100'] = float(i[3])
    mn_dict[i[0]]['Q2000'] = float(i[4])
    mn_dict[i[0]]['U100'] = float(i[5])
    mn_dict[i[0]]['U2000'] = float(i[6])
 

coaddweights = None
coadd = {}
#some bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []
bads = []

print(g3_bundles['sorted'])

for fle in g3_bundles['sorted']:
    fr = list(core.G3File(map_dir+fle+'.g3'))[0]
    fl_nm = map_dir+fle+'.g3'
    if 'Id' not in fr or fr['Id'] != 'ra0hdec-57.5-150GHz':
        continue
    if fl_nm not in list(mn_dict.keys()):
        print('not in dictionary')
        bads = np.append(bads,fl_nm)
        continue
    if mn_dict[fl_nm]['Q100']> 1700:
        print('bad map')
        bads = np.append(bads,fl_nm)
        continue
    if mn_dict[fl_nm]['U100']> 1700:
        print('bad map')
        bads = np.append(bads,fl_nm)
        continue

    if mn_dict[fl_nm]['Q2000']> 650:
        print('bad map')
        bads = np.append(bads,fl_nm)
        continue

    if mn_dict[fl_nm]['U2000']> 650:
        print('bad map')
        bads = np.append(bads,fl_nm)
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
    if mean_weight > 5000:
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
with open('bundle_defs_%s.pkl'%name, 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




