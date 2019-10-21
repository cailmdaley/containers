import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters

'''
name = 'chron_weights_remove_noisy_maps'

map_dir = '/spt/user/nwhitehorn/sptpol/lowell/maps/'

base_maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*')

dirslist_pre = [d.split('/')[-1] for d in base_maps]

dirslist = [d.split('.')[0] for d in dirslist_pre]

g3_bundles = {}

g3_bundles['sorted'] = sorted(dirslist)

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

norms = pickle.load(open("normalization_files_final.pkl","rb")) 

coaddweights = None
coadd = {}
#some bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []
bads = []
for fle in g3_bundles['sorted'][90:100]:
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
    if np.mean(np.asarray(fr['Wpol'].TT)) == 0.0:
        print('bad map (zeros)')
        continue

    num_files = num_files + 1
    fil_list = np.append(fil_list,fle)
    defs_dict[bun_num] = fil_list
    #calculate the renormalization constants
    nc = np.mean([norms[map_dir+fle+'.g3']['Q'],norms[map_dir+fle+'.g3']['U']])
    if np.isfinite(nc) == False:
        print('Infinite Integral of Power Spec')
        continue
    #calculate the means of every piece
    const = {}
    const['T'] = nc/np.mean(np.asarray(fr['Wpol'].QQ))
    const['Q'] = nc/np.mean(np.asarray(fr['Wpol'].QQ)) 
    const['U'] = nc/np.mean(np.asarray(fr['Wpol'].QQ))

    new_weights = core.G3SkyMapWeights()
    new_weights.weight_type = core.WeightType.Wpol
    new_weights.TT = (fr['Wpol'].TT/np.mean(np.asarray(fr['Wpol'].QQ)))*nc
    new_weights.TU = (fr['Wpol'].TU/np.mean(np.asarray(fr['Wpol'].QQ)))*nc
    new_weights.TQ = (fr['Wpol'].TQ/np.mean(np.asarray(fr['Wpol'].QQ)))*nc
    new_weights.UU = (fr['Wpol'].UU/np.mean(np.asarray(fr['Wpol'].QQ)))*nc
    new_weights.QU = (fr['Wpol'].QU/np.mean(np.asarray(fr['Wpol'].QQ)))*nc
    new_weights.QQ = (fr['Wpol'].QQ/np.mean(np.asarray(fr['Wpol'].QQ)))*nc

    if coaddweights is None:
        for p in ('T', 'Q', 'U'):
            coadd[p] = fr[p]*const[p]
        coaddweights = new_weights
    else:
        for p in ('T', 'Q', 'U'):
            coadd[p] += fr[p]*const[p]
        coaddweights += new_weights
    mean_weight = np.mean(np.asarray(coaddweights.QQ))
    print('bundle', bun_num, 'mean weight of coadd = ', mean_weight, '; number of files = ', num_files) 
    #if mean_weight > 5000:

f = core.G3Frame(core.G3FrameType.Map)
for p in ('T', 'Q', 'U'):
    f[p] = coadd[p]
f['Wpol'] = coaddweights
core.G3Writer( '/big_scratch/javva/reweight_coadds/100_small_bundle_reweighted.g3')(f)
num_files = 0 
defs_dict[str(bun_num)+'_wts'] = mean_weight
bun_num = bun_num + 1
coaddweights = None
coadd = {}
fil_list = []

import pickle
with open('bundle_defs_100_small_map_coadd.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


