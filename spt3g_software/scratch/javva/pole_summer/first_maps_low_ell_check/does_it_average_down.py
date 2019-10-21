import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util
import copy

'''
this is a stock code for making jackknives for various parameters
Author : Jessica Avva, 2018
'''
coaddweights = {}
coadd = {}

#some silly bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []
bads = []

rt = 'adam_gain_match'
all_mps = glob.glob('/spt/user/javva/maps_2019/*/%s*/7*/*.g3'%rt)
print('we are going to add this number of maps', len(all_mps))
dts = {}
for en,m in enumerate(sorted(all_mps)):
    if int(std_processing.utils.obsid_to_g3time(m.split('/')[-2]).mjd) not in dts.keys():
        dts[int(std_processing.utils.obsid_to_g3time(m.split('/')[-2]).mjd)] = []
    k = copy.copy(dts[int(std_processing.utils.obsid_to_g3time(m.split('/')[-2]).mjd)])
    dts[int(std_processing.utils.obsid_to_g3time(m.split('/')[-2]).mjd)] = np.append(k, m)
fle = all_mps[0]
fr = list(core.G3File(fle))[-6:]
for i in fr:
    coadd[i['Id']] = {}
    coaddweights[i['Id']] = {}
bd_files = []
#bundle the maps chronologically with ~equal Q weights
for d in sorted(list(dts.keys())):
    for fle in dts[d]:
        try:
            fr = list(core.G3File(fle))[-6:]
        except:
            print('had an error')
            bd_files.append(fle)
            continue
        print('adding in '+fle)
        fl_nm = fle
        num_files = num_files + 1
        fil_list = np.append(fil_list,fle)
        defs_dict[bun_num] = fil_list
        if bun_num == 20:
            print('found 20')
            continue
        if bun_num == 21:
            print('found 21')
            continue
        if bun_num == 26:
           print('found 26')
           continue
        if bun_num == 8:
            print('found 8')
            continue
        for i in fr:
            if 'Id' not in i:
                print('had an error')
                if fle not in bd_files:
                    bd_files.append(fle)
                continue
            if len(list(coadd[i['Id']].keys())) <2 :
                print('creating frames for',i['Id'])
                for p in ('T', 'Q', 'U'):
                    coadd[i['Id']][p] = i[p]
                coaddweights[i['Id']] = i['Wpol']
            else:
                for p in ('T', 'Q', 'U'):
                    coadd[i['Id']][p] += i[p]
                coaddweights[i['Id']] += i['Wpol']
            mean_weight = np.mean(np.asarray(coaddweights[i['Id']].QQ))
        print('bundle', bun_num, 'mean weight of coadd = ', mean_weight, '; number of files = ', num_files) 
    g = core.G3Writer('/big_scratch/javva/maps_2019_lowell/%s_strict_cuts_LR_bundle_%s.g3'%(rt,str(bun_num)))
    for k in coadd.keys():
        f = core.G3Frame(core.G3FrameType.Map)
        for p in ('T', 'Q', 'U'):
            f[p] = coadd[k][p]
        f['Wpol'] = coaddweights[k]
        g(f)
    print('Wrote out bundle %s'%str(bun_num))
    num_files = 0 
    defs_dict[str(bun_num)+'_wts'] = mean_weight
    bun_num = bun_num + 1
    fil_list = []
#write out the definitions of what maps went into each bundle 
import pickle
with open('/big_scratch/javva/maps_2019_lowell/%s_strict_cuts_defs_150.pkl'%rt, 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('There were %s bad files'%len(bd_files))
print(bd_files)



