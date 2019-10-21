import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters
Author : Jessica Avva, 2018
'''
'''
#Used this code originally to read in the data and categorize as lead or trail

def read_in(fle):
    return list(core.G3File(fle))[0]

maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*')
lt_dict = {}
for m in sorted(maps):
    a = read_in(m)
    lt_dict[m] = {}
    lt_dict[m]['q_mean_weight'] = np.mean([i for i in (a['Wpol'].QQ) if i> 0.])
    if a['Q'][800,2000]!=0:        
        lt_dict[m]['direction'] = 'right'
    if a['Q'][800,1000]!=0:        
        lt_dict[m]['direction'] = 'left'
    if (a['Q'][800,1000]!=0 and a['Q'][800,2000]!=0):        
        lt_dict[m]['direction'] = 'not-lead-trail'
    print(m)

    import pickle
    with open('lead_or_trail.pkl', 'wb') as handle:
        pickle.dump(lt_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''
rights = []
lefts = []

a = pickle.load(open('lead_or_trail.pkl', 'rb'))
for i in a.keys():
    if 'direction' not in a[i].keys():
        continue
    if a[i]['direction'] == 'right':
        rights = np.append(rights,i)
    if a[i]['direction'] == 'left':
        lefts = np.append(lefts, i)
            



coaddweights = None
coadd = {}

#some silly bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []
bads = []

print('we are going to add %s maps'%len(rights)*2) #sanity check, make sure this isn't zero

#bundle the maps chronologically with ~equal Q weights
for i,fle in enumerate(sorted(rights)):
    fr = list(core.G3File(fle))[0]
    print('adding in '+ fle)
    fl_nm = fle
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
    fr = list(core.G3File(sorted(lefts)[i]))[0]
    print('adding in '+ fle)
    fl_nm = fle
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
    core.G3Writer( '/spt/user/javva/lowell/maps/paired_lt/lt_pair_%s.g3'%(bun_num))(f)
    print('Read out /spt/user/javva/lowell/maps/paired_lt/lt_pair_%s.g3'%(bun_num))
    num_files = 0 
    defs_dict[str(bun_num)+'_wts'] = mean_weight
    bun_num = bun_num + 1
    coaddweights = None
    coadd = {}
    fil_list = []

#write out the definitions of what maps went into each bundle 
import pickle
with open('/spt/user/javva/lowell/maps/paired_lt/lead_trail_defs_dict.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
