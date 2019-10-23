import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters
Author : Jessica Avva, 2018
'''
#define where the maps come from, and where you want the bundles to go
name = 'chron_weights_remove_noisy_maps_stringent_cuts_poly4_correct'

map_dir = '/spt/user/nwhitehorn/sptpol/lowell/maps/'

base_maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*')

dirslist_pre = [d.split('/')[-1] for d in base_maps]

#This file is the 1D fourier transforms of all of the maps. Basic metrics for cutting bad maps.
ep = pickle.load(open("/spt/user/nwhitehorn/sptpol/lowell/scripts/mapnoise.pkl","rb"))

g3_bundles = {}

file_list = sorted(glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*'))
ids = [f0.split('/')[-1].split('.')[0] for f0 in file_list]
real_ids = [i for i in ids if glob.glob('/spt/user/javva/lowell/maps/ra0hdec-57.5/maps_renormalized_gain/'+i+'/*') != []] #same syntax as g3_bundles                                                                                          
maps_to_bun = ['/spt/user/javva/lowell/maps/ra0hdec-57.5/maps_renormalized_gain/'+i+'/maps_renormalized_\gain_'+i+'.g3' for i in real_ids]


#get rid of bad maps

dirslist1 = [d.split('.')[0] for d in dirslist_pre]
print(len(dirslist1))

dirslist = [d for d in dirslist1 if d in real_ids]
g3_bundles = {}
print(len(dirslist))
g3_bundles['sorted'] = sorted(dirslist)


#Make a dictionary of all of the noise values for various ell so you can do statistics for all the maps to find outliers. Note tha tnumber of points in PS is hard coded in
pts_dict = {}
for q in range(80):
    q_list = [ep['Q'][n][q] for n in list(ep['Q'].keys())]
    pts_dict[q] = q_list

#select the outliers. Note that this is a slightly motivated choice I made by looking at some histograms.
def is_outlier(data):
    q25,q75 = np.percentile(data, 25), np.percentile(data, 75)
    iqr = q75 - q25
    cut_off = iqr * 1.5
    lower, upper = q25 - cut_off, q75 + cut_off
    return lower, upper
outliers_dict = {}
#make the outlier cutoff for each ell bin
for k in range(80):
    outliers_dict[k] = is_outlier(pts_dict[k]) 

bad_list = []                             
good_list = []                      
#for each map, make sure that none of the ell bins are outliers. If so mark the map in the bad maps list
for fl in list(ep['Q'].keys()):
    pts = ep['Q'][fl]
    for k in range(80):
         if fl in bad_list:
             continue
         if pts[k]> outliers_dict[k][1]:
             print(ep['bins'][k], pts[k],outliers_dict[k][1])
             bad_list = np.append(bad_list, fl)
             continue
    if fl not in bad_list:
        good_list = np.append(good_list, fl)

coaddweights = None
coadd = {}

#some silly bookkeeping stuff for matching bundles later on
num_files = 0
bun_num = 1 
defs_dict = {}
fil_list = []
bads = []

print('we are going to add %s maps'%len(good_list)) #sanity check, make sure this isn't zero

#bundle the maps chronologically with ~equal Q weights
for fle in sorted(good_list):
    fr = list(core.G3File(map_dir+fle+'.g3'))[0]
    print('adding in '+ map_dir+fle+'.g3')
    fl_nm = map_dir+fle+'.g3'
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
    if mean_weight > 3500: #note 3500 is an arbrtrary cutoff gotten by eyeballing the coadd weight
        f = core.G3Frame(core.G3FrameType.Map)
        for p in ('T', 'Q', 'U'):
            f[p] = coadd[p]
        f['Wpol'] = coaddweights
        core.G3Writer( '/big_scratch/javva/filtering_bundles/%s/bundles_3gpipe_%s.g3'%(name,bun_num))(f)
        num_files = 0 
        defs_dict[str(bun_num)+'_wts'] = mean_weight
        bun_num = bun_num + 1
        coaddweights = None
        coadd = {}
        fil_list = []

#write out the definitions of what maps went into each bundle 
import pickle
with open('bundle_defs_chron_remove_noisy_maps_stringent_cuts_poly4_correct_pickle.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




