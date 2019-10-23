import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters
Author : Jessica Avva, 2018
'''
#define where the maps come from, and where you want the bundles to go
name = 'bundles_with_lt'

map_dir = '/spt/user/nwhitehorn/sptpol/lowell/finalmaps/'

#This file is the 1D fourier transforms of all of the maps. Basic metrics for cutting bad maps.

path = "/spt/user/nwhitehorn/sptpol/lowell/finalmaps/mapnoise"
print('reading in file')
import csv
with open(path) as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)
print('read in noise file')
g3_bundles = {}


#Make a dictionary of all of the noise values for various ell so you can do statistics for all the maps to find outliers. Note tha tnumber of points in PS is hard coded in
pts_dict = {}
for i, pt in enumerate(d[0]):
    try:
        pts_dict[pt] = [float(m[i]) for m in d[1:]]
    except:
        pts_dict[pt] = [m[i] for m in d[1:]]
#select the outliers. Note that this is a slightly motivated choice I made by looking at some histograms.
def is_outlier(data):
    q25,q75 = np.percentile(data, 25), np.percentile(data, 75)
    print(q25, q75)
    iqr = q75 - q25
    print(iqr)
    #before I had iqr * 1.5
    cut_off = iqr * 4.
    #cut_off = iqr* 1.5
    print(cut_off)
    lower, upper = q25 - cut_off, q75 + cut_off
    print(lower, upper)
    return lower, upper
outliers_dict = {}
#make the outlier cutoff for each ell bin
for k in d[0][1:]:
    print(k)
    outliers_dict[k] = is_outlier(pts_dict[k]) 
print('made a dictionary of the outliers')
bad_list = []                             
good_list = []                      
#for each map, make sure that none of the ell bins are outliers. If so mark the map in the bad maps lis
for fl in d[1:]:
    for i, k in enumerate(d[0]):
        if k == 'File':
            continue
        if fl[0] in bad_list:
            continue
        if float(fl[i])> outliers_dict[k][1]:
            #print(fl[i],'bigger than',outliers_dict[k][1])
            bad_list = np.append(bad_list, fl[0])
            continue
        if float(fl[i])< 200.0:
            #print(fl[i],'smaller than',200.0)
            bad_list = np.append(bad_list, fl[0])
    if fl[0] not in bad_list:
        good_list = np.append(good_list, fl[0])
print('This cut got rid of %s percent of the data'%(len(bad_list)/(len(good_list)+len(bad_list))))
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
    fr = list(core.G3File(map_dir+fle))[0]
    print('adding in '+ map_dir+fle)
    fl_nm = map_dir+fle
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
core.G3Writer( '/big_scratch/javva/all_data_bundles/coadd.g3'%(name,bun_num))(f)
num_files = 0 
defs_dict[str(bun_num)+'_wts'] = mean_weight
bun_num = bun_num + 1
coaddweights = None
coadd = {}
fil_list = []

#write out the definitions of what maps went into each bundle 
import pickle
with open('/big_scratch/javva/all_data_bundles/coadd_defs.pkl', 'wb') as handle:
    pickle.dump(defs_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




