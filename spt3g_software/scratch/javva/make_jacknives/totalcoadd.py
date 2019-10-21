import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
this is a stock code for making jackknives for various parameters

'''

name = 'chron_weights'

bundles_list = glob.glob('/big_scratch/javva/jackknives/chron_weights/*')


'''
This deals with the fact that with lead-trail 
start obsids were 10 seconds off of what 3gepoch defines them to be

'''

g3_bundles = {}
g3_bundles['sorted'] = bundles_list


#Now that we have our bundles, we need to coadd them.

k = 'sorted'

for k in list(g3_bundles.keys()):
    try:
        files_list = g3_bundles[k]
        print(k,len(files_list))
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename = files_list)
        pipe.Add(util.framecombiner.MapFrameCombiner)
        pipe.Add(core.G3Writer, filename = '/big_scratch/javva/jackknives/%s/total_coadd.g3'%name)
        pipe.Run(profile=True)
        print("bundled %s"%k)
    except:
        print("bundle %s not done" %k)

