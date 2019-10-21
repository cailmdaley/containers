import numpy as np
import pickle
import glob
from spt3g import core,dfmux,std_processing, util

'''
This code takes in the bundles from JT's 500dBB analysis
and bundles the lowell 3g analysis bundles into this. 

'''
defs = pickle.load(open('BB_bundle_defs_150GHz.pkl','rb'))

dirs = glob.glob('/spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/*')

dirslist = [d.split('/')[-1] for d in dirs]

'''
This deals with the fact that with lead-trail 
start obsids were 10 seconds off of what 3gepoch defines them to be

'''

dirslist_cal = [int(m) +10 for m in dirslist] 

obsdates = [str(std_processing.obsid_to_g3time(d)) for d in dirslist_cal]

jtformat = [m.split(':')[0].split('-')[2]+m.split(':')[0].split('-')[1]
            +m.split(':')[0].split('-')[0]+'_'+m.split(':')[1]
            +m.split(':')[2]+m.split(':')[3][0:2] for m in obsdates]

jtf = {}

for k,i in enumerate(jtformat):
    if 'Jan' in i:
        m = i.replace("Jan","01")
    if 'Feb' in i:
        m = i.replace("Feb","02")
    if 'Mar' in i:
        m = i.replace("Mar","03")
    if 'Apr' in i:
        m = i.replace("Apr","04")
    if 'May' in i:
        m = i.replace("May","05")
    if 'Jun' in i:
        m = i = i.replace("Jun","06")
    if 'Jul' in i:
        m = i.replace("Jul","07")
    if 'Aug' in i:
        m = i.replace("Aug","08")
    if 'Sep' in i:
        m = i.replace("Sep","09")
    if 'Oct' in i:
        m = i.replace("Oct","10")
    if 'Nov' in i:
        m = i.replace("Nov","11")
    if 'Dec' in i:
        m = i.replace("Dec","12")
    jtf[m] = dirslist[k]

d = defs

g3_bundles = {}

for bun_num in d.keys():
    try:
        g3_bundles[bun_num] = [jtf[m] for m in d[bun_num]]
        print len(g3_bundles[bun_num]), len(d[bun_num])
        print "did one %s" %bun_num
    except:
        continue


#Now that we have our bundles, we need to coadd them.

for k in g3_bundles.keys():
    try:
        files_list = g3_bundles[k]
        dirslist_cal = [int(m) -10 for m in files_list] 
        files_list_g3 = ['/spt/user/nwhitehorn/sptpol/lowell/maps/'+i+'.g3' for i in files_list]
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename = files_list_g3)
        pipe.Add(util.framecombiner.MapFrameCombiner)
        pipe.Add(core.G3Writer, filename = '/big_scratch/javva/bundles3g/bundles_3gpipe_0%s.g3'%k)
        pipe.Run(profile=True)
        print "bundled %s"%k
    except:
        print "bundle %s not done" %k

