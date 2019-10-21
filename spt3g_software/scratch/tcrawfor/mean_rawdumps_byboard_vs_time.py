import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core, std_processing, gcp
import os
import glob

#dirs = glob.glob('/poleanalysis/pydfmux_output/2017*/*take_rawdump*/')
#dirs.sort()
#
#rddict = {}
#rddict['005'] = {}
#rddict['015'] = {}
#rddict['016'] = {}
#for key in rddict.keys():
#    for i in np.arange(16):
#        rddict[key][i+1] = {}
#        for j in np.arange(2):
#            rddict[key][i+1][j+1] = {}
#            for k in np.arange(4):
#                rddict[key][i+1][j+1][k+1] = {}
#            
#ndone = 0
#for thisdir in dirs:
#    ndone += 1
#    if (ndone/20)*20 == ndone: 
#        print(str(ndone)+" dirs done.\n")
#    files = glob.glob(thisdir + '/data/*.pkl')
#    for file1 in files:
#        ff1 = filter(None,file1.split('/'))
#        ff2 = filter(None,ff1[3].split('_'))
#        ff3 = filter(None,ff1[5].split('.'))
#        thisckey = filter(None,ff3[0].split('_'))[1]
#        thisbkey = np.int(filter(None,ff3[1].split('_'))[1])
#        thismzkey = np.int(filter(None,ff3[2].split('_'))[1])
#        thismdkey = np.int(filter(None,ff3[3].split('_'))[1])
#        thistime = core.G3Time(ff2[0] + '_' + ff2[1])
#        d1 = pickle.load(open(file1))
#        rddict[thisckey][thisbkey][thismzkey][thismdkey][thistime] = d1['freq_domain']['y']
#
#rddict = pickle.load(open('rddict_temp_10feb17.pkl'))
#
#for i in np.arange(len(dirs)-365)+365:
#    ndone = i
#    thisdir = dirs[i]
#    if (ndone/20)*20 == ndone: 
#        print(str(ndone)+" dirs done.\n")
#    files = glob.glob(thisdir + '/data/*.pkl')
#    for file1 in files:
#        try:
#            ff1 = filter(None,file1.split('/'))
#            ff2 = filter(None,ff1[3].split('_'))
#            ff3 = filter(None,ff1[5].split('.'))
#            thisckey = filter(None,ff3[0].split('_'))[1]
#            thisbkey = np.int(filter(None,ff3[1].split('_'))[1])
#            thismzkey = np.int(filter(None,ff3[2].split('_'))[1])
#            thismdkey = np.int(filter(None,ff3[3].split('_'))[1])
#            thistime = core.G3Time(ff2[0] + '_' + ff2[1])
#            d1 = pickle.load(open(file1))
#            rddict[thisckey][thisbkey][thismzkey][thismdkey][thistime] = d1['freq_domain']['y']
#        except:
#            pass
#
#freqs = d1['freq_domain']['x']
#
#ddict = {}
#ddict['freqs'] = freqs
#ddict['rddict'] = rddict
#pickle.dump(ddict,open('ddict_temp_10feb17.pkl','w'))

#ddict = pickle.load(open('ddict_temp_10feb17.pkl'))

rddict = ddict['rddict']

allmjds = []
for key in rddict.keys():
    for key2 in rddict[key].keys():
        for key3 in rddict[key][key2].keys():
            for key4 in rddict[key][key2][key3].keys():
                for time1 in rddict[key][key2][key3][key4].keys():
                    allmjds.append(time1.mjd)

allmjds = np.asarray(allmjds)
umjds = []
for mjd1 in allmjds:
    if len(umjds) == 0:
        umjds.append(mjd1)
    else:
        dmjd = np.min(np.abs(np.asarray(umjds) - mjd1))
        if dmjd > 1./24.:
            umjds.append(mjd1)
umjds.sort()
umjds = np.asarray(umjds)

cratemeans = np.zeros([3,len(umjds),24999])
craten = np.zeros([3,len(umjds)])
for i in np.arange(3):
    key = (rddict.keys())[i]
    for key2 in rddict[key].keys():
        for key3 in rddict[key][key2].keys():
            for key4 in rddict[key][key2][key3].keys():
                for time1 in rddict[key][key2][key3][key4].keys():
                    if len(rddict[key][key2][key3][key4][time1]) == 24999:
                        whmjd = np.argmin(np.abs(umjds-time1.mjd))
                        craten[i,whmjd] += 1
                        cratemeans[i,whmjd,:] += rddict[key][key2][key3][key4][time1]
 
for i in np.arange(3):
    for j in np.arange(len(umjds)):
        if craten[i,j] > 0.:
            cratemeans[i,j,:] /= craten[i,j] 

#rdave = np.zeros([len(dirs),24999])
#for i in np.arange(len(dirs)):
#    tfiles = glob.glob(dirs[i]+'data/IceCrate*.pkl')
#    for tfile in tfiles:
#        d1temp = pickle.load(open(tfile))
#        rdave[i,:] += d1temp['freq_domain']['y']/np.float(len(tfiles))
