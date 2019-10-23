import numpy as np
import scipy
from scipy import ndimage
import pickle
from spt3g import core
import os
import glob

dirs = glob.glob('/scratch/pydfmux_output/2017*/*take_rawdump*/')
dirs.sort()

rddict = {}
rddict['005'] = {}
rddict['015'] = {}
rddict['016'] = {}
            
ndone = 0
for thisdir in dirs:
    ndone += 1
    if (ndone/20)*20 == ndone: 
        print(str(ndone)+" dirs done.\n")
    files = glob.glob(thisdir + '/data/*.pkl')
    for file1 in files:
#        try:
        if 'IceCrate' in file1:
            ff1 = filter(None,file1.split('/'))
            ff2 = filter(None,ff1[3].split('_'))
            ff3 = filter(None,ff1[5].split('.'))
            thisckey = filter(None,ff3[0].split('_'))[1]
            thisbkey = np.int(filter(None,ff3[1].split('_'))[1])
            thismzkey = np.int(filter(None,ff3[2].split('_'))[1])
            thismdkey = np.int(filter(None,ff3[3].split('_'))[1])
            if thisbkey == 4 and thismzkey == 2 and thismdkey == 2:
                try:
                    thistime = core.G3Time(ff2[0] + '_' + ff2[1])
                    d1 = pickle.load(open(file1))
                    if len(d1['freq_domain']['y']) == 24999: 
                        rddict[thisckey][thistime] = d1['freq_domain']['y']            
                except:
                    pass

#        except:
#            pass
        else:
            ff1 = filter(None,file1.split('/'))
            ff2 = filter(None,ff1[3].split('_'))
            ff3 = filter(None,ff1[5].split('.'))
            thisckey = '016'
            thisbkey = 4
            thismzkey = np.int(filter(None,ff3[1].split('_'))[1])
            thismdkey = np.int(filter(None,ff3[2].split('_'))[1])
            if thismzkey == 2 and thismdkey == 2:
                try:
                    thistime = core.G3Time(ff2[0] + '_' + ff2[1])
                    d1 = pickle.load(open(file1))
                    if len(d1['freq_domain']['y']) == 24999: 
                        rddict[thisckey][thistime] = d1['freq_domain']['y']
                except:
                    pass

freqs = d1['freq_domain']['x']

#ddict = {}
#ddict['freqs'] = freqs
#ddict['rddict'] = rddict
#ickle.dump(ddict,open('ddict_temp_10feb17_fewboards.pkl','w'))

