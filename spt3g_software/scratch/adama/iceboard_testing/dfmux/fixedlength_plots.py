import numpy as np
import matplotlib.pyplot as plt
import dfmuxtools
import cPickle as pickle
import sptpol_arcy_bolos
import os
import hwm3channels


filestr = 'dfmux_idf_20151213_095735_to_20151213_120041.g3'
hwm_num = 3
os.mkdir('figures_%s' % filestr)
scans = dfmuxtools.extract_fixed_scans('/data/sptdaq/iceboard_testing/dfmux/processed/' + filestr, 10)
nmodules = 4

if hwm_num == 3:
    chaninfo = hwm3channels.info()
elif hwm_num == 2:
    boloinfo = sptpol_arcy_bolos.get_arcy_bolos()
    chaninfo = dict()
    for b in boloinfo:
        for jmodule in range(nmodules):
            chaninfo[b[0]+'_'+str(jmodule+1)] = b[1]

frequencies = np.sort(np.unique(np.asarray(chaninfo.values())))

for jbolo,id in enumerate(scans.keys()):
    # get bolometer number and module number from channel key
    jmodule = int(id[-1])
    jchannel = int(id.split('_')[1])
     
    flat_map = scans[id][~np.isnan(scans[id])]
    colormin = np.percentile(flat_map, 2)
    colormax = np.percentile(flat_map, 98)
    colormax = np.max([np.abs(colormin), np.abs(colormax)])
    colormin = -colormax
    
    # frequency sorting
    ind = np.arange(64)
    jfreq = ind[frequencies == chaninfo[id]]

    f = plt.figure(jmodule, figsize=(12,6), dpi=300)
    f.subplots_adjust(left=.01, bottom=.01, right=.99, top=.99, wspace=.0, hspace=.0)
    sub = plt.subplot(8, 8, jfreq+1)
    plt.imshow(scans[id], vmin=colormin, vmax=colormax, aspect='auto', interpolation='none')
    plt.text(0.65, 0.1, '%s\n%d Hz' % (id, chaninfo[id]), fontsize=4, transform=plt.gca().transAxes)
    sub.set_xticks([])
    sub.set_yticks([])
    
    plt.xlim([round(scans[id].shape[1]*0.1), round(scans[id].shape[1]*0.95)])
    plt.ylim([1, scans[id].shape[0]])
    
for jmodule in np.arange(1,nmodules+1):
    fig = plt.figure(jmodule)
    plt.savefig('figures_%s/map_module=%d.pdf' % (filestr, jmodule), dpi=2000)
