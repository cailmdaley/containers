import numpy as np
import matplotlib.pyplot as plt
import dfmuxtools
import cPickle as pickle
import sptpol_arcy_bolos
import os
import hwm3channels

bolo_info = sptpol_arcy_bolos.get_arcy_bolos()
hwm_num = 3
filestr = 'dfmux_idf_20151213_095735_to_20151213_120041.g3'
#filestr = 'dfmux_idf_20151213_121402_to_20151213_142309.g3'
os.mkdir('figures_%s' % filestr)
maps_left, maps_right = dfmuxtools.make_bolo_maps('/data/sptdaq/iceboard_testing/dfmux/processed/' + filestr, 20)
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

for jdirection,map in enumerate([maps_left, maps_right]):
    for jbolo,id in enumerate(map.keys()):
        # get bolometer number and module number from channel key
        jmodule = int(id[-1])
        jchannel = int(id.split('_')[1])
        
        flat_map = map[id][~np.isnan(map[id])]
        colormin = np.percentile(flat_map, 2)
        colormax = np.percentile(flat_map, 98)
        colormax = np.max([np.abs(colormin), np.abs(colormax)])
        colormin = -colormax

        # frequency sorting
        ind = np.arange(64)
        jfreq = ind[frequencies == chaninfo[id]]

        f = plt.figure(jdirection*10 + jmodule, figsize=(12,6), dpi=300)
        f.subplots_adjust(left=.01, bottom=.01, right=.99, top=.99, wspace=.0, hspace=.0)
        nplots = len(f.get_axes())
        sub = plt.subplot(8, 8, jfreq+1)
        plt.imshow(map[id], vmin=colormin, vmax=colormax, aspect='auto', interpolation='none')
        plt.text(0.65, 0.1, '%s\n%d Hz' % (id, chaninfo[id]), fontsize=4, transform=plt.gca().transAxes)
        sub.set_xticks([])
        sub.set_yticks([])

        plt.xlim([round(map[id].shape[1]*0.1), round(map[id].shape[1]*0.95)])
        plt.ylim([1, map[id].shape[0]])
    
    for jmodule in np.arange(1,nmodules+1):
        fig = plt.figure(jdirection*10 + jmodule)
        plt.savefig('figures_%s/map_module=%d_direction=%d.pdf' % (filestr, jmodule, jdirection), dpi=2000)
