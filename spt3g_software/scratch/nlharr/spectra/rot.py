from spt3g import core, todfilter, calibration, coordinateutils
from spt3g import mapmaker, todfilter, frbutils, util, mapspectra
from spt3g.mapspectra import basicmaputils, apodmask
import matplotlib.pyplot as plt
import numpy as np
import pickle


ind = 1
#a = list(core.G3File('/home/javva/spt3g/spt3g_software/scratch/javva/pca_polyfilt/no_lpf_400_tryagain.g3'))
#a = list(core.G3File('/home/javva/spt3g/spt3g_software/scratch/javva/pca_polyfilt/lpf_400_tryagain.g3'))

a = list(core.G3File('/home/nlharr/tmp/hpf/actuallpf_withcodefix.g3'))
t,q,u = mapmaker.mapmakerutils.remove_weight(a[ind]['T'], a[ind]['Q'], a[ind]['U'], a[ind]['Wpol'])
w_pol = a[ind]['Wpol']

if 0:
    print("generating apodization mask")
    convolution_size = 60
    kernel = np.asarray(np.matrix(np.blackman(convolution_size)).transpose() *
                np.matrix(np.blackman(convolution_size)))
    apod_mask = apodmask.generate_apodization_mask(
        w_pol.TT, kernel, use_square_smoothing = True, 
        bad_threshold_fraction = 0.03,
        post_mask_erosion = 0, pre_mask_erosion = 0)
    pickle.dump( apod_mask, open('apod_mask.pkl', 'w'))
else:
    apod_mask = pickle.load(open('apod_mask.pkl'))

print("estimating ps")

cen_ells = np.arange( 10, 20000, 20)
ell_bins = basicmaputils.get_reg_spaced_ell_bins(cen_ells)


qf, uf = basicmaputils.flatten_pol(q,u)
t_cls,e_cls,b_cls = basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins)

plt.semilogy(cen_ells, cen_ells**2 * t_cls)



a = list(core.G3File('/home/nlharr/tmp/hpf/lpf_withcodefix.g3'))
tn,qn,un = mapmaker.mapmakerutils.remove_weight(a[ind]['T'], a[ind]['Q'], a[ind]['U'], a[ind]['Wpol'])
w_pol = a[ind]['Wpol']
qn, un = basicmaputils.flatten_pol(qn,un)
tn_cls,en_cls,bn_cls = basicmaputils.get_map_cls(tn,qn,un, apod_mask, ell_bins)

plt.semilogy(cen_ells, cen_ells**2 * tn_cls)

plt.show()

