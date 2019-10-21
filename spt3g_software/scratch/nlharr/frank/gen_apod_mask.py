from spt3g import core, todfilter, calibration, coordinateutils
from spt3g import mapmaker, todfilter, frbutils, util, mapspectra
from spt3g.mapspectra import basicmaputils, apodmask
import matplotlib.pyplot as plt
import numpy as np
import pickle


pntsrc = pickle.load(open('point_source_mask.pkl'))


ind = 0
#a = list(core.G3File('/home/javva/spt3g/spt3g_software/scratch/javva/pca_polyfilt/no_lpf_400_tryagain.g3'))
#a = list(core.G3File('/home/javva/spt3g/spt3g_software/scratch/javva/pca_polyfilt/lpf_400_tryagain.g3'))

a = list(core.G3File('/home/nlharr/tmp/maps500d/sptpol_season.g3'))
t,q,u = mapmaker.mapmakerutils.remove_weight(a[ind]['T'], a[ind]['Q'], a[ind]['U'], a[ind]['Wpol'])
w_pol = a[ind]['Wpol']

print("generating apodization mask")
convolution_size = 360
kernel = np.asarray(np.matrix(np.hanning(convolution_size)).transpose() *
                    np.matrix(np.hanning(convolution_size)))
apod_mask = apodmask.generate_apodization_mask(
    w_pol.TT, kernel, use_square_smoothing = False,
    bad_threshold_fraction = 0.03,
    #point_source_mask = pntsrc,

    post_mask_erosion = 60, pre_mask_erosion = 0)

v, msg = apodmask.validate_apodization_mask(apod_mask, w_pol)
if not v:
    core.log_warn(msg)
pickle.dump( apod_mask, open('widestest_inner_apod_mask.pkl', 'w'))

plt.imshow(apod_mask)
plt.show()
