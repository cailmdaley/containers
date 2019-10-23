from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import numpy as np
import pylab
import scipy.ndimage as ndimage
try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError

digest = False
coadddiff = True
coadd_file = '/big_scratch/javva/filtering_bundles/chron_weights_remove_noisy_maps_stringent_cuts_poly4_correct/kmask_bun_1.g3'
map_dir = '/big_scratch/javva/filtering_bundles/chron_weights_remove_noisy_maps_stringent_cuts_poly4_correct/'


def ps_from_map(fr, apod_mask, qu=False, coadd=None,ell_weights_2d_t = None, ell_weights_2d_e = None, ell_weights_2d_b = None):
        apod_good, msg  = mapspectra.apodmask.validate_apodization_mask(apod_mask, fr['Wpol'])
        if not apod_good:
                core.log_error(msg)

        t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'], fr['Wpol'])
        if coadd is not None:
                t -= coadd[0]
                q -= coadd[1]
                u -= coadd[2]

        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)

        cen_ells = numpy.linspace(10, 2500, 200)
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
        t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu,ell_weights_2d_t = ell_weights_2d_t, ell_weights_2d_e = ell_weights_2d_e, ell_weights_2d_b = ell_weights_2d_b)

        return cen_ells, t_cls, e_cls, b_cls

apod= pickle.load(open("/spt/user/javva/lowell/good_apodization_mask.pkl","rb"))
numpy.asarray(apod)[numpy.asarray(apod) < 0] = 0

coadd = None
if coadddiff:
	m = None
	w = None
	for i in [1]:
                i = coadd_file
                for f in core.G3File(i):
                        if m is None:
                                m = [f['T'], f['Q'], f['U']]
                                w = f['Wpol']
                        else:
                                m[0] += f['T']
                                m[1] += f['Q']
                                m[2] += f['U']
                                w += f['Wpol']
	coadd = mapmaker.mapmakerutils.remove_weight(m[0], m[1], m[2], w) 


stf = pickle.load(open('simkspacetf_100_correctapod_median_pluscor.pkl','rb'))
q_ft_sn = stf['q']
u_ft_sn = stf['u']
t_ft_sn = stf['t']

tcls1 = {}
qcls1 = {}
ucls1 = {}
for i in range(1,19):
        f = list(core.G3File(map_dir+'bundles_3gpipe_'+str(i)+'.g3'))[0]
        print(i)
        c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True, coadd = coadd)
        tcls1[i] = t
        qcls1[i] = q
        ucls1[i] = u

tcls = {}
qcls = {}
ucls = {}
for i in range(1,19):
        f = list(core.G3File(map_dir+'bundles_3gpipe_'+str(i)+'.g3'))[0]

        ep = pickle.load(open("/spt/user/javva/lowell/2dpsd/bun"+str(i)+".pkl","rb"))
        k_weights_t = ep['t']                                        
        k_weights_q = ep['q']                                          
        k_weights_u = ep['u']          
        
        c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True,ell_weights_2d_t=(1./t_ft_sn)*(np.abs(k_weights_t**2)),ell_weights_2d_e = (1./q_ft_sn)*(np.abs(k_weights_q**2)),ell_weights_2d_b = (1./u_ft_sn)*(np.abs(k_weights_u**2)), coadd = coadd)

        tcls[i] = t
        qcls[i] = q
        ucls[i] = u

def plot_cls_tf_accounted(cls,title,color, label , tff):
        a = pickle.load(open(tff,'rb'))
        if title == 'Q':
                cll = 'q'
        if title == 'U':
                cll = 'u'
        if title == 'T':
                cll = 't'
        for i in args:
                pylab.plot(c, (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/a[cll],color)

        pylab.semilogx()
        pylab.xlabel('Ell')
        pylab.ylabel('Depth (uK-arcmin)')
        pylab.title('Noise Power '+title)

p4 = 'poly4TF.pkl'
p4cor = 'poly4_ind_bunds_tfcorr.pkl' #THIS IS MAYBE NOT Most RECENT RERUN CODE JA                        
plot_cls(qcls, 'Q', 'r', 'k-weighted', p4corr)
plot_cls(qcls1, 'Q', 'b', 'vanilla', p4)

