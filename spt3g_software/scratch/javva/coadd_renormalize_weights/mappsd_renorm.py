from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import glob
import matplotlib.pyplot as plt
import numpy as np
def apod_mask(fr, r=100):
	edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
	tt = copy.copy(fr['Wpol'].TT)
	return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)

def ps_from_map(fr, apod_mask, qu=False):
#	apod_good, msg  = mapspectra.apodmask.validate_apodization_mask(apod_mask, fr['Wpol'])
#	if not apod_good:
#		core.log_fatal(msg)

	t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'], fr['Wpol'])
	qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)

	cen_ells = numpy.linspace(10, 2500, 200)
	ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
	t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu)

	return cen_ells, t_cls, e_cls, b_cls

def plot_cls(cls, title, clum):
        print("plotting...")
        #2000 is at index 159                                                                            
        for i in clum:
                plt.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)), alpha = 0.5) 
 
        plt.semilogx()
#        plt.legend()
        plt.xlim(50,3000)
        plt.xlabel('Ell')
        plt.ylabel('Depth (uK-arcmin)')
        plt.title(title)
        plt.savefig('/home/javva/spt3g_software/scratch/javva/coadd_renormalize_weights/plots/pwsp/many_buns%s.png'%title)

        plt.clf()
print('Making apodization mask...')
f = list(core.G3File('/big_scratch/javva/reweight_coadds/10bundle_reweighted.g3'))[0]
try: 
        apod = pickle.load(open('apod.pkl', 'rb'))
        print('Found apod')
except:
        apod = apod_mask(f)

for i in [1]:
        files = []
        for i in np.arange(10,100,10):
                files = np.append('/big_scratch/javva/reweight_coadds/%sbundle_reweighted.g3'%i, files)
        files = np.append('/big_scratch/javva/reweight_coadds/100_small_bundle_reweighted.g3',files)
        clump = files
       
        print(clump)
        tcls = {}
        qcls = {}
        ucls = {}
        ccls = {}
        for cl in clump:
                f = list(core.G3File(cl))[0]
                c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True)
                ccls[cl] = c
                tcls[cl] = t
                qcls[cl] = q
                ucls[cl] = u
        plot_cls(tcls,'T',clump)
        plot_cls(qcls,'Q',clump)
        plot_cls(ucls,'U',clump)

