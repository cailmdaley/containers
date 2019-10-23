from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import glob
import matplotlib.pyplot as plt

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
        i = clum[0]
        plt.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.
                                          uK)), label='No bad maps')
        i = clum[1]
        plt.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.
                                          uK)), label='With bad maps')


        plt.semilogx()
        z = clum[0]
        znom = z.split('/')[-1].split('.')[0].split('_')[-1]
        plt.legend()
        plt.xlim(50,3000)
        plt.xlabel('Ell')
        plt.ylabel('Depth (uK-arcmin)')
        plt.title(title+ i.split('/')[-1])
        plt.savefig('/home/javva/spt3g_software/scratch/javva/make_jacknives/plots/chron_weights/comp_total_coadd_remove_bad_maps_%s.png'%title)

        plt.clf()

print('Making apodization mask...')
f = list(core.G3File('/big_scratch/javva/jackknives/chron_weights_remove_noisy_maps/total_coadd.g3'))[0
]
apod = apod_mask(f)

for i in [1]:
        clump = ['/big_scratch/javva/jackknives/chron_weights_remove_noisy_maps/total_coadd.g3','/big_scratch/javva/jackknives/chron_weights/total_coadd.g3']
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

