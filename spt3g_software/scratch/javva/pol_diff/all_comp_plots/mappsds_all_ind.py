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
uK)), label='glitchless cleaned')
        i = clum[1]
        plt.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.
uK)), label='uncleaned')


        plt.semilogx()
        plt.legend()
        plt.xlim(50,3000)
        plt.xlabel('Ell')
        plt.ylabel('Depth (uK-arcmin)')
        plt.title(title+ i.split('/')[-1])
        plt.savefig('/home/javva/spt3g_software/scratch/javva/pol_diff/all_comp_plots/plots/'+i.split('/'\
)[-1].split('.')[0]+'_'+title+'_.png')
        plt.clf()

maps = glob.glob('/spt/user/nwhitehorn/sptpol/lowell/maps/*')
f_names = [i.split('/')[-1] for i in maps]

print('Making apodization mask...')
f = list(core.G3File(maps[0]))[0]
apod = apod_mask(f)

for z in f_names:
	clump = ['/spt/user/nwhitehorn/sptpol/lowell/maps/'+z, '/spt/user/nwhitehorn/sptpol/lowell/uncleanedmaps/'+z]
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

