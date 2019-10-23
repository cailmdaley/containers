from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy


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

f = list(core.G3File(sys.argv[1]))[0]


#print("APOD MIN", np.min(apod))
tcls = {}
qcls = {}
ucls = {}
ccls = {}
for i in sys.argv[1:]:
	f = list(core.G3File(i))[0]
	apod = apod_mask(f)
	c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True)
	ccls[i] = c
	tcls[i] = t
	qcls[i] = q
	ucls[i] = u



def plot_cls(cls, title):
	pylab.figure()

	#2000 is at index 159
	i = sys.argv[4]
	#pylab.plot(ccls[i], cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), l#abel=' unnormalized JT')
	nm = [cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)][0][159]
	print 'renormalizing', i
	a,b = renorm(ccls[i],cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))
	sf = a[159]
	ssf = nm/sf
	pylab.plot(b[3:], a[3:], label='TF Normalized JT')
	
	i = sys.argv[5]
	sf2 = [.96*cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)][0][159]
	ssf2 = nm/sf2
        pylab.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/.96, label='Jessica/Nathan')


	#pylab.loglog()
	pylab.semilogx()
	pylab.legend()
	pylab.xlim(50,3000)
	pylab.xlabel('Ell')
	pylab.ylabel('Depth (uK-arcmin)')
	pylab.title(title+' bundle'+ i.split('_')[-1].split('.')[0])
	pylab.savefig('renorm'+title+i.split('_')[-1].split('.')[0]+'.png')
plot_cls(tcls, 'T')
plot_cls(qcls, 'Q')
plot_cls(ucls, 'U')

import pickle
import os
cwd = os.getcwd()
os.chdir('/home/javva/')
n = pickle.load(open('/home/javva/spt3g_software/scratch/javva/bundles/bb500d_TF_CC_preprocessed.pkl'))
os.chdir(cwd)

def renorm(cen_ells,cls):
	normed = []
	clz = []
	for en1, i in enumerate(cen_ells):
		for en, m in enumerate(n['150']['ell_bins']):
			if (i >= m[0] and i < n['150']['ell_bins'][en+1][0]):
				normed = np.append(normed,cls[en1]/n['150']['BB'][en])
				clz = np.append(clz, i)
	return normed, clz
