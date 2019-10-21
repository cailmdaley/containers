from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy

sptpols = ['jt_bundles/g3counterparts/bundles_082_sptpol.g3','jt_bundles/g3counterparts/bundles_089_sptpol.g3']

g3s = ['/big_scratch/javva/bundles3g/bundles_3gpipe_082.g3','/big_scratch/javva/bundles3g/bundles_3gpipe_089.g3']

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


def apod_mask(fr, r=100):
	edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
	tt = copy.copy(fr['Wpol'].TT)
	return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)

def ps_from_map(t, q,u, apod_mask,qu=False):
#	apod_good, msg  = mapspectra.apodmask.validate_apodization_mask(apod_mask, fr['Wpol'])
#	if not apod_good:
#		core.log_fatal(msg)

#	t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'], fr['Wpol'])
	qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)

	cen_ells = numpy.linspace(10, 2500, 200)
	ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
	t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu)

	return cen_ells, t_cls, e_cls, b_cls

tcls = {}
qcls = {}
ucls = {}
ccls = {}


grp = sptpols

f1 = list(core.G3File(grp[0]))[0]
f2 = list(core.G3File(grp[1]))[0]

#apod_pol = apod_mask(f1,r=100)

t1,q1,u1 = mapmaker.mapmakerutils.remove_weight(f1['T'], f1['Q'], f1['U'], f1['Wpol'])
t2,q2,u2 = mapmaker.mapmakerutils.remove_weight(f2['T'], f2['Q'], f2['U'], f2['Wpol'])

tnu_pol = t1-t2
qnu_pol = q1-q2
unu_pol = u1-u2


grp = g3s

f1 = list(core.G3File(grp[0]))[0]
f2 = list(core.G3File(grp[1]))[0]

#apod_3 = apod_mask(f1,r=100)

t1,q1,u1 = mapmaker.mapmakerutils.remove_weight(f1['T'], f1['Q'], f1['U'], f1['Wpol'])
t2,q2,u2 = mapmaker.mapmakerutils.remove_weight(f2['T'], f2['Q'], f2['U'], f2['Wpol'])

tnu = t1-t2
qnu = q1-q2
unu = u1-u2

'''
i = 'g3s'

c,t,q,u = ps_from_map(tnu,qnu,unu, numpy.asarray(apod_3), qu=True)
ccls[i] = c
tcls[i] = t
qcls[i] = q
ucls[i] = u


def plot_cls(cls, title):
	pylab.figure()

	#2000 is at index 159
	i = 'pol'
	#pylab.plot(ccls[i], cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), l#abel=' unnormalized JT')
	nm = [cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)][0][159]
	print 'renormalizing', i
	a,b = renorm(ccls[i],cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))
	sf = a[159]
	ssf = nm/sf
	pylab.plot(b[3:], a[3:]/np.sqrt(2), label='TF Normalized JT')
	
	i = 'g3s'
	sf2 = [.96*cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)][0][159]
	ssf2 = nm/sf2
        pylab.plot(ccls[i], ((cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/.96)/np.sqrt(2), label='Jessica/Nathan')


	pylab.loglog()
	pylab.semilogx()
	pylab.legend()
	pylab.xlim(50,3000)
	pylab.xlabel('Ell')
	pylab.ylabel('Depth (uK-arcmin)')
	pylab.title(title+' bundle 82-89')
	pylab.savefig('plots/renorm'+title+'82_89.png')
plot_cls(tcls, 'T')
plot_cls(qcls, 'Q')
plot_cls(ucls, 'U')

'''
