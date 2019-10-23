from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy

def apod_mask(fr, r=100):
	edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
	tt = copy.copy(fr['Wpol'].TT)
	return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)

def ps_from_map(fr, apod_mask, qu=False):
        '''
	code to make a power spectrum from a map
        '''
        t,q,u = mapmaker.mapmakerutils.remove_weight(fr['T'], fr['Q'], fr['U'], fr['Wpol'])
        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)

        cen_ells = numpy.linspace(10, 2500, 200)
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
        t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu)

        return cen_ells, t_cls, e_cls, b_cls


def ps_from_map_double(fr1,fr2, apod_mask, qu=False):
        '''
	code to make a power spectrum from the difference of two maps.
        '''

        t1,q1,u1 = mapmaker.mapmakerutils.remove_weight(fr1['T'], fr1['Q'], fr1['U'], fr1['Wpol'])
        t2,q2,u2 = mapmaker.mapmakerutils.remove_weight(fr2['T'], fr2['Q'], fr2['U'], fr2['Wpol'])
        t = (t1-t2)/2.
        q = (q1-q2)/2.
        u = (u1-u2)/2.
        qf, uf = mapspectra.basicmaputils.flatten_pol(q,u)
        cen_ells = numpy.linspace(10, 2500, 200)
        ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(cen_ells)
        t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu)

        return cen_ells, t_cls, e_cls, b_cls

#Define the two groups we will use for the subtraction

#index 0 = JT
#index 1 = J/N/cleaned
#index2 = J/N/uncleaned

bun_1 = ['/big_scratch/javva/dontmessup2/bundles/jt_bundles/g3counterparts/bundles_082_sptpol.g3', '/big_scratch/javva/bundles_glitchless/bundles_3gpipe_082.g3', '/big_scratch/javva/bundles3g/uncleaned/bundles_3gpipe_082.g3']
bun_2 = ['/big_scratch/javva/dontmessup2/bundles/jt_bundles/g3counterparts/bundles_089_sptpol.g3', '/big_scratch/javva/bundles_glitchless/bundles_3gpipe_089.g3', '/big_scratch/javva/bundles3g/uncleaned/bundles_3gpipe_089.g3']

#gather power spectra of the differences of these groups

tcls = {}
qcls = {}
ucls = {}
ccls = {}
for i in [0,1,2]:
        print('do apod')
        f1 = list(core.G3File(bun_1[i]))[0]
        f2 = list(core.G3File(bun_2[i]))[0]
        apod = apod_mask(f1)
        print('calc cls')
        c,t,q,u = ps_from_map_double(f1,f2, numpy.asarray(apod), qu=True)
        ccls[i] = c
        tcls[i] = t
        qcls[i] = q
        ucls[i] = u

#grab JT's transfer function so I can divide it out for a low ell comparison
import pickle
n = pickle.load(open('new_bb_tf.pkl', "rb"), encoding = 'latin1')
def renorm(cen_ells,cls):
        normed = []
        clz = []
        for en1, i in enumerate(cen_ells):
                for en, m in enumerate(n['ell']):
                        if (i >= m[0] and i < n['ell'][en+1][0]):
                                normed = np.append(normed,cls[en1]/n['BB'][en])
                                clz = np.append(clz, i)
        return normed, clz

def plot_cls(cls, title):
        pylab.figure()

        i = 0
        print('renormalizing', i)
        a,b = renorm(ccls[i],cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))
        pylab.plot(b[3:], a[3:], label='TF Normalized JT')

        i = 1
        pylab.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)), label='Jessica/Nathan cleaned glitchless')
        i = 2
        pylab.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)), label='Jessica/Nathan uncleaned')

        pylab.semilogx()
        pylab.legend()
        pylab.xlim(50,3000)
        pylab.xlabel('Ell')
        pylab.ylabel('Depth (uK-arcmin)')
        pylab.title(title+' subtract bundle 82- 89')
        pylab.savefig('plots/renorm_bunsubtract_'+title+'89_82.png')
plot_cls(tcls, 'T')
plot_cls(qcls, 'Q')
plot_cls(ucls, 'U')


def plot_cls_ratio(cls, title):
        pylab.figure()

        #2000 is at index 159                                                                            
        i = 0
        print('renormalizing', i)
        a,b = renorm(ccls[i],cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))
        ta,tb = renorm(ccls[i],tcls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))
        pylab.plot(b[3:], a[3:]/ta[3:], label='TF Normalized JT')
        i = 1
        pylab.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/(tcls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)), label='Jessica/Nathan cleaned glitchless')

        i = 2
        pylab.plot(ccls[i], (cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/(tcls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)), label='Jessica/Nathan uncleaned')


        #pylab.loglog()                                                                                  
        pylab.semilogx()
        pylab.legend()
        pylab.xlim(50,3000)
        pylab.xlabel('Ell')
        pylab.ylabel('ratio %s/tcls'%title)
        pylab.title(title+' bundle 82-89 ratio')
        pylab.savefig('plots/ratio_bunsubtract_'+title+'.png')
plot_cls_ratio(qcls, 'Q')
plot_cls_ratio(ucls, 'U')
