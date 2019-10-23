
from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy, pylab
from spt3g.mapspectra import map_analysis
import matplotlib.pylab as plt
import glob
try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError
rt = 'adam_gain_match'
digest = False
coadddiff = False
args = []
for i in sys.argv[1:]:
	if i == '-digest':
		digest = True
		continue
	if i == '-coadddiff':
		coadddiff = True
		continue
	args.append(i)
'''
def apod_mask(fr, r=100):
	edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
	tt = copy.copy(fr['Wpol'].TT)
	return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)
'''
#all_mps = glob.glob('/big_scratch/javva/maps_2019_lowell/%s*'%rt)
all_mps = glob.glob('/spt/user/javva/maps_2019/*/%s*/*/*.g3'%rt)
tcls = {}
qcls = {}
ucls = {}

for enu, i in enumerate((all_mps)):
        print('plotting %s out of %s'%(enu, len(all_mps)))
        try:
                f = list(core.G3File(i))[-6:]
        except:
                continue
        tcls[i] = {}
        qcls[i] = {}
        ucls[i] = {}
        for m1,m2 in [(0,3),(1,4),(2,5)]:
                apod = mapspectra.apodmask.make_border_apodization((f[m1]['Wpol']),
                        apod_type='cos', radius_arcmin=90.)
                t1,q1,u1 = mapmaker.mapmakerutils.remove_weight(f[m1]['T'], f[m1]['Q'], f[m1]['U'],f[m1]['Wpol'])
                t2,q2,u2 = mapmaker.mapmakerutils.remove_weight(f[m2]['T'], f[m2]['Q'], f[m2]['U'],f[m2]['Wpol'])
                numpy.asarray(apod)[numpy.asarray(apod) < 0] = 0
                mfr = core.G3Frame(core.G3FrameType.Map)
                mfr['T'] = t1-t2
                mfr['Q'] = q1-q2
                mfr['U'] = u1-u2
                ma = map_analysis.calculateCls(mfr, apod_mask = numpy.asarray(apod), qu=True, ell_min=10, ell_max = 2500, delta_ell = 12)
                bands = ['90','150','220']
                tcls[i][bands[m1]] = ma['TT']/2.
                qcls[i][bands[m1]] = ma['QQ']/2.
                ucls[i][bands[m1]] = ma['UU']/2.
        c = ma['ell']
        clz = {}
        clz['T'] = tcls
        clz['Q'] = qcls
        clz['U'] = ucls
        clz['c'] = ma['ell']
        with open('ps_summary_singles_%s_strict_cuts.pkl'%rt, 'wb') as handle:
                pickle.dump(clz, handle, protocol=pickle.HIGHEST_PROTOCOL)

#make summary plot
a = clz 
for en, i in enumerate(a['Q'].keys()):
        plt.plot(en, a['Q'][i][8]**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko')
        plt.plot(en, a['Q'][i][50]**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')
