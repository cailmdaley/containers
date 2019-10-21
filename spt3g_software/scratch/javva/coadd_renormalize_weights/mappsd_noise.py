from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy
import numpy as np
import pylab

try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError

digest = False
coadddiff = True
coadd_file = '/big_scratch/javva/reweight_coadds/100bundle_reweighted.g3'
files = []
for i in np.arange(10,100,10):
        files = np.append('/big_scratch/javva/reweight_coadds/%sbundle_reweighted.g3'%i, files)
files = np.append('/big_scratch/javva/reweight_coadds/100_small_bundle_reweighted.g3',files)
args = files

'''
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

def ps_from_map(fr, apod_mask, qu=False, coadd=None):
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
	t_cls,e_cls,b_cls = mapspectra.basicmaputils.get_map_cls(t,qf,uf, apod_mask, ell_bins, qu=qu)

	return cen_ells, t_cls, e_cls, b_cls

f = list(core.G3File('/big_scratch/javva/reweight_coadds/100_small_bundle_reweighted.g3'))[0]

try:
	apod = pickle.load(open('apod.pkl', 'rb'))
except FileNotFoundError:
	apod = apod_mask(f)
	pickle.dump(apod, open('apod.pkl', 'wb'))
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
	coadd = mapmaker.mapmakerutils.remove_weight(m[0], m[1], m[2], w) # XXX assumes #maps >> 1

#print("APOD MIN", numpy.min(apod))
tcls = {}
qcls = {}
ucls = {}
for i in args:
	f = list(core.G3File(i))[0]
	c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True, coadd=coadd)
	tcls[i] = t
	qcls[i] = q
	ucls[i] = u

if coadd is not None:
	f = core.G3Frame()
	f['T'] = m[0]
	f['Q'] = m[1]
	f['U'] = m[2]
	f['Wpol'] = w
	c,t,q,u = ps_from_map(f, numpy.asarray(apod), qu=True)
	tcls['Coadd'] = t
	qcls['Coadd'] = q
	ucls['Coadd'] = u

def plot_cls(cls, title):
        pylab.figure()
        for i in args:
                pylab.plot(c, cls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), label=i)
                #pylab.loglog()
        pylab.semilogx()
#        pylab.legend()
        pylab.xlabel('Ell')
        pylab.ylabel('Depth (uK-arcmin)')
        pylab.title('Noise Power '+title)
        pylab.savefig('plots/pwsp/noise/not_rw_%s.png'%title)

avgs = {}

def plot_cls_avg(cls, title):
        pylab.figure()
        tot_avg = []
        for idx, c_i  in enumerate(c):
                avg = []
                for i in args:
                        avg = np.append(cls[i][idx]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), avg) 
                tot_avg = np.append( tot_avg, np.mean(avg))
        avgs[title] = tot_avg
        pylab.plot(c, tot_avg, label = "Renormalized")
        old_np = pickle.load(open("avgps_nonnorm.pkl","rb"))
                #pylab.loglog()
        pylab.plot(c, old_np[title], label = 'Baseline')
        pylab.semilogx()
        pylab.legend()
        pylab.xlabel('Ell')
        pylab.ylabel('Depth (uK-arcmin)')
        pylab.title('Noise Power '+title)
        pylab.savefig('plots/pwsp/noise/avg_comp_%s.png'%title)


if digest:
	print('File\tT100\tT2000\tQ100\tQ2000\tU100\tU1000')
	bin100 = numpy.argmin(numpy.abs(c - 100))
	bin2000 = numpy.argmin(numpy.abs(c - 2000))
	for i in args:
		t = tcls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
		q = qcls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
		u = ucls[i]**0.5 / (core.G3Units.arcmin * core.G3Units.uK)
		print('%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f' % (i, t[bin100], t[bin2000], q[bin100], q[bin2000], u[bin100], u[bin2000]))
else:
	plot_cls_avg(tcls, 'T')
	plot_cls_avg(qcls, 'Q')
	plot_cls_avg(ucls, 'U')


import pickle
with open('avgps.pkl', 'wb') as handle:
    pickle.dump(avgs, handle, protocol=pickle.HIGHEST_PROTOCOL)

