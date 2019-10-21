from spt3g import core, coordinateutils, mapmaker, mapspectra
from spt3g.mapspectra import map_analysis
import numpy, pickle, sys, copy

try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError

qu = True 
args = []
for i in sys.argv[1:]:
	args.append(i)

def apod_mask(fr, r=100):
	edge_kernel = numpy.asarray(numpy.matrix(numpy.blackman(r)).transpose()*numpy.matrix(numpy.blackman(r)))
	tt = copy.copy(fr['Wpol'].TT)
	return mapspectra.apodmask.generate_apodization_mask(tt, edge_kernel=edge_kernel)

#c = numpy.linspace(10**.5, 5000**.5, 200)**2
c = numpy.arange(25, 5000, 25)
ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(c)

f = list(core.G3File(args[0]))[0]

try:
	apod = pickle.load(open('apod.pkl', 'rb'))
except FileNotFoundError:
	apod = apod_mask(f)
	pickle.dump(apod, open('apod.pkl', 'wb'))
numpy.asarray(apod)[numpy.asarray(apod) < 0] = 0

def addtocoadd(coadd, newmap):
	if coadd is None:
		return [copy.copy(newmap[0]), copy.copy(newmap[1]), copy.copy(newmap[2]), copy.copy(newmap[3])]
	coadd[0] += newmap[0]
	coadd[1] += newmap[1]
	coadd[2] += newmap[2]
	coadd[3] += newmap[3]
	return coadd
def frame_from_maps(t,q,u,w):
	f = core.G3Frame()
	f['T'] = t
	f['Q'] = q
	f['U'] = u
	f['Wpol'] = w
	return f
def explode_ps_dict(d):
	if qu:
		return d['TT'], d['QQ'], d['UU']
	else:
		return d['TT'], d['EE'], d['BB']

coadd = None
ms = []
for i in args:
	for f in core.G3File(i):
		print(i)
		ms.append((f['T'], f['Q'], f['U'], f['Wpol']))
		coadd = addtocoadd(coadd, ms[-1])
tauto,qauto,uauto = explode_ps_dict(map_analysis.calculateCls(frame_from_maps(*coadd), apod_mask=numpy.asarray(apod), qu=qu, ell_bins=ell_bins))
coadd = mapmaker.mapmakerutils.remove_weight(*coadd)

noisecls = []
xcls = []
pairs = [(i, j) for i in range(len(ms)) for j in range(len(ms)) if i != j]
#pairs = [(i, j) for i in range(len(ms)) for j in range(len(ms)) if i != j and i % 10 == 0 and j % 10 == 0]
for i,j in pairs:
	print(i,j)
	f1 = frame_from_maps(*ms[j])
	f2 = frame_from_maps(*ms[i])
	a = mapmaker.mapmakerutils.remove_weight(*ms[j])
	b = mapmaker.mapmakerutils.remove_weight(*ms[i])
	dt, dq, du = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
	dqf, duf = mapspectra.basicmaputils.flatten_pol(dq,du)
	tps, qps, ups = mapspectra.basicmaputils.get_map_cls(dt,dqf,duf, numpy.asarray(apod), ell_bins, qu=qu)
	crosst, crossq, crossu = explode_ps_dict(map_analysis.calculateCls(f1, f2, apod_mask=numpy.asarray(apod), qu=qu, ell_bins=ell_bins))
	noisecls.append((tps, qps, ups))
	xcls.append((crosst, crossq, crossu))
	
noisecls = numpy.asarray(noisecls)/2/len(ms) # Each half has ~ twice the variance
noise = numpy.mean(noisecls, axis=0)
xcls = numpy.asarray(xcls)

meanpower = numpy.mean(xcls, axis=0)
t,q,u = meanpower
noisestdpower = numpy.std(xcls, axis=0)/numpy.sqrt(len(xcls))

pylab.figure()
#pylab.plot(c, noise[0]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), label='T')
pylab.plot(c, noise[1]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), label='Q')
pylab.plot(c, noise[2]**0.5 / (core.G3Units.arcmin * core.G3Units.uK), label='U')
#pylab.loglog()
pylab.semilogx()
pylab.legend()
pylab.xlabel('Ell')
pylab.ylabel('Depth (uK-arcmin)')

pylab.figure()
# XXX units incorrect
#pylab.errorbar(c, t *c*(c+1), noisestdpower[0]*c*(c+1), label='T')
pylab.errorbar(c, q / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), noisestdpower[1] / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), label='Q')
pylab.errorbar(c, u / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), noisestdpower[2] / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1),label='U')
pylab.semilogx()
pylab.legend()
pylab.xlabel('Ell')
pylab.ylabel('C_ell')
pylab.grid()
