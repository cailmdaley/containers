import scipy.signal
import numpy, sys
from spt3g import core, todfilter, std_processing, calibration

'''
This is code to do an svd on data such that the modes can be used 
for cleaning future data.

run it with:

python pca_filt_solve.py output.g3 inputts.npy [calibration.g3]

'''
print('Reading data')

use_poldiff = len(sys.argv) > 3

if use_poldiff:
    for fr in core.G3File(sys.argv[3]):
        if 'BolometerProperties' in fr:
            bpm = fr['BolometerProperties']
        break

tss = numpy.load(sys.argv[2]).item()
lens = {t[0]: len(t[1]) for t in tss.items() if numpy.isfinite(t[1]).all()}
max_len = numpy.max(list(lens.values()))
for i,l in lens.items():
	if l != max_len:
		del tss[i]

print('Solving SVD')
if not use_poldiff:
	# Solve on all TOD
	tod_solver_matrix = numpy.column_stack([t for t in tss.values() if numpy.isfinite(t).all()])
else:
	# Solve on explicitly pair-differenced TOD
	# XXX: Any risk of mean corruption here? Should we subtract off the common mode?
	tod_solver_matrix = []
	invnamemapping = {bpm[b].physical_name: b for b in tss.keys()}
	for b,t in tss.items():
		if not bpm[b].physical_name.endswith('.X'):
			continue
		partner = bpm[b].physical_name[:-2] + '.Y'
		if partner not in invnamemapping or invnamemapping[partner] not in tss:
			continue
		tod_solver_matrix.append(t - tss[invnamemapping[partner]])
	tod_solver_matrix = numpy.column_stack([t for t in tod_solver_matrix if numpy.isfinite(t).all()])
svd = numpy.linalg.svd(tod_solver_matrix, full_matrices=False)

# svd now contains time-domain modes

#randbasis = numpy.random.normal(size=svd[0].shape)

def varred(N):
	# This function solves for spatial modes
	orthobasis = svd[0][:,:N]
	#orthobasis = randbasis[:,:N]
	#numpy.savetxt('svd-templates', orthobasis)

	# First bulk-solve the good detectors
	coeff = numpy.linalg.lstsq(orthobasis, numpy.column_stack([t for t in tss.values() if numpy.isfinite(t).all()]))[0].transpose()
	coeffs = dict()
	for i,k in enumerate([k for k in tss.keys() if numpy.isfinite(tss[k]).all()]):
		coeffs[k] = coeff[i]

	# Then mask out NaNs and solve the bad ones
	for k,t in tss.items():
		if numpy.isfinite(t).all():
			continue
		if not numpy.isfinite(t).any():
			continue
		coeffs[k] = numpy.linalg.lstsq(orthobasis[numpy.isfinite(t),:], t[numpy.isfinite(t)])[0]
	return coeffs

c = varred(1)

print('Median %.2f, sigma: %.2f' % (numpy.median(numpy.transpose([i[0] for i in c.values()])), numpy.std(numpy.transpose([i[0] for i in c.values()]))))

def coeffhist(n, coeffs, bins=numpy.linspace(-2,2,80)):
	hist = numpy.histogram([i[n] for i in coeffs.values()], bins=bins)
	return ((hist[1][:-1] + hist[1][1:])/2, hist[0])

def psdstuff(N):
	orthobasis = svd[0][:,:N]
	psd = lambda d: (numpy.fft.rfft(d)*numpy.fft.rfft(d).conj()).real
	bandvar = lambda d: numpy.sum(psd(d)[band])
	tss0 = list(tss.values())[0]
	freq = numpy.fft.fftfreq(tss0.size, 100./191)[:len(numpy.fft.rfft(tss0))]
	#band = numpy.logical_and(freq < 0.1, freq > 0.04)
	band = numpy.logical_and(freq < 0.5, freq > 0.2)
	c = varred(N)
	vars = numpy.asarray([bandvar(tss[i] - numpy.dot(c[i], orthobasis.T))/bandvar(tss[i]) for i in tss.keys()])
	print('Average variance reduction for %d components: %e' % (N, numpy.mean(vars)))
	return numpy.mean(vars)

c = varred(15)

# Write to cal frame
out = core.G3MapVectorDouble()
for d in c.keys():
	out[d] = c[d]
f = core.G3Frame(core.G3FrameType.Calibration)
f['NonsenseModes'] = out
core.G3Writer(sys.argv[1])(f)

