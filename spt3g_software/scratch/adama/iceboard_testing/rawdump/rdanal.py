import numpy, scipy.fftpack

data = numpy.memmap('outfile', numpy.int16, 'r')

def makepsd(start, stop, fftlen=65536):
	'''Compute PSD from sample start to stop, averaging PSDs of length fftlen samples. Returns the frequencies of the bins (in Hz) and the PSD'''
	freq = scipy.fftpack.rfftfreq(fftlen, 1./20e6)
	psd = numpy.zeros(len(freq))
	i = int(start)
	while i + fftlen <= stop:
		fftstub = scipy.fftpack.rfft(data[i:i+fftlen])
		psd += (numpy.conj(fftstub)*fftstub).real
		i += fftlen
	psd /= (i - start)/fftlen
	return freq, psd

def makespectrogram(start, fftlen=65536, slicelen=0.5, seconds=10):
	freq = scipy.fftpack.rfftfreq(fftlen, 1./20e6)
	psds = []

	i = 0
	while i < seconds:
		psds.append(makepsd(start + i*20e6, start + i*20e6 + slicelen*20e6, fftlen=fftlen)[1])
		i += slicelen
	
	return freq, numpy.asarray(psds)

