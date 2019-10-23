import pylab, sys
from spt3g import core

rcw38cal = {}

p = core.G3Pipeline()
p.Add(core.G3Reader, filename=sys.argv[1:])
def accumulatecal(fr):
	global rcw38cal
	if 'RCW38FluxCalibration' not in fr:
		return
	for k,c in fr['RCW38FluxCalibration'].items():
		if k not in rcw38cal:
			rcw38cal[k] = []
		rcw38cal[k].append(c)
p.Add(accumulatecal)
p.Run()

scans = max([len(c) for c in rcw38cal.values()])
k = numpy.asarray([c for c in rcw38cal.values() if len(c) == scans])
k2 = k/numpy.nanmedian(k, axis=0)

pylab.hist([numpy.nanstd(numpy.abs(c))/numpy.nanmedian(numpy.abs(c)) for c in k2], bins=numpy.linspace(0,.1,100))
pylab.xlabel('RMS Fractional Fluctuation in Relative RCW38 Flux Calibration')
pylab.ylabel('Detectors per Bin')
