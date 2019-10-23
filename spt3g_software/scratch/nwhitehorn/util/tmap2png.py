from spt3g import core, coordinateutils
import matplotlib
import sys
matplotlib.use('Agg', warn=False)
import pylab

for fr in core.G3File(sys.argv[1]):
	if 'Id' not in fr or fr['Id'] != sys.argv[2]:
		continue

	pylab.imshow(fr['T']/fr['Wpol'].TT, vmin=-2e-3, vmax=2e-3, interpolation='None')
	pylab.savefig(sys.argv[3])
	pylab.clf()
	
