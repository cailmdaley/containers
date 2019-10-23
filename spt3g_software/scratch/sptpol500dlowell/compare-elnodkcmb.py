from spt3g import core
import sys, numpy

out = {}

for i,f in enumerate(sys.argv[1:]):
	fr = core.G3File(f).next()
	for b,s in fr['ElnodKcmbSlopes'].iteritems():
		if b not in out:
			out[b] = numpy.zeros(len(sys.argv[1:]))
		out[b][i] = s

numpy.savetxt('compare', list(out.values()))
