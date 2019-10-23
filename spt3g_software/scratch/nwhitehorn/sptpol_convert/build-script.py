from __future__ import print_function
import time, sys, os, numpy
from spt3g import core

arcfiles = {}
for f in os.listdir(sys.argv[1]):
	arcfiles[core.G3Time(f.split('.')[0]).time] = os.path.join(sys.argv[1], f)
sortedtimes = numpy.asarray(sorted(arcfiles.keys()))

out = {}
files = {}

for fn in sys.argv[2:]:
	f = open(fn, 'r')
	for line in f.readlines():
		if line[0] == '#':
			continue
		times = line.split()

		if (len(times)) > 2:
			# Source observation
			start = times[1]
			stop = times[2]
		else:
			# Feature bit
			start = times[0]
			stop = times[1]

		# Cut down to times of interest: 500d survey
		#if 'Jun-2016' not in start and 'Oct-2015' not in start:
		if '2013' not in start and '2014' not in start and '2015' not in start and '2016' not in start:
			continue
		if 'Jan-2013' in start or 'Feb-2013' in start or 'Mar-2013' in start or 'Apr-2013' in start:
			continue

		obsid = int((core.G3Time(start).time - core.G3Time('20170101_000000').time)/core.G3Units.s)
		out[obsid] = (start, stop)

		arcindices = numpy.searchsorted(sortedtimes, [core.G3Time(start).time, core.G3Time(stop).time])
		files[obsid] = [arcfiles[sortedtimes[i]] for i in range(arcindices[0]-1, arcindices[1])]

for obsid, times in out.items():
	print('JOB\t%s\tconvert.sub' % obsid)
	#print('SCRIPT POST')
	print('VARS\t%s\tObsStart="%s"' % (obsid, times[0]))
	print('VARS\t%s\tObsStop="%s"' % (obsid, times[1]))
	print('VARS\t%s\tInputFiles="%s"' % (obsid, ' '.join(files[obsid])))
