import time, numpy
from spt3g import core

pipe = core.G3Pipeline()

scanlength = 100000

pipe.Add(core.G3InfiniteSource(core.G3FrameType.Scan))
def buildscan(fr):
	fr['ScanStart'] = core.G3Time(int(time.time()*core.G3Units.s))
	ts = core.G3TimestreamMap()
	ts.start = core.G3Time(int(time.time()*core.G3Units.s))
	for bolo in ['Fake1', 'Fake2', 'Fake3', 'Fake4', 'Fake5']:
		ts[bolo] = core.G3Timestream(numpy.random.normal(size=scanlength))
		ts[bolo].units = core.G3TimestreamUnits.Tcmb
	fr['BoresightAz'] = core.G3Timestream(numpy.linspace(0, 360*core.G3Units.deg, scanlength))
	fr['BoresightEl'] = core.G3Timestream(50*core.G3Units.deg*numpy.ones(scanlength))
	ts.stop = core.G3Time(int(time.time()*core.G3Units.s))
	fr['CalTimestreams'] = ts
pipe.Add(buildscan)
def printfr(fr):
	#print numpy.asarray(fr['BoresightAz'])
	print(fr)
pipe.Add(printfr)
pipe.Add(core.G3Writer('/tmp/test2.g3'))

pipe.Run()
