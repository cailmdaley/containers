'''
Script to make NETs from some random data. Relies on existing bolometer properties map.
'''
from spt3g import core, std_processing, dfmux, calibration, xtalk, todfilter, mapmaker
import scipy.stats
import os, numpy, sys

# Usage: nets.py <input files.g3> 

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=sys.argv[1:])

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Cut uncalibratable detectors early
pipe.Add(calibration.elnod_analysis.RotateIQ, i_rotate_key='RawTimestreamsRotated')
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, input='RawTimestreamsRotated', output='FilteredTimestreams')

# Invert crosstalk and convert to Watts.
pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='FilteredTimestreams', output='CalTimestreams', units=core.G3TimestreamUnits.Power, ignore_missing=True)

# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering, ts_in_key='CalTimestreams',
    ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=True,
    poly_order=4)

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams'])

variances = {}
bands = {}
frequencies = {}

# Quick and dirty weight calculation
class vars(object):
    def __init__(self):
        pass
    def __call__(self, fr):
        if 'WiringMap' in fr:
            self.wiringmap = fr['WiringMap']
        if 'CalibratorResponseSN' in fr:
            self.calsn = fr['CalibratorResponseSN']
        if 'BolometerProperties' in fr:
            for k in fr['BolometerProperties'].keys():
                bands[k] = fr['BolometerProperties'][k].band/core.G3Units.GHz
        if 'PolyFilteredTimestreams' not in fr:
            return

        for k,ts in fr['PolyFilteredTimestreams'].iteritems():
            if not numpy.isfinite(ts).all():
                continue
            elif (ts == 0).all():
                continue
            elif k not in self.calsn or self.calsn[k] < 20:
                continue
            else:
                if k not in variances:
                    variances[k] = []
                # Store noise level assuming it is flat to Nyquist
                variances[k].append(
                  numpy.var(scipy.stats.sigmaclip(ts, 2.5, 2.5).clipped)/
                            (ts.sample_rate/core.G3Units.Hz/2))
            frequencies[k] = dfmux.HousekeepingForBolo(fr['DfMuxHousekeeping'], self.wiringmap, k).carrier_frequency
pipe.Add(vars)

pipe.Add(core.Dump)

pipe.Run()

nets = {k: (numpy.sqrt(numpy.nanmean(variances[k]) / 2), bands[k], frequencies[k]) for k in variances.keys()}

# Make hists
data = numpy.asarray(list(nets.values()))
for band in [90, 150, 220]:
    pylab.hist(data[:,0][data[:,1] == band]*1e18, label=str(band), bins=numpy.linspace(0, 200, 100))
pylab.legend()

pylab.xlabel('NEP [aW-rt(s)]')
pylab.ylabel('# of bolos')

pylab.figure()
for band in [90, 150, 220]:
	pylab.scatter(data[:,2][data[:,1] == band]/core.G3Units.MHz, data[:,0][data[:,1] == band]*1e18, marker='.', s=1, label=str(band))

pylab.legend()
pylab.ylabel('NEP [aW-rt(s)]')
pylab.xlabel('Readout Frequency (MHz)')

