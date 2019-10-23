'''
Script to make CMB field maps.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
from spt3g.timestreamflagging.glitchfinding import get_num_glitches
from spt3g.pointing import offline_pointing as op
import scipy.stats
import os, numpy

# Usage: makemaps.py <input files.g3> -o outputmaps.g3
P = ap.ArgumentParser(description='Maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
           help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
P.add_argument('-s', '--source', action='store', 
           default=None, help='name of source')
P.add_argument('-r', '--res', action='store', 
           default=1.0, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', 
           default=45, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', 
           default=25, help='map height [deg]')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# Grab the observation time in case we need it for planets
starttime = None
for fname in args.input_files:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            starttime = frame['RawTimestreams_I'].start
            if args.source is None:
                args.source = frame['SourceName']
                break
    if starttime is not None:
        break

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Cut uncalibratable detectors early
pipe.Add(calibration.elnod_analysis.RotateIQ, i_rotate_key='RawTimestreamsRotated')
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, input='RawTimestreamsRotated', output='FilteredTimestreams')

# Invert crosstalk and convert to Watts.
pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='FilteredTimestreams', output='TimestreamsWatts', units=core.G3TimestreamUnits.Power, ignore_missing=True)

# Next to source-relative units
pipe.Add(calibration.ApplyTCalibration, Input='TimestreamsWatts', Output='CalTimestreams')

# Redo pointing
if True:
	pipe.Add(core.Delete, keys=['RawBoresightAz', 'RawBoresightEl', 'OnlineBoresightAz', 'OnlineBoresightEl', 'OnlineBoresightRa', 'OnlineBoresightDec'])
	pipe.Add(std_processing.pointing.NaiveBoresightPointing)
pipe.Add(op.CorrectBoresightPointing, model='OnlinePointingModel',
         flags=['az_tilts', 'el_tilts', 'flexure','collimation',
                'refraction']) #, 'thermolin'])
pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
         az_timestream='OnlineBoresightAz', el_timestream='OnlineBoresightEl',
         ra_timestream='OnlineBoresightRa', dec_timestream='OnlineBoresightDec')

# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering, ts_in_key='CalTimestreams',
    ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=True,
    # XXX: dynamic source filter is turbobad!
    poly_order=4,
    filters_are_ell_based=True, lpf_filter_frequency=6600,
    mhpf_cutoff=0.1, boresight_az_key='OnlineBoresightAz',
    boresight_el_key='OnlineBoresightEl')

def notch_filter(fr):
    # XXX Should do this properly, but it's annoying to insert this in the
    # middle of the above.
    if fr.type != core.G3FrameType.Scan:
        return
    base = fr['PolyFilteredTimestreams'].values()[0]
    line = 1.53*core.G3Units.Hz
    bw = 0.1*core.G3Units.Hz

    filt = numpy.ones(len(base))
    freqs = numpy.fft.fftfreq(len(base), 1.0/base.sample_rate)
    filt[numpy.abs(freqs - line) < bw/2] = 0
    filt[numpy.abs(freqs - 2*line) < bw/2] = 0
    filt[numpy.abs(freqs - 3*line) < bw/2] = 0

    out_map = core.G3TimestreamMap()
    todfilter.fft_filter_mem_friendly(fr['PolyFilteredTimestreams'], filt,
      out_map, False, None)
    del fr['PolyFilteredTimestreams'] # XXX: in-place I think is fine here? could invent new name
    fr['PolyFilteredTimestreams'] = out_map
pipe.Add(notch_filter)

# Split by band
pipe.Add(calibration.SplitTimestreamsByBand, input='PolyFilteredTimestreams',
    output_root='PolyFilteredTimestreams')

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams'])

# Quick and dirty weight calculation
class weights(object):
    def __init__(self, indata, outdata):
        self.indata = indata
        self.outdata = outdata
    def __call__(self, fr):
        if 'CalibratorResponseSN' in fr:
            self.calsn = fr['CalibratorResponseSN']
        if 'PolyFilteredTimestreams' not in fr:
            return
        w = core.G3MapDouble()
        for k,ts in fr[self.indata].iteritems():
            # Cut detectors with too many glitches, all zero, bad cal S/N
            # XXX: find and kill ones with bad lines?
            if not numpy.isfinite(ts).all():
                w[k] = 0
            elif (ts == 0).all():
                w[k] = 0
            elif k not in self.calsn or self.calsn[k] < 20:
                w[k] = 0
            else:
                glitches = get_num_glitches(ts, [5,9])
                if glitches[0] > 4 or glitches[1] > 0:
                    w[k] = 0
                    continue
                w[k] = \
                  1./numpy.var(scipy.stats.sigmaclip(ts, 2.5, 2.5).clipped)
                if not numpy.isfinite(w[k]):
                    w[k] = 0
        # Clamp weights to twice the median; crazy high weights are probably
        # bogus
        # XXX: kill detectors with high weights rather than clamping?
        medweight = numpy.median(w.values())
        for k,v in w.iteritems():
            if v > 3*medweight:
                w[k] = 0 #2*medweight
            if v < 0.5*medweight:
                w[k] = 0 #0.5*medweight
            
        fr[self.outdata] = w
        print(self.indata, numpy.sum(numpy.asarray(w.values()) != 0))
pipe.Add(weights, indata='PolyFilteredTimestreams150GHz', outdata='TodWeights150GHz')
pipe.Add(weights, indata='PolyFilteredTimestreams90GHz', outdata='TodWeights90GHz')

pipe.Add(core.Dump)

# Kick off maps
for band in ['90', '150']: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id=args.source + '-%sGHz' % band,
      maps_lst=[smstub], is_stub=True, make_polarized=True, do_weight=True)
pipe.Add(mapmaker.CalculatePointing, map_id=args.source + '-150GHz',
  is_healpix = False, pointing_store_key = 'PixelPointing',
  ts_map_key = 'PolyFilteredTimestreams',
  boresight_ra_key='OnlineBoresightRa', boresight_dec_key='OnlineBoresightDec')
for band in ['90', '150']: # XXX should be automatic
    pipe.Add(mapmaker.BinMap, map_id=args.source + '-%sGHz' % band,
      ts_map_key='PolyFilteredTimestreams%sGHz' % band,
      pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights%sGHz' % band)

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

