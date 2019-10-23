'''
Script to make CMB field maps.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, todfilter, coordinateutils
from spt3g.pointing import offline_pointing as op
import scipy.stats
import pdb
import os, numpy, copy

# Usage: makemaps.py <input files.g3> -o outputmaps.g3

P = ap.ArgumentParser(description='Maps for a CMB field (SPTpol 500d by default)',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[], help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
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
            if source is None:
                source = frame['SourceName']
                break
    if starttime is not None:
        break

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(std_processing.CalibrateRawTimestreams)

# Calculate detector pointing early and add point source mask

# Make an empty flatsky map for the source filtering
fs_stub = std_processing.CreateSourceMapStub(
    args.source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)
ps_map = mapmaker.pointsourceutils.make_point_source_map(fs_stub, '/spt/user/javva/lowell/ptsrc_config_ra0hdec-57p5_both_50mJy.txt') # XXX example file, use one appropriate for your field!
pipe.Add(mapmaker.MapInjector, map_id='ps_flatskymap',maps_lst=[fs_stub,], is_stub=False)

pipe.Add(mapmaker.CalculatePointing, map_id = 'ps_flatskymap', is_healpix=False,
         pointing_store_key = 'PixelPointing', boresight_ra_key = 'OnlineBoresightRa',
         boresight_dec_key = 'OnlineBoresightDec', ts_map_key = 'CalTimestreams')

# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering, ts_in_key='CalTimestreams',
    ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=False,
    poly_order=4, point_source_mask_id = 'ps_flatskymap',
    point_source_pointing_store_key = 'PixelPointing',
    filters_are_ell_based = True, lpf_filter_frequency=6600,
    mhpf_cutoff=20, boresight_az_key='OnlineBoresightAz',
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

# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'TimestreamsWatts', 'CalTimestreams'])

# Calculate Weights
pipe.Add(std_processing.weighting.AddPSDWeights)

# Split by band (could do other things: wafer, bolos, etc.)
pipe.Add(calibration.SplitByBand, input='PolyFilteredTimestreams',
    output_root='PolyFilteredTimestreams')

pipe.Add(core.Dump)

# Kick off maps
for band in ['90', '150', '220']: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id=args.source + '-%sGHz' % band,
      maps_lst=[smstub], is_stub=True, make_polarized=True, do_weight=True)
    pipe.Add(mapmaker.BinMap, map_id=args.source + '-%sGHz' % band,
      ts_map_key='PolyFilteredTimestreams%sGHz' % band,
      pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()


