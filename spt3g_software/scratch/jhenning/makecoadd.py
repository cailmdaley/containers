#!/usr/bin/env python

'''
Script to make coadded map.  Can turn on/off online/offline pointing, 
and specify which pointing model flags to use.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, todfilter, coordinateutils
import spt3g.pointing.offline_pointing as op
from spt3g.mapmaker.mapmakerutils import MakeMap
import scipy.stats
import os, numpy

# Usage: makebolomaps.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
		      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
	       help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
	       help='Output filename')
P.add_argument('--input-ts', action='store',
	       choices=['RawTimestreams_I', 'RawTimestreams_Q'],
	       default='RawTimestreams_I', help='Input timestream type')
P.add_argument('-s', '--source', action='store', 
	       default='rcw38', help='name of source')
P.add_argument('-c', '--calpath', action='store', default=None,
               help='path to directory of calibrator results')
P.add_argument('-x', '--xtalkpath', action='store', default=None,
               help='path to directory of crosstalk results')
P.add_argument('-r', '--res', action='store', 
	       default=0.25, help='resolution [arcmin]')
P.add_argument('--online', action='store', default=True, help='Online pointing')

P.add_argument('--flags', action='store', nargs='+',
               default=['az_tilts' 'el_tilts' 'flexure' 'collimation' 'refraction'], 
               help='Pointing model flags.')

args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# Guess the list of bolos to use
bolos = None
if not bolos:
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if 'CalibratorResponseSN' in frame:
                bolos = [k for k in frame['CalibratorResponseSN'].keys() if frame['CalibratorResponseSN'][k] > 20]
            if args.input_ts in frame:
                starttime = frame[args.input_ts].start
                break
print('Using %d bolometers passing the calibrator SN cut' % len(bolos))

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = 3.*core.G3Units.deg/res,
    y_len = 3.*core.G3Units.deg/res, res = res,
    proj = mapmaker.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)
pipe.Add(core.DeduplicateMetadata)

# Combine our various in-progress cal data
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

# Invert crosstalk and convert to Watts
if args.xtalkpath is not None:
    pipe.Add(xtalk.CrosstalkCleanedTimestreams, input=args.input_ts, output='BoloMapTimestreams', units=core.G3TimestreamUnits.Power)
else:
    pipe.Add(dfmux.ConvertTimestreamUnits, Input=args.input_ts, Output='BoloMapTimestreams', Units=core.G3TimestreamUnits.Power)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, input='BoloMapTimestreams', output='BoloMapTimestreamsFilt')
def weights(fr):
	if 'BoloMapTimestreamsFilt' not in fr:
		return
	w = core.G3MapDouble()
	for k,ts in fr['BoloMapTimestreamsFilt'].iteritems():
		if not numpy.isfinite(ts).all():
			w[k] = 0
		if (ts == 0).all():
			w[k] = 0
		w[k] = 1./numpy.var(scipy.stats.sigmaclip(ts, 2.5, 2.5).clipped)
		if k not in bolos or not numpy.isfinite(w[k]):
			w[k] = 0
	fr['TodWeights'] = w
	print(numpy.sum(numpy.asarray(w.values()) != 0))
pipe.Add(weights)

#Add an offline pointing field (online should already be included in new scanified bolodata).
if args.online == 'False':
    pipe.Add(op.CorrectBoresightPointing, online=args.online,
             flags=list(args.flags))

    pipe.Add(coordinateutils.azel.LocalToAstronomicalPointing,
             az_timestream='OfflineBoresightAz', el_timestream='OfflineBoresightEl',
             ra_timestream='OfflineBoresightRa', dec_timestream='OfflineBoresightDec')

pipe.Add(core.Dump)

#Set which Az/El fields to use
if args.online == 'True':
    print 'Using online pointing...'
    boresight_ra_key = 'OnlineBoresightRa'
    boresight_dec_key = 'OnlineBoresightDec'
elif args.online == 'False':
    print 'Using offline pointing...'
    boresight_ra_key = 'OfflineBoresightRa'
    boresight_dec_key = 'OfflineBoresightDec'

# Kick off maps
pipe.Add(MakeMap,
         # Timestream filtering -- these really shouldn't be part of MakeMap
         poly_order = 4,
         use_dynamic_source_filter = True, # Enable dynamic PS filtering

         # Actual map-making options
         map_in = smstub,
         map_id = 'bsmap',
         ts_in_key = 'BoloMapTimestreamsFilt',

         make_polarized = False,
         do_weight = True,
         timestream_weight_key = 'TodWeights',
         boresight_ra_key = boresight_ra_key,
         boresight_dec_key = boresight_dec_key)

pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run(profile=True)

