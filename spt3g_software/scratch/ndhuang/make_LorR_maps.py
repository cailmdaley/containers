'''
Script to make individual bolometer maps with boresight pointing for use in
pointing calibration (e.g. for RCW38).
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, mapmaker
import os

# Usage: makebolomaps.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
		      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
	       help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
	       help='Output filename')
P.add_argument('-s', '--source', action='store', 
	       default=None, help='name of source (defaults to SourceName)')
P.add_argument('-r', '--res', action='store', 
	       default=0.5, help='resolution [arcmin]')
P.add_argument('-L', '--left', action = 'store_true', 
               help = 'Make left only maps (default is right only), where left is not necessarily left, but right is never left')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin

# Guess the list of bolos to use and other metadata for configuration
bolos = None
if not bolos:
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if 'RawTimestreams_I' in frame:
                bolos = frame['RawTimestreams_I'].keys()
                starttime = frame['RawTimestreams_I'].start
                if args.source is None:
                    # Grab source from data, removing any '-pixelraster' from
                    # the end that indicates a fast point (every pixel rastered
                    # over the source)
                    args.source = frame['SourceName'].replace('-pixelraster', '')
                break
        if bolos is not None:
            break

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    args.source, x_len = 3.*core.G3Units.deg/res,
    y_len = 3.*core.G3Units.deg/res, res = res,
    proj = mapmaker.mapmakerutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)
pipe.Add(core.DeduplicateMetadata)

# Combine our various in-progress cal data
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# cut left or right
class lrcutter(object):
    def __init__(self, left = False):
        self.scan_count = -1
        self.left = int(left)
        
    def __call__(self, fr):
        if fr.type != core.G3FrameType.Scan:
            return
        self.scan_count += 1
        return bool((self.scan_count + self.left) % 2)
pipe.Add(lrcutter, left = args.left)

# Convert to Watts, applying various corrections
pipe.Add(std_processing.CalibrateRawTimestreams, output='BoloMapTimestreams',
  units=core.G3TimestreamUnits.Power)

# Reset all bolometer pointing that we think we already know to 0 to get boresight-pointing maps
pipe.Add(std_processing.MakeBoresightBolometerProperties)

# pipe.Add(core.Dump, added_message = 'Starting mapmaking')

# Kick off maps
map_id = 'bsmap'
ptng_key = 'MapPointing'
pipe.Add(mapmaker.MapInjector, map_id = map_id, maps_lst = [smstub],
         is_stub = True, make_polarized = False, do_weight = True)
pipe.Add(mapmaker.CalculatePointing,
         map_id = map_id, is_healpix = False, pointing_store_key = ptng_key,
         ts_map_key = 'RawTimestreams_I', 
         use_boresight_pointing = True,
         boresight_ra_key = 'OnlineBoresightRa',
         boresight_dec_key = 'OnlineBoresightDec',
         boresight_az_key = 'OnlineBorseightAz',
         boresight_el_key = 'OnlineBoresightEl')
pipe.Add(mapmaker.TodFiltering, 
         ts_in_key = 'RawTimestreams_I', 
         ts_out_key = 'FilteredBoloMapTimestreams',
         delete_input_ts = True,
         poly_order = 4, use_dynamic_source_filter = True)
pipe.Add(mapmaker.BinMap, map_id = map_id, 
         ts_map_key = 'FilteredBoloMapTimestreams',
         use_unity_weights = True,
         use_boresight_pointing = True,
         individual_bolos_to_map = bolos,
         pointing_store_key = ptng_key)

pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

