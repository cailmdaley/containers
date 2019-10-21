#!/usr/bin/env python
'''
Script to make individual bolometer maps with boresight pointing for use in
pointing calibration (e.g. for RCW38).
'''
import argparse as ap
import pickle
from spt3g import core, std_processing, mapmaker, dfmux, calibration, xtalk, coordinateutils, todfilter
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
from spt3g.std_processing.pointing import CalculateCoordTransRotations, CalculateLocalOffsetPointing

import os, numpy

# Usage: makebolomaps.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
		      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[],
	       help='Input files')
P.add_argument('-o', '--output', action='store', default='output.g3',
	       help='Output filename')
P.add_argument('-s', '--source', action='store', 
	       default=None, help='name of source (defaults to SourceName)')
P.add_argument('-r', '--res', action='store', type=float,
	       default=0.5, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float,
               default=3.0, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
               default=3.0, help='map height [deg]')
E = P.add_mutually_exclusive_group()
E.add_argument('-t', '--kcmb', action='store_true',
               default=False, help='calibrate in K_cmb')
E.add_argument('-k', '--source-relative', action='store_true',
               default=False, help='calibrate in source-relative units rather than K_cmb')
args = P.parse_args()

res = float(args.res) * core.G3Units.arcmin


def DeconvolutionFilter(tcs, sample_rate, all_bolos):
    out = {}
    klen = int(max(tcs.values())*sample_rate*3.)*2 + 1
    print("Using deconvolution kernel of %d samples (3x longest tc of %.2f)"%(klen, max(tcs.values())/core.G3Units.ms))
    ktimes = (numpy.arange(klen) - klen // 2)/sample_rate
    assert(ktimes[klen // 2] == 0.)
    for k in all_bolos:
        if k in tcs.keys():
            tc = tcs[k]
        else:
            print("Did not find tc for %s, substituting 10 us"%k)
            tc = 0.01*core.G3Units.ms
        kern = numpy.where(ktimes >= 0, numpy.exp(-ktimes/tc), 0.)
        ikern = numpy.fft.ifft(1./numpy.fft.fft(kern), len(kern))
        ikern = ikern[2:] #ad hoc fix to center the maximum                                                                                                                           
        ikern = numpy.real(ikern) / sum(numpy.real(ikern))
        if not numpy.argmax(ikern) == len(ikern) // 2:
            print(numpy.argmax(ikern), len(ikern) // 2)
        assert(ikern[-1] < 0.01*numpy.max(ikern))
        assert(numpy.argmax(ikern) == len(ikern) // 2)
        assert(numpy.abs(sum(ikern) - 1.) < 0.0001)
        out[k] = ikern
    some_key = list(out.keys())[123]
    print("Kernel length check: %d"%len(out[some_key]))
    return out

def InsertTC(frame, timeconstants, bolos, missing_val):
    if frame.type != core.G3FrameType.Scan:
        return
    else:
        tc = core.G3MapDouble()
        for bolo in bolos:
            if bolo in timeconstants.keys():
                tc[bolo] = timeconstants[bolo]
            else:
                tc[bolo] = missing_val
        frame['TimeConstants'] = tc
        return



# Guess the list of bolos to use and other metadata for configuration
bolos = None
if not bolos:
    for fname in args.input_files:
        for frame in core.G3File(fname):
            if 'RawTimestreams_I' in frame:
                bolos = frame['RawTimestreams_I'].keys()
                starttime = frame['RawTimestreams_I'].start
                sample_rate = frame['RawTimestreams_I'].sample_rate
                all_bolos = bolos
                if args.source is None:
                    # Grab source from data, removing any '-pixelraster' from
                    # the end that indicates a fast point (every pixel rastered
                    # over the source)
                    args.source = frame['SourceName'].replace('-pixelraster', '')
                break
        if bolos is not None:
            break

# Generate map stub
smstub = coordinateutils.FlatSkyMap(
    x_len = int(args.xlen*core.G3Units.deg/res),
    y_len = int(args.ylen*core.G3Units.deg/res),
    alpha_center = 0, delta_center = 0, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    pol_type = core.MapPolType.T,
    coord_ref = core.MapCoordReference.Local)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)

tcs_file = '/home/sguns/tau_stat_2019.pkl'

tcs_dict = pickle.load(open(tcs_file, 'rb'))
tcs_dict = {k: v['avg']*core.G3Units.s for k,v in tcs_dict.items()}

pipe.Add(InsertTC, timeconstants=tcs_dict, bolos=all_bolos, missing_val = 6.*core.G3Units.ms)

pipe.Add(core.DeduplicateMetadata)

# Combine our various in-progress cal data
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Interpolate over NaNs
pipe.Add(std_processing.InterpOverNans)

# Add Online pointing to scans.
# Clean up pre-existing timestreams
pipe.Add(
    core.Delete,
    keys=['OnlineBoresightAz', 'OnlineBoresightEl',
          'OnlineBoresightRa', 'OnlineBoresightDec', 'OnlineRaDecRotation']
)
pipe.Add(
    CalculateCoordTransRotations,
    raw_az_key='RawBoresightAz',
    raw_el_key='RawBoresightEl',
    output='OnlineBoresight',
    transform_store_key='OnlineRaDecRotation',
    model='OnlinePointingModel',
    flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
    #, 'thermolin'] # Thermoline broken as of 4/13/17
)

pipe.Add(CalculateLocalOffsetPointing, source=args.source,
         x_offset_key='AzOffset', y_offset_key='ElOffset',
         trans_key='OnlineRaDecRotation', ts_ref_key='OnlineBoresightAz',
         max_throw=args.xlen * core.G3Units.deg * 1.66)

pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OffsetRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
         bs_ra_key = 'AzOffset', bs_dec_key = 'ElOffset',
         do_bad_transform = True)

# Convert to Watts, applying various corrections
pipe.Add(std_processing.CalibrateRawTimestreams, output='TimestreamsWatts',
         units=core.G3TimestreamUnits.Current)
'''
if args.kcmb or args.source_relative:
    pipe.Add(calibration.ApplyTCalibration, InKCMB=not args.source_relative,
        OpacityCorrection=not args.source_relative, Input='TimestreamsWatts', Output='CalTimestreams')
else:
'''
pipe.Add(core.Rename, keys={'TimestreamsWatts': 'CalTimestreams'})

# Drop calibration frames so that we don't pollute the output file
pipe.Add(lambda fr: fr.type != core.G3FrameType.Calibration)

# Reset all bolometer pointing that we think we already know to 0 to get boresight-pointing maps
pipe.Add(std_processing.MakeBoresightBolometerProperties)

pipe.Add(core.Dump)
pipe.Add(todfilter.convolve_filter.ConvolveFilter, filter=DeconvolutionFilter(tcs = tcs_dict, sample_rate=sample_rate, all_bolos=all_bolos),keys=['CalTimestreams'], key_suffix='Deconvolved')

pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreamsDeconvolved',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = True,
    poly_order = 4)

pipe.Add(core.Delete, keys = ['CalTimestreams','CalTimestreamsDeconvolved'], type = core.G3FrameType.Scan)


pipe.Add(mapmaker.MapInjector, map_id='bsmap',
         maps_lst=[smstub], is_stub=True, make_polarized=False, do_weight=True)

pipe.Add(mapmaker.mapmakerutils.CalculateBoresightPointing, map_id='bsmap',
         pointing_store_key='PixelPointing', trans_key='OffsetRotation',
         expand=True, bolo_map_key='PolyFilteredTimestreams')

pipe.Add(mapmaker.BinMap, map_id='bsmap',
         ts_map_key='PolyFilteredTimestreams',
         pointing_store_key='PixelPointing', 
         trans_key='OffsetRotation',
         use_boresight_pointing = True,
         use_unity_weights = True,
         individual_bolos_to_map = bolos)


pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan) # Drop TOD
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

