'''
Script to make coadded maps of sources for calibration purposes, such as for
very-fast-point observations. Relies on existing bolometer properties map.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging
import scipy.stats
import os, numpy
import glob
import pickle


# Usage: makecoadd.py <input files.g3> -o outputmaps.g3 -s rcw38
P = ap.ArgumentParser(description='Single bolometer maps with boresight pointing',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
#P.add_argument('input_files', action='store', nargs='+', default=[],
#           help='Input files')
P.add_argument('-on', '--obsnum', action='store', default='output.g3',
           help='Output filename')
'''

P.add_argument('-o', '--output', action='store', default='output.g3',
           help='Output filename')
P.add_argument('-s', '--source', action='store', 
           default=None, help='name of source')
P.add_argument('-k', '--source-relative', action='store_true',
           default=False, help='calibrate in source-relative units rather than K_cmb')
P.add_argument('-p', '--polarized', action='store_true',
           default=False, help='make polarized maps')
P.add_argument('-r', '--res', action='store', type=float,
           default=0.1, help='resolution [arcmin]')
P.add_argument('-x', '--xlen', action='store', type=float,
           default=3, help='map width [deg]')
P.add_argument('-y', '--ylen', action='store', type=float,
           default=3, help='map height [deg]')
'''
args = P.parse_args()

obsnum = args.obsnum
'''
res = float(args.res) * core.G3Units.arcmin
ifs = args.input_files
src = args.source
op = args.ooutput
xlen = args.xlen
ylen = args.ylen
ispol = args.polarized
srcrel = args.source-relative


#This bit for hard coded
obsnum = '64146996' #oroginal good one for rcw38                                                    '''    
#bd = '/big_scratch/javva/crosstalk/templates/'
bd = '/spt/user/javva/crosstalk/templates/'
#res = float(.1)*core.G3Units.arcmin
res = float(.5)*core.G3Units.arcmin 
ifs = ['/spt/user/production/calibration/calframe/RCW38-pixelraster/%s.g3'%obsnum]+ glob.glob('/spt/data/bolodata/downsampled/RCW38-pixelraster/%s/0*'%obsnum)
src = 'rcw38'
op = bd+'SingleWaferMapsAmps_%s_time_constant_deconvolved_pt5_arcmin.g3'%obsnum
xlen = 3
ylen = 3
ispol = False
srcrel = False

#End hard coded bit

#tcs_dict = pickle.load(open('/spt/user/panz/data/time_constant_results/2019/tau_stat_2019.pkl','rb'))

tcs_dict = pickle.load(open('/home/sguns/tau_stat_2019.pkl','rb'))
tcs_dict = {k: v['avg']*core.G3Units.s for k,v in tcs_dict.items()}

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


# Grab the observation time in case we need it for planets
starttime = None
for fname in ifs:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            starttime = frame['RawTimestreams_I'].start
            sample_rate = frame['RawTimestreams_I'].sample_rate
            all_bolos = frame['RawTimestreams_I'].keys()
            if src is None:
                # Grab source from data, removing any '-pixelraster' from
                # the end that indicates a fast point (every pixel rastered
                # over the source)
                src = frame['SourceName'].replace('-pixelraster', '')
                break
    if starttime is not None:
        break

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    src, x_len = float(xlen)*core.G3Units.deg/res,
    y_len = float(ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=ifs)

# Deal with partially-complete calibration frames made during autoprocessing
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
pipe.Add(core.DeduplicateMetadata)

# Cut turnarounds
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

# Add Online pointing to scans.
# Clean up pre-existing timestreams
pipe.Add(
    core.Delete,
    keys=['OnlineBoresightAz', 'OnlineBoresightEl',
          'OnlineBoresightRa', 'OnlineBoresightDec', 'OnlineRaDecRotation']
)
pipe.Add(
    std_processing.CalculateCoordTransRotations,
    raw_az_key='RawBoresightAz',
    raw_el_key='RawBoresightEl',
    output='OnlineBoresight',
    transform_store_key='OnlineRaDecRotation',
    model='OnlinePointingModel',
    flags=['az_tilts', 'el_tilts', 'flexure', 'collimation', 'refraction']
    #, 'thermolin'] # Thermoline broken as of 4/13/17
)

# do some very nice early flagging
pipe.Add(std_processing.flagsegments.FieldFlaggingPreKcmbConversion, ts_key='RawTimestreams_I')

# Next to source-relative units (XXX: hardcode list of sources here), without the opacity
# correction this is likely meant to find (XXX: what about other uses?)
# Stop CalibrateRawTimestreams at Watts to (optionally) avoid the opacity correction.
pipe.Add(std_processing.CalibrateRawTimestreams, units=core.G3TimestreamUnits.Current,
    output='CalTimestreams')
#pipe.Add(calibration.ApplyTCalibration, InKCMB=not srcrel,
#    OpacityCorrection=not srcrel, Input='TimestreamsWatts', Output='CalTimestreams')

pipe.Add(todfilter.convolve_filter.ConvolveFilter, filter=DeconvolutionFilter(tcs = tcs_dict, sample_rate=sample_rate, all_bolos=all_bolos),keys=['CalTimestreams'], key_suffix='Deconvolved')


# Basic timestream filtering with dynamic source filter to handle bright point sources with
# unreliable pointing
pipe.Add(mapmaker.TodFiltering, ts_in_key = 'CalTimestreamsDeconvolved',
    ts_out_key = 'PolyFilteredTimestreams', use_dynamic_source_filter = True,
    poly_order = 4)

# Standard bolometer weighting
pipe.Add(std_processing.weighting.AddSigmaClippedWeight)

# do some very nice flagging
pipe.Add(std_processing.flagsegments.FlagNonResponsive, flag_key = 'Flags')
pipe.Add(timestreamflagging.flaggingutils.SigmaclipFlagGroupG3MapValue,
         m_key = 'TodWeights', low = 3.0, high = 3.0, flag_reason = 'BadWeight')
# remove flagged detectors
pipe.Add(timestreamflagging.RemoveFlagged, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'DeflaggedTimestreams')

# Split by band
pipe.Add(calibration.SplitByBand, input='DeflaggedTimestreams',
    output_root='DeflaggedTimestreams')

for band in ['90', '150', '220']:
    pipe.Add(calibration.SplitByWafer, input = 'DeflaggedTimestreams%sGHz' % band) 
# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'FilteredTimestreams', 'TimestreamsWatts', 'CalTimestreams', 'PolyFilteredTimestreams'])

ws = ['W172','W181','W206','W204','W203','W180','W174','W176','W188', 'W177']

tots = ['90GHz'+wi for wi in ws] + ['150GHz'+wi for wi in ws] + ['220GHz'+wi for wi in ws]

#pipe.Add(core.Dump)

# Kick off maps
for band in tots: # XXX should be automatic
    pipe.Add(mapmaker.MapInjector, map_id=src + '-%s' % band,
      maps_lst=[smstub], is_stub=True, make_polarized=ispol, do_weight=True)
pipe.Add(mapmaker.CalculatePointing, map_id=src + '-90GHzW180',
         pointing_store_key = 'PixelPointing',
         ts_map_key = 'DeflaggedTimestreams', trans_key='OnlineRaDecRotation')
for band in tots: # XXX should be automatic
    pipe.Add(mapmaker.BinMap, map_id=src + '-%s' % band,
             ts_map_key='DeflaggedTimestreams%s' % band,
             pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights',
             trans_key='OnlineRaDecRotation')
pipe.Add(core.Dump)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.Dump)

pipe.Add(core.G3Writer, filename=op)
pipe.Run()
