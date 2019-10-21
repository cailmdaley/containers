'''
Script to make CMB field maps.
'''
import argparse as ap
from spt3g import core, std_processing, mapmaker, calibration, todfilter, coordinateutils, timestreamflagging
from spt3g.std_processing import flagsegments
from spt3g.coordinateutils.coordsysmodules import FillCoordTransRotations
import os, numpy, copy
import scipy.stats

# Usage: makemaps.py <input files.g3> -o outputmaps.g3

#Jessica example:python makemaps_pieces.py /spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/-37748841.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-37748841/0000.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-37748841/0001.g3  -o /home/javva/spt3g/spt3g_software/scratch/javva/javva/jessicatrail.g3

#Example with sims python makemaps.py /spt/user/nwhitehorn/sptpol/autoproc/calibration/calframe/ra0hdec-57.5/-44660745.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0000.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0001.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0002.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0003.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0004.g3 /spt/user/nwhitehorn/sptpol/fullrate/ra0hdec-57.5/-44660745/0005.g3 -v -s /home/javva/spt3g_software/scratch/javva/simulations/sim_map_TEB_0.fits -o sim_maps.g3

#core.set_log_level(core.G3LogLevel.LOG_DEBUG)

P = ap.ArgumentParser(description='Maps with boresight pointing',
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
P.add_argument('-v', '--no_vector_cleaning', action='store_true',
           help='make maps without vector cleaning')
P.add_argument('-lt', '--lead_trail', action='store_true',
           help='Flag if using lead-trail data, to do correct ell based weighting')
P.add_argument('-S', '--accumulate_ts_for_analysis', action='store',
           help='path to low-frequency timestreams to save')
P.add_argument('-s', '--sims_path', action='store',
           help='path to sim .fits map')
P.add_argument('--left', action='store_true', help='use only left-going scans')
P.add_argument('--right', action='store_true', help='use only right-going scans')

args = P.parse_args()

if args.left and args.right: # Both left and right is the default again
    args.left = False
    args.right = False

res = float(args.res) * core.G3Units.arcmin

# Grab the observation time in case we need it for planets
starttime = None
source = None
for fname in args.input_files:
    for frame in core.G3File(fname):
        if 'RawTimestreams_I' in frame:
            starttime = frame['RawTimestreams_I'].start
            if source is None:
                source = frame['SourceName']
                break
    if starttime is not None:
        break

print(source)
print(starttime)

# Generate map stub
smstub = std_processing.CreateSourceMapStub(
    source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=args.input_files)
pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)
# Select which scan directions we want
pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
if args.left:
	pipe.Add(std_processing.SelectScanDirection, left=True)
if args.right:
	pipe.Add(std_processing.SelectScanDirection, left=False)

# Do calibration and cut broken/glitchy detectors
pipe.Add(std_processing.CalibrateRawTimestreams, q_output='CalTimestreamsQ')
pipe.Add(flagsegments.FieldFlaggingPreKcmbConversion, ts_key='RawTimestreams_I')
pipe.Add(flagsegments.FieldFlaggingPostKcmbConversion, ts_key='CalTimestreamsQ')
pipe.Add(core.Delete, keys=['GlitchesNumberOf', 'GlitchesThresholds'])
pipe.Add(flagsegments.FieldFlaggingPostKcmbConversion, ts_key='CalTimestreams')
pipe.Add(timestreamflagging.FlagHighQWithPolyFilter, variance_threshold=5)
pipe.Add(timestreamflagging.flaggingutils.RemoveFlaggedTimestreams,
         input_ts_key = 'CalTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'FlaggedTimestreams_I'
         )
pipe.Add(timestreamflagging.GenerateFlagStats, flag_key = 'Flags')
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 'CalTimestreams'])
pipe.Add(core.Rename, keys={'FlaggedTimestreams_I': 'CalTimestreams'})

# Calculate detector pointing early and add point source mask

# Make an empty flatsky map for the source filtering
fs_stub = std_processing.CreateSourceMapStub(
    source, x_len = float(args.xlen)*core.G3Units.deg/core.G3Units.arcmin,
    y_len = float(args.ylen)*core.G3Units.deg/core.G3Units.arcmin, res = core.G3Units.arcmin,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)
ps_map = mapmaker.pointsourceutils.make_point_source_map(fs_stub, 'ptsrc_config_ra0hdec-57p5_both_50mJy.txt')
pipe.Add(mapmaker.MapInjector, map_id='ps_flatskymap',maps_lst=[fs_stub,], is_stub=False)

#only needed for data before feb 2018
pipe.Add(FillCoordTransRotations,
         transform_store_key = 'OnlineRaDecRotation',
         bs_az_key = 'RawBoresightAz', bs_el_key = 'RawBoresightEl',
         bs_ra_key = 'OnlineBoresightRa', bs_dec_key = 'OnlineBoresightDec',
         do_bad_transform = True)

pipe.Add(mapmaker.CalculatePointing,map_id = 'ps_flatskymap', pointing_store_key = 'FinePixelPointing', ts_map_key = 'CalTimestreams', trans_key='OnlineRaDecRotation')

class VectorModeCleaner(object):
    def __init__(self, InputI='FilteredTS_I150GHz', InputQ='FilteredTS_Q150GHz', OutputI='VectorCleanedI', OutputQ='VectorCleanedQ', AngularCutoff=10*core.G3Units.deg, ModeMedianAmplitudeCutoff=2.5, AltSolverInputI=None):
        self.inputi = InputI
        if AltSolverInputI is not None:
            self.modesi = AltSolverInputI
        else:
            self.modesi = InputI
        self.inputq = InputQ
        self.outputi = OutputI
        self.outputq = OutputQ
        self.calsn = None
        self.modeamp = ModeMedianAmplitudeCutoff
        self.angular_cutoff = AngularCutoff
    def __call__(self, fr):
        if 'CalibratorResponseSN' in fr:
            self.calsn = fr['CalibratorResponseSN']
        if self.inputi not in fr or self.modesi not in fr:
            return
        ts = {}
        ms = {}
        for i in fr[self.inputi].keys():
            if i not in fr[self.inputq]:
                continue
            if i not in self.calsn or self.calsn[i] < 20:
                continue
            
            tsi = numpy.asarray(fr[self.inputi][i]) + 1j*numpy.asarray(fr[self.inputq][i])
            msi = numpy.asarray(fr[self.modesi][i]) + 1j*numpy.asarray(fr[self.inputq][i])
            if not numpy.isfinite(tsi).all() or not numpy.isfinite(msi).all():
                continue
            ts[i] = tsi
            ms[i] = msi
        if len(ts) < 300:
            return False
        tod_solver_matrix = numpy.column_stack([t for t in ms.values()])
        svd = numpy.linalg.svd(tod_solver_matrix.imag, full_matrices=False)
        orthobasis = svd[0][:,:tod_solver_matrix.shape[1]]
        orthobasis /= numpy.std(orthobasis, axis=0) # All modes have unit variance

        # Get phase information for all modes
        modescoeff = numpy.linalg.lstsq(orthobasis, tod_solver_matrix)[0].transpose()
        if self.modesi != self.inputi:
            tod_solver_matrix = numpy.column_stack([t for t in ts.values()])
            coeff = numpy.linalg.lstsq(orthobasis, tod_solver_matrix)[0].transpose()
        else:
            coeff = modescoeff

        # Get template of "bad" modes by nulling good ones
        coeff[(numpy.abs(numpy.arctan(modescoeff.imag/modescoeff.real)) < self.angular_cutoff/core.G3Units.rad) | (numpy.abs(modescoeff) < self.modeamp*numpy.median(numpy.abs(modescoeff)))] = 0
        print('Fraction of modes remaining: ', numpy.sum(coeff == 0)/numpy.sum(numpy.isfinite(coeff)))
        newdata = tod_solver_matrix - numpy.matrix(orthobasis)*numpy.matrix(coeff).transpose()
        print('Starting Q/I variance: %.2f' % (numpy.var(tod_solver_matrix.imag)/numpy.var(tod_solver_matrix.real)))
        print('Ending Q/I variance: %.2f' % (numpy.var(newdata.imag)/numpy.var(newdata.real)))
        print('Ending/starting I variance: %.2f' % (numpy.var(newdata.real)/numpy.var(tod_solver_matrix.real)))

        outi = core.G3TimestreamMap()
        outq = core.G3TimestreamMap()
        for i,k in enumerate([t for t in ts.keys() if numpy.isfinite(ts[t]).all()]):
            outtsi = core.G3Timestream(newdata[:,i].real.flatten())
            outtsq = core.G3Timestream(newdata[:,i].imag.flatten())
            outtsi.units = outtsq.units = fr[self.inputi].values()[0].units
            outtsi.start = outtsq.start = fr[self.inputi].start
            outtsi.stop = outtsq.stop = fr[self.inputi].stop
            outi[k] = outtsi
            outq[k] = outtsq

        fr[self.outputi] = outi
        if self.outputq is not None:
            fr[self.outputq] = outq

class RemoveNonsense(object):
    def __init__(self, min_order=0, max_order=15, data='NonsenseModes',
                 input='CalTimestreams', output='CleanedTS',
                 silent_failure=False, anti_common_mode=False):
        self.coeffs = None
        self.data = data
        self.min_order = min_order
        self.max_order = max_order
        self.input = input
        self.output = output
        self.silent_failure = silent_failure
        self.anti_common_mode = anti_common_mode
    def __call__(self, frame):
        if self.data in frame:
            self.coeffs = frame[self.data]
        if 'BolometerProperties' in frame:
            bpm = frame['BolometerProperties']
            pixlist = {}
            for bolo, bp in bpm.iteritems():
                if not bp.physical_name.endswith('.Y') and not bp.physical_name.endswith('.X'):
                    continue
                pixname = bp.physical_name[:-2]
                if pixname not in pixlist:
                    pixlist[pixname] = [None, None]
                if bp.physical_name.endswith('.X'):
                    pixlist[pixname][0] = bolo
                else:
                    pixlist[pixname][1] = bolo
            self.pixlist = {a: b for a,b in pixlist.items() if None not in b}
        if self.input not in frame:
            return
        if self.coeffs is None and self.silent_failure:
            frame[self.output] = frame[self.input]
            return

        #global modes, d, s, dstart, last_frame, outts, ints, pixel, pair

        pixs_in_play = set()
        dets_in_play = set()
        for pixel,pair in list(self.pixlist.items()):
            if pair[0] is None or pair[1] is None:
                continue
            if pair[0] not in frame[self.input]:
                continue
            if pair[1] not in frame[self.input]:
                continue
            pixs_in_play.add(pixel)
            dets_in_play.add(pair[0])
            dets_in_play.add(pair[1])

        keys_to_use = (set(frame[self.input].keys()) & set(self.coeffs.keys())) - dets_in_play
        print('Functioning detectors at start: ', len(frame[self.input].keys()))
        print('Fitting modes from %d detectors' % len(keys_to_use))
        if len(keys_to_use) == 0:
            return False

        if self.data + 'TOD' not in frame:
            c = 1 if self.anti_common_mode else 0
            data_array = numpy.asarray([numpy.asarray(frame[self.input][k]) for k in keys_to_use])
            templates = numpy.ones(shape=(len(data_array),
                len(list(self.coeffs.values())[0])+c))
            for i,b in enumerate(keys_to_use):
                if b in self.coeffs:
                    templates[i,c:] = self.coeffs[b]

            # Solve for time-domain modes
            last_frame = frame
            modes = numpy.linalg.lstsq(templates, data_array)[0]
            outmodes = core.G3TimestreamMap()
            for i,m in enumerate(modes):
                outmodes[str(i)] = core.G3Timestream(m)
            frame[self.data + 'TOD'] = outmodes
        else:
            modes = numpy.asarray(frame[self.data + 'TOD'])


        out = core.G3TimestreamMap()
        for pixel in pixs_in_play:
            pair = self.pixlist[pixel]
            ts = [frame[self.input][p] for p in pair]
            s = ts[0] + ts[1]
            d = ts[0] - ts[1]
            dstart = ts[0] - ts[1]
            x = numpy.linalg.lstsq(modes.transpose(), d)[0]
            d -= numpy.asarray(numpy.matrix(x)*numpy.matrix(modes)).reshape(len(d))
            assert(numpy.var(d) < numpy.var(dstart))
            assert((numpy.asarray(d) != 0).all())
            out[pair[0]] = (s + d)/2.
            out[pair[1]] = (s - d)/2.
        frame[self.output] = out
        outts = out

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

tsmappings = {'CalTimestreams': 'CleanedTS'}
# Simulate detectors from which we are going to make maps
if args.sims_path:
    m = mapmaker.mapmakerutils.load_spt3g_map(args.sims_path)
    t_map = m['T']
    q_map = m['Q']
    u_map = m['U']
    pipe.Add(mapmaker.MapInjector, map_id='sim_map',maps_lst=[t_map,q_map,u_map], is_stub=False)

    def MoveTsKeysToList(frame, tsm_key, lst_key):
        if frame.type != core.G3FrameType.Scan:
            return
        frame[lst_key] = core.G3VectorString(frame[tsm_key].keys())
    pipe.Add(MoveTsKeysToList, tsm_key = 'CalTimestreams', lst_key = 'SimTimestreamKeys')

    interp_sim = True
    if not interp_sim:
        pipe.Add(mapmaker.CalculatePointing,map_id = 'sim_map', 
                 pointing_store_key = 'FinePixelPointing_sim', 
                 ts_map_key = 'CalTimestreams',
                 trans_key='OnlineRaDecRotation')

    pipe.Add(mapmaker.mapmakerutils.FillSimTodSegment,
             out_ts_key = 'SimTimestreams',
             ts_to_get_sample_rate_key = 'OnlineBoresightRa',             
             interp_sim = interp_sim,
             map_is_healpix = True,
             trans_key='OnlineRaDecRotation',
             sim_map_id = 'sim_map', 
             valid_ids = 'SimTimestreamKeys',
             sim_pointing_key = 'FinePixelPointing_sim' if not interp_sim else '',
         )
    pipe.Add(core.Dump)
    tsmappings['SimTimestreams'] = 'CleanedSimTS'
    print('Using simulated timestreams')


pipe.Add(lambda fr: 'CalTimestreams' not in fr or len(fr['CalTimestreams']) > 0)
pipe.Add(todfilter.MaskedPolyHpf, in_ts_map_key='CalTimestreamsQ',
    out_ts_map_key='DetrendedTS_Q', poly_order=4) # XXX: should have PS mask applied
pipe.Add(todfilter.MaskedPolyHpf, in_ts_map_key='CalTimestreams',
    out_ts_map_key='DetrendedTS_IData', poly_order=4) # XXX: should have PS mask applied
pipe.Add(calibration.SplitTimestreamsByBand, input='DetrendedTS_Q', output_root='DetrendedTS_Q')

for i,o in tsmappings.items():
    print('Filtering %s -> %s' % (i, o))
    pipe.Add(todfilter.MaskedPolyHpf, in_ts_map_key=i,
        out_ts_map_key='DetrendedTS', poly_order=4) # XXX: should have PS mask applied

    if not args.no_vector_cleaning:
        pipe.Add(calibration.SplitTimestreamsByBand, input='DetrendedTS', output_root='DetrendedTS_I')
        pipe.Add(VectorModeCleaner, InputI='DetrendedTS_I150GHz', InputQ='DetrendedTS_Q150GHz', AltSolverInputI='DetrendedTS_IData')
        pipe.Add(core.Delete, keys=['DetrendedTS'])
        pipe.Add(core.Rename, keys={'VectorCleanedI': 'DetrendedTS'})

    pipe.Add(core.Delete, keys=['DetrendedTS_I150GHz', 'DetrendedTS_I90GHz'])

    # Basic timestream filtering
    pipe.Add(mapmaker.TodFiltering, ts_in_key='DetrendedTS',
        ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=False,
        poly_order=4, point_source_mask_id = 'ps_flatskymap',
        point_source_pointing_store_key = 'FinePixelPointing',
        filters_are_ell_based = True, lpf_filter_frequency=6600,
        boresight_az_key='OnlineBoresightAz', boresight_el_key='OnlineBoresightEl')
    pipe.Add(core.Delete, keys=['DetrendedTS', 'DetrendedTS_Q', 'VectorCleanedQ'])

    pipe.Add(notch_filter)

    pipe.Add(lambda fr: fr.type != core.G3FrameType.Scan or 'PolyFilteredTimestreams' in fr) # Cut frames where there are no timestreams left

    pipe.Add(RemoveNonsense, input='PolyFilteredTimestreams', output=o,
        silent_failure=True)
    pipe.Add(core.Delete, keys=[i, 'PolyFilteredTimestreams'])
pipe.Add(core.Delete, keys=[i, 'DetrendedTS_IData'])

# Calculate Weights
class CalculateExplicitPairDifferences(object):
    def __init__(self):
        self.pixlist = None
    def __call__(self, fr):
        if 'BolometerProperties' in fr:
            bpm = fr['BolometerProperties']

            # Collect list of detectors sharing a pixel
            pixlist = {}
            for bolo, bp in bpm.iteritems():
                if not bp.physical_name.endswith('.Y') and not bp.physical_name.endswith('.X'):
                    continue
                pixname = bp.physical_name[:-2]
                if pixname not in pixlist:
                    pixlist[pixname] = [None, None]
                if bp.physical_name.endswith('.X'):
                    pixlist[pixname][0] = bolo
                else:
                    pixlist[pixname][1] = bolo
            self.pixlist = {a: b for a,b in pixlist.items() if None not in b}
        if 'CleanedTS' in fr:
            out = core.G3TimestreamMap()
            ints = fr['CleanedTS']
            for pix in self.pixlist.values():
                if len(pix) != 2 or pix[0] not in ints or pix[1] not in ints:
                    continue
                diffts = ints[pix[0]] - ints[pix[1]]
                out[pix[0]] = diffts
                out[pix[1]] = diffts
            fr['PairDifferencedTimestreams'] = out
pipe.Add(CalculateExplicitPairDifferences)
pipe.Add(lambda fr: 'PairDifferencedTimestreams' not in fr or len(fr['PairDifferencedTimestreams']) > 0) # Give up on scans with no complete detector pairs

if not args.lead_trail:
    print('Doing weighting ranged for normal observing cadence')
    pipe.Add(std_processing.weighting.AddPSDWeights, low_f=.1*core.G3Units.Hz, high_f=1*core.G3Units.Hz, input='PairDifferencedTimestreams', output='PairDiffWeights')

if args.lead_trail:
    print('Doing weighting ranged for lead-trail data scan speeds')
    pipe.Add(std_processing.weighting.AddPSDWeights, low_f=.05*core.G3Units.Hz, high_f=.5*core.G3Units.Hz, input='PairDifferencedTimestreams', output='PairDiffWeights')

# Clean up detritus
pipe.Add(core.Delete, keys=['PolyFilteredTimestreams', 'DetrendedTS', 'PairDifferencedTimestreams'])

if args.sims_path is not None:
    pipe.Add(core.Delete, keys=['CleanedTS'])
    pipe.Add(core.Rename, keys={'CleanedSimTS': 'CleanedTS'})

if args.accumulate_ts_for_analysis is not None:
    ts = {}
    calsn = None
    baddets = set() # List of detectors flagged in at least one scan
    def collect_ts(fr):
        if 'CalibratorResponseSN' in fr:
            global calsn
            calsn = fr['CalibratorResponseSN']
        if 'PolyFilteredTimestreams' not in fr:
            return
        for i in fr['PolyFilteredTimestreams'].keys():
            prefactor = 1
            if i not in calsn or calsn[i] < 20:
                if i not in baddets:
                    baddets.add(i)
                prefactor = numpy.nan
            if i not in ts:
                ts[i] = []
            ts[i].append(prefactor*numpy.asarray(fr['PolyFilteredTimestreams'][i]))
    pipe.Add(collect_ts)

def sigmaclipweights(fr):
	if 'PairDiffWeights' not in fr:
		return
	pdweights = fr['PairDiffWeights']
	clipped, low, high = scipy.stats.sigmaclip(pdweights.values(), 3, 2)
	out = core.G3MapDouble()
	for bolo, w in pdweights.items():
		if w >= low and w <= high:
			out[bolo] = w
		else:
			out[bolo] = 0
	fr['TodWeights'] = out
pipe.Add(sigmaclipweights)

# Split by band
pipe.Add(calibration.SplitTimestreamsByBand, input='CleanedTS',
    output_root='CleanedTS')

# Add maps for actual results
pipe.Add(core.Delete, keys=['FinePixelPointing', 'DetectorAlphaPointing', 'DetectorDeltaPointing'])
pipe.Add(lambda fr: fr.type != core.G3FrameType.Map)
sm_stub = std_processing.CreateSourceMapStub(
    source, x_len = float(args.xlen)*core.G3Units.deg/res,
    y_len = float(args.ylen)*core.G3Units.deg/res, res = res,
    proj = coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea,
    at_time = starttime)

pipe.Add(core.Dump)

def printweights(fr):
	if 'TodWeights' not in fr:
		return
	print('Functioning detectors: ', numpy.sum(numpy.asarray(list(fr['TodWeights'].values())) != 0))
pipe.Add(printweights)

# Kick off maps
pipe.Add(mapmaker.MapInjector, map_id=source + '-150GHz',
  maps_lst=[smstub], is_stub=True, make_polarized=True, do_weight=True)
pipe.Add(mapmaker.CalculatePointing, map_id=source + '-150GHz',
  pointing_store_key = 'PixelPointing', ts_map_key = 'CleanedTS',
  trans_key='OnlineRaDecRotation')
pipe.Add(mapmaker.BinMap, map_id=source + '-150GHz',
  ts_map_key='CleanedTS150GHz',
  pointing_store_key='PixelPointing', timestream_weight_key = 'TodWeights',
  trans_key='OnlineRaDecRotation')

pipe.Add(lambda fr: fr.type == core.G3FrameType.Map) # Drop TOD

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()

if args.accumulate_ts_for_analysis is not None:
    print('Concatenating and downsampling timestreams')
    def downsample(scan, N=100):
        # XXX: This is a terrible downsampler and should be replaced
        data = scan[:scan.size - (scan.size % N):N].copy()
        for i in range(1, N):
                data += scan[i:scan.size - (scan.size % N):N]
        data /= N
        return data
    tss = {i[0]: downsample(numpy.concatenate(i[1])) for i in ts.items()}
    numpy.save(args.accumulate_ts_for_analysis, tss)


