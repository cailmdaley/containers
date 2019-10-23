import numpy, sys, os
from spt3g import core, dfmux, calibration, std_processing, todfilter, mapmaker, util
import argparse as ap

P = ap.ArgumentParser(description='Analyze calibrator data',
                      formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', default='output.g3', action='store',
               help='Output filename')
P.add_argument('--no-nei', action='store_true',
               help='Do not calculate NEI.')
P.add_argument('--no-nep', action='store_true',
               help='Do not calculate NEP.')
P.add_argument('--no-net', action='store_true',
               help='Do not calculate NET.')
P.add_argument('--n-frames', action='store', type=int, default=None,
               help='Calculate noise for the specified number of frames in '
               'the raw data file.')
P.add_argument('--no-turnarounds', action='store_true',
               help='Do not calculate noise for turnaround frames.')
P.add_argument('--cal-sn-threshold', action='store', default=None, type=float,
               help='Minimum calibrator S/N of bolometers to use in NET '
               'calculation. Default is no threshold.')
P.add_argument('--no-opacity-correction', action='store_true',
               help='Do not apply the opacity correction when converting to '
               'units of Kcmb.')
args = P.parse_args()


class StopAfterNFrames(object):
    def __init__(self, n_frames, type_to_skip=None):
        self.j_frame = 0
        self.n_frames = n_frames
        self.type_to_skip = type_to_skip

    def __call__(self, frame):
        if self.type_to_skip != None and frame.type == self.type_to_skip:
            self.j_frame += 1
        if self.type_to_skip == None:
            self.j_frame += 1

        if self.n_frames != None and self.j_frame > self.n_frames:
            core.G3Pipeline.halt_processing()


pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=args.input)

pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames)

if args.n_frames != None:
    pipe.Add(StopAfterNFrames, n_frames=args.n_frames,
             type_to_skip=core.G3FrameType.Scan)

if args.no_turnarounds == True:
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])

pipe.Add(core.Delete, keys=['RawTimestreams_Q'])
pipe.Add(core.Dump)
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties,
         input='RawTimestreams_I', output='TimestreamsWithProperties')
pipe.Add(core.Delete, keys=['RawTimestreams_I'])

pipe.Add(mapmaker.TodFiltering, ts_in_key='TimestreamsWithProperties',
         ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=False,
         poly_order=5)
pipe.Add(core.Delete, keys=['TimestreamsWithProperties'])

if args.no_nei == False:
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='PolyFilteredTimestreams',
              Output='TimestreamsAmps', Units=core.G3TimestreamUnits.Current)
    pipe.Add(calibration.noise_analysis.AnalyzeNoise,
             freq_ranges=[(0.1 * core.G3Units.Hz, 0.5 * core.G3Units.Hz),
                          (1 * core.G3Units.Hz, 2 * core.G3Units.Hz),
                          (3 * core.G3Units.Hz, 5 * core.G3Units.Hz),
                          (10 * core.G3Units.Hz, 15 * core.G3Units.Hz),
                          (30 * core.G3Units.Hz, 40 * core.G3Units.Hz)],
             noise_types=['NEI'],
             nei_input='TimestreamsAmps')
    pipe.Add(core.Delete, keys=['TimestreamsAmps'])

if args.no_nep == False or args.no_net == False:
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='PolyFilteredTimestreams',
              Output='TimestreamsWatts', Units=core.G3TimestreamUnits.Power)

if args.no_nep == False:
    pipe.Add(calibration.noise_analysis.AnalyzeNoise,
             freq_ranges=[(0.1 * core.G3Units.Hz, 0.5 * core.G3Units.Hz),
                          (1 * core.G3Units.Hz, 2 * core.G3Units.Hz),
                          (3 * core.G3Units.Hz, 5 * core.G3Units.Hz),
                          (10 * core.G3Units.Hz, 15 * core.G3Units.Hz),
                          (30 * core.G3Units.Hz, 40 * core.G3Units.Hz)],
             noise_types=['NEP'],
             nep_input='TimestreamsWatts')
    pipe.Add(core.Delete, keys=['PolyFilteredTimestreams'])

if args.no_net == False:
    pipe.Add(calibration.ApplyTCalibration, Input='TimestreamsWatts',
             Output='TimestreamsKcmb', Source=['RCW38', 'MAT5A'],
             OpacityCorrection=(not args.no_opacity_correction))
    pipe.Add(calibration.noise_analysis.AnalyzeNoise,
             freq_ranges=[(0.1 * core.G3Units.Hz, 0.5 * core.G3Units.Hz),
                          (1 * core.G3Units.Hz, 2 * core.G3Units.Hz),
                          (3 * core.G3Units.Hz, 5 * core.G3Units.Hz),
                          (10 * core.G3Units.Hz, 15 * core.G3Units.Hz),
                          (30 * core.G3Units.Hz, 40 * core.G3Units.Hz)],
             noise_types=['NET'],
             net_input='TimestreamsKcmb',
             min_cal_sn=args.cal_sn_threshold)
    pipe.Add(core.Delete, keys=['TimestreamsKcmb'])

pipe.Add(calibration.noise_analysis.MakeNoiseFrame)

pipe.Add(core.G3Writer, filename=args.output)
pipe.Run()
