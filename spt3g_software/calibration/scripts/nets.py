description = '''Script to calculate NETs from some random data. Noise is taken from power in some band (3.0 to 5.0 Hz by default).  NETs are in uK.'''
import sys
import argparse
try:
    import cPickle as pickle
except ImportError:
    import pickle
import scipy.stats
import numpy as np
from spt3g import core, std_processing, calibration, xtalk
from spt3g import todfilter, mapmaker, timestreamflagging

parser = argparse.ArgumentParser(description = description)
parser.add_argument('calframe', 
                      help = 'Calibration frame for the observation')
parser.add_argument('timestreams', nargs = '+',
                      help = 'Detector timestreams to measure NET')
parser.add_argument('--lowf', type = float, default = 3.0,
                    help = 'The low end of the frequency band')
parser.add_argument('--highf', type = float, default = 5.0,
                    help = 'The high end of the frequency band')
parser.add_argument('--output', '-o', default = 'nets.pkl',
                  help = 'Output filename')
parser.add_argument('--verbose', '-v', help = 'Print all frames', 
                    action = 'store_true', default = False)
parser.add_argument('--keep-turnarounds', action = 'store_true', 
                    default = False, 
                    help = 'Do not drop turnaround frames.  Useful for noise stares')
parser.add_argument('--plot', action = 'store_true', default = False,
                    help = 'Plot a histogram of NETs by band at the end')
parser.add_argument('--num-frames', type = int, default = 4,
                    help = 'Maximum number of Scan frames to process.  0 means no limit')
args = parser.parse_args()

def weight_to_net(fr, weight_key = 'TodWeights', low_f = 3.0,
                    high_f = 5.0):
    '''
    Convert a bolometer weight (inverse integrated noise power) to NET
    '''
    if weight_key not in fr:
        return
    net = core.G3MapDouble()
    bandwidth = high_f - low_f
    for bolo in fr[weight_key].keys():
        net[bolo] = 1. / np.sqrt(bandwidth * fr[weight_key][bolo] * 2) * 1e6
    fr['NET'] = net

class NETGobbler(object):
    '''
    Collect NETs across multiple scans
    '''
    def __init__(self, net_key = 'NET'):
        self.net = dict()
        self.input = net_key
        self.boloprops = None

    def __call__(self, fr):
        if 'BolometerProperties' in fr:
            self.boloprops = fr['BolometerProperties']
        if self.input not in fr:
            return
        for bolo in fr[self.input].keys():
            if bolo not in self.net:
                self.net[bolo] = set()
            self.net[bolo].add(fr[self.input][bolo])
            
    def get_final_net(self):
        '''
        Average (median) all NETs and return a dictionary with the results
        '''
        out = {}
        for bolo in self.net.keys():
            out[bolo] = (np.median([n for n in self.net[bolo]]), 
                         self.boloprops[bolo].band)
        return out

class FrameLimiter(object):
    def __init__(self, maxframes = 4):
        self.maxframes = maxframes
        self.nframes = 0

    def __call__(self, fr):
        if self.nframes >= self.maxframes:
            core.G3Pipeline.halt_processing()
        if fr.type == core.G3FrameType.Scan:
            self.nframes += 1
            return

if not isinstance(args.timestreams, list):
    inputfiles = [args.calframe, args.timestreams]
else:
    inputfiles = [args.calframe,] + args.timestreams
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = inputfiles)
# Cut turnarounds
if not args.keep_turnarounds:
    pipe.Add(lambda fr: 'Turnaround' not in fr or not fr['Turnaround'])
# Limit the number of Scan frames we use
if args.num_frames > 0:
    pipe.Add(FrameLimiter, maxframes = args.num_frames)
# Cut uncalibratable detectors early
pipe.Add(calibration.elnod_analysis.RotateIQ, 
         i_rotate_key='RawTimestreamsRotated')
pipe.Add(todfilter.util.CutTimestreamsWithoutProperties, 
         input='RawTimestreamsRotated', output='FilteredTimestreams')
# Invert crosstalk and convert to Watts.
pipe.Add(xtalk.CrosstalkCleanedTimestreams, input='FilteredTimestreams', 
         output='TimestreamsWatts', units=core.G3TimestreamUnits.Power,
         ignore_missing=True)
# Next to source-relative units (XXX: hardcode RCW38 here)
pipe.Add(calibration.ApplyTCalibration, Input='TimestreamsWatts', 
         Output='CalTimestreams')
# Basic timestream filtering
pipe.Add(mapmaker.TodFiltering, ts_in_key='CalTimestreams',
    ts_out_key='PolyFilteredTimestreams', use_dynamic_source_filter=True,
    poly_order=4)
# Clean up detritus
pipe.Add(core.Delete, keys=['RawTimestreams_I', 'RawTimestreams_Q', 
                            'FilteredTimestreams', 'TimestreamsWatts', 
                            'CalTimestreams'])
# Flag timestreams for calibrator SN and glitches
# Deglitching will also remove real sources
pipe.Add(std_processing.flagsegments.FlagNonResponsive,
         flag_key = 'Flags')
pipe.Add(timestreamflagging.glitchfinding.AddNumGlitches, thresholds = [5, 10],
         input_ts_key = 'PolyFilteredTimestreams')
pipe.Add(timestreamflagging.glitchfinding.FlagGlitches, 
         min_num_glitch_map = {5: 3, 10:1})
pipe.Add(timestreamflagging.RemoveFlaggedTimestreams, 
         input_ts_key = 'PolyFilteredTimestreams',
         input_flag_key = 'Flags',
         output_ts_key = 'FlaggedTimestreams')
# We use TodWeights for NETs, since the AddPSDWeights module
# calculates almost exactly what we want.
pipe.Add(std_processing.weighting.AddPSDWeights, 
         low_f = args.lowf * core.G3Units.Hz,
         high_f = args.highf* core.G3Units.Hz, 
         input = 'FlaggedTimestreams')
pipe.Add(weight_to_net, low_f = args.lowf, high_f = args.highf)
# Grab the NETs from each frame
ng = NETGobbler()
pipe.Add(ng)
if args.verbose:
    pipe.Add(core.Dump)
pipe.Run(profile = args.verbose)

nets = ng.get_final_net()
if len(nets.keys()) == 0:
    print('No bolometers passed cuts!')
    sys.exit()
bands = [res[1] for res in nets.values()]
nets_arr = np.array([res[0] for res in nets.values()])
bolos = np.array(list(nets.keys()))
nets_by_band = {}
for band in np.unique(bands):
    inds = np.where(bands == band)[0]
    nets_by_band[int(band / core.G3Units.GHz)] = (nets_arr[inds], bolos[inds])
with open(args.output, 'wb') as f:
    pickle.dump(nets_by_band, f, protocol = 2, fix_imports = True)

if args.plot:
    from matplotlib import pyplot as pl
    pl.figure()
    bins = np.linspace(0, 4000, num = 51)
    for b in sorted(nets_by_band.keys()):
        nets, bolos = nets_by_band[b]
        inds = np.where(np.isfinite(nets))
        pl.hist(nets[inds], bins = bins, alpha = .6,
                label = '{:d} GHz'.format(b), histtype = 'stepfilled')
    pl.legend()
    pl.title('NETs by band')
    pl.xlabel('NET ($\mu K \sqrt{s}$)')
    pl.ylabel('Number of bolometers')
    pl.show()
