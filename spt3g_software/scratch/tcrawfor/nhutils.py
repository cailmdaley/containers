import os
import datetime
import glob
import re
from collections import deque
import numpy as np
from spt3g import core, mapmaker, gcp, timestreamflagging, todfilter
from spt3g import std_processing, calibration
from spt3g.mapmaker import mapmakerutils as maputils
from spt3g.todfilter import util
from spt3g.pointing import focus

class Accumulator(object):
    '''
    A handy module for grabbing a key (and possibly subkeys) from frames.
    Extracted data is available in Accumulator.values
    '''
    def __init__(self, keys, callback = None, frametype = None):
        '''
        keys: str or list of str
            The keys (and subkeys) you want to extract.  For example,
            keys = ['key', 'subkey1', 'subkey2']
            will result in Accumulator.values containing the value of
            fr['key']['subkey1']['subkey2'] for each frame.
        callback: callable
            `callback` is run on the object extracted from each frame.
            For example, callback = lambda az: az / core.G3Units.deg
            would convert azimuth to degrees before storing it in
            Accumulator.values.
        frametype: core.G3FrameType
            Only extrac objecst from frames of this type.        
        '''
        self.keys = keys
        if isinstance(keys, str):
            self.keys = [keys]
        else:
            self.keys = keys
        self.values = deque()
        if callback is None:
            self.callback = lambda v: v
        else:
            self.callback = callback
        self.frametype = frametype

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        if self.frametype is not None and fr.type != self.frametype:
            return
        try:
            thingy = fr[self.keys[0]]
        except KeyError:
            return fr
        for k in self.keys[1:]:
            try:
                thingy = thingy[k]
            except KeyError:
                return fr
        self.values.append(self.callback(thingy))
        return fr

def get_pointing(files, ptg_type = 'azel', units = 1):
    ptg_keys = {'azel': ['RawBoresightAz', 'RawBoresightEl'],
                'radec': ['BoresightRa', 'BoresightDec'],}
    az_key, el_key = ptg_keys[ptg_type]
    acc_az = Accumulator(az_key, frametype = core.G3FrameType.Scan)
    acc_el = Accumulator(el_key, frametype = core.G3FrameType.Scan)
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename = files)
    pipe.Add(acc_az)
    pipe.Add(acc_el)
    pipe.Run()
    az = np.concatenate(acc_az.values) / units
    el = np.concatenate(acc_el.values) / units
    return az, el

def convert_G3_types(a):
    try:
        if np.iterable(a):
            a = np.array([v.value for v in a])
        else:
            a = a.value
    except AttributeError:
        pass
    return a        

def extract_keys(file, keys, frametype = None, verbose = False, 
                 arc_extract = True):
    '''
    This is a wrapper around Accumulator to extract one or more keys
    from a list of files.

    INPUTS
    ------
    file: str or list
        The filename or list of filenames from which to extract information
    keys: str, list or dict
        If `keys` is a string or list, only that key (or the specified subkey)
        are extracted from the frames.  The output will be the concatenated
        data from the specified key.
        If `keys` is a dict, each key in `keys` should map to a string
        or list of strings that specify the key or subkey(s) to be extracted.
        The return from `exract_keys` will be a dictionary, with the 
        concatenated data from each key.
        For example, if 
        keys = {'az': 'OnlineBoresightAz', 
                'bolo': ['RawTimestreams_I', <boloname>]}
        Then the output will be a dictionary mapped as follows:
        {'az': the azimuth for all frames,
         'bolo': the timestream for <boloname>}
    frametype: core.G3FrameType
        Only extrac objecst from frames of this type.
    verbose: bool
        If True, print every frame
    arc_extract: bool
        If the input files are archive files, and `arc_extract` is True,
        run extra processing (gcp.ARCExtractor) to make the output slightly 
        nicer.

    TODO: respect G3 types instead of dumping into np.array
    '''
    if isinstance(file, str):
        if file.endswith('.g3'):
            arcfile = False
            reader = core.G3Reader(file)
        elif file.endswith('.dat') or file.endswith('.dat.gz'):
            arcfile = True
            reader = gcp.ARCFileReader(file)
        else:
            raise ValueError('Extension not recognized on {}'.format(file))
    elif np.iterable(file):
        if file[0].endswith('.g3'):
            arcfile = False
            reader = core.G3Reader(file)
        elif file[0].endswith('.dat') or file.endswith('.dat.gz'):
            arcfile = True
            reader = gcp.ARCFileReader(file)
        else:
            raise ValueError('Extension not recognized on {}'.format(file[0]))
    else:
        arcfile = False
        reader = file
    if isinstance(keys, dict):
        accums = {}
        for k in keys:
            accums[k] = Accumulator(keys[k], frametype = frametype)
    else:
        accums = {'value': Accumulator(keys, frametype = frametype)}
    pipe = core.G3Pipeline()
    pipe.Add(reader)
    if arcfile and arc_extract:
        pipe.Add(gcp.ARCExtract)
    if verbose:
        pipe.Add(core.Dump, added_message = 'Before accumulate')
    for a in accums.values():
        pipe.Add(a)
    if verbose:
        pipe.Add(core.Dump, added_message = 'After accumulate')
    pipe.Run()
    out = {}
    for k, acc in accums.items():
        if (isinstance(acc.values[0], core.G3Timestream) or 
            isinstance(acc.values[0], core.G3TimestreamMap)):
            arr = type(acc.values[0]).concatenate(acc.values)
        else:
            arr = np.array(acc.values).flatten()
            if arr.dtype == np.dtype('O'):
                for i, v in enumerate(arr):
                    arr[i] = convert_G3_types(v)
            shape = np.shape(acc.values[0])
            if len(shape) > 0 and shape[0] != 1:
                if len(shape) == 1:
                    if shape[0] == 100:
                        out[k] = arr
                        continue
                    else:
                        # arr = arr.reshape((len(arr) // shape[0], shape[0])).T
                        arr = np.concatenate(arr)
                else:
                    arr = arr.reshape((len(arr) // (shape[0] * shape[1]),
                                       shape[1] * shape[0])).T
        out[k] = arr
    if 'value' in out.keys() and len(out.keys()) == 1:
        out = out['value']
    return out

def guess_obsid(txt):
    if txt.endswith('.g3'):
        return int(re.findall('([0-9]+)', txt)[-2])
    return int(re.findall('([0-9]+)', txt)[-1])

def extract_bolos(files, bolos):
    out = {}
    for b in bolos:
        out[b] = deque()
    for file in files:
        print(file)
        file = core.G3File(file)
        for fr in file:
            if fr.type != core.G3FrameType.Scan:
                continue
            for b in bolos:
                out[b].append(fr['RawTimestreams_I'][b])
    for b in bolos:
        out[b] = core.G3Timestream.concatenate(out[b])
    return out

def obsid_to_time(obsid):
    obsid = int(obsid)
    mjd0 = core.G3Time('01-Jan-2017:00:00:00').mjd
    time_out = core.G3Time()
    time_out.mjd = mjd0 + obsid/86400.
    return time_out

def find_cal(obsid, dir = '/spt/user/production/calibration/calibrator/'):
    if isinstance(obsid, str):
        obsid = int(obsid[:-3])
    fobs = glob.glob(os.path.join(dir, '*.g3'))
    calids = []
    for f in fobs:
        match = re.match('(\d+)(\.g3)?', os.path.basename(f))
        calids.append(int(match.group(1)))
    calids = sorted(calids)
    if calids[0] > obsid:
        raise RuntimeError("Can't find calibration before {:d}".format(obsid))
    if calids[-1] < obsid:
        return fobs[-1]
    for i, cal in enumerate(calids):
        if cal > obsid:
            break
    return os.path.join(dir, str(calids[i - 1]) + '.g3')

def trim_bolos(frame, bolos_to_use = [], ts_in_key = 'RawTimestreams_I', 
               ts_out_key = 'TrimmedTimestreams'):
    if ts_in_key not in frame:
        return
    outmap = core.G3TimestreamMap()
    for b in bolos_to_use:
        if b not in frame[ts_in_key]:
            continue
        outmap[b] = frame[ts_in_key][b]
    frame[ts_out_key] = outmap

def map_from_file(filename, map_params, ts_in_key, 
                  bs_x_key = 'BoresightRa', bs_y_key = 'BoresightDec',
                  bolo_props_key = 'BolometerProperties',
                  map_id = 'quickmap', flag_reasons = None, 
                  bolos_to_use = None, 
                  freq_low = 11.0, freq_high = 19.0, **kwargs):
    psd_bins = [(0.0*core.G3Units.hz, 1.0*core.G3Units.hz), 
                (1.0*core.G3Units.hz, 2.0*core.G3Units.hz),
                (2.0*core.G3Units.hz, 3.0*core.G3Units.hz), 
                (3.0*core.G3Units.hz, 6.0*core.G3Units.hz), 
                (6.0*core.G3Units.hz, 9.0*core.G3Units.hz),
                (9.0*core.G3Units.hz, 11.0*core.G3Units.hz),
                (11.0*core.G3Units.hz, 15.0*core.G3Units.hz),
                (15.0*core.G3Units.hz, 19.0*core.G3Units.hz),
                (19.0*core.G3Units.hz, 30.0*core.G3Units.hz),
                (30.0*core.G3Units.hz, 40.0*core.G3Units.hz),
                (40.0*core.G3Units.hz, 50.0*core.G3Units.hz)]

    def replace_key(frame, oldkey = None, newkey = None, 
                    type = core.G3FrameType.Scan):
        if frame.type != type:
            return
        frame.pop(oldkey, None)
        frame[oldkey] = frame.pop(newkey, None)

    map_extractor = maputils.ExtractTheMaps()
    pipe = core.G3Pipeline()
    r = core.G3Reader(filename)
    pipe.Add(core.G3Reader, filename = filename)
    pipe.Add(core.DeduplicateMetadata)
    pipe.Add(util.CutTimestreamsWithoutProperties, map = bolo_props_key,
             input = ts_in_key, output = ts_in_key + '_trimmed')
    pipe.Add(replace_key, oldkey = ts_in_key, newkey = ts_in_key + '_trimmed')
    if bolos_to_use is not None:
        pipe.Add(trim_bolos, bolos_to_use = bolos_to_use, 
                 ts_out_key = ts_in_key + '_trimmed')
        pipe.Add(replace_key, oldkey = ts_in_key, newkey = ts_in_key + '_trimmed')
    pipe.Add(lambda fr: not ('Turnaround' in fr and fr['Turnaround']))
    # pipe.Add(std_processing.ScanPreprocessing.PreprocessNoiseInfo, 
    #          ts_map_key = ts_in_key,
    #          psd_bins = psd_bins,
    #          prefilter_poly_order = 1,
    #          deriv_presmooth_scale = 10,
    #          glitch_thresholds = [10, 20],
    #          glitch_filter_scale = 11,
    #          max_derivative_key = 'TimestreamMaxDerivative',
    #          variance_key = 'TimestreamVariance',
    #          mad_variance_key = 'TimestreamMadVariance',
    #          glitch_num_key = 'GlitchesNumberOf',
    #          glitch_thresholds_key = 'GlitchesThresholds')
    # pipe.Add(timestreamflagging.flaggingutils.GroupAverageG3MapValue,
    #          per_band = True,
    #          per_wafer = True,
    #          per_squid = True,
    #          input_g3_map_key = 'TimestreamMaxDerivative',
    #          output_g3_map_key = 'AvMaxTimestreamDerivative',
    #          average_func = np.median)
    # pipe.Add(timestreamflagging.glitchfinding.FlagGlitches,
    #          min_num_glitch_map = {float(10):1},
    #          flag_key = 'Flags')
    # pipe.Add(timestreamflagging.flaggingutils.FlagGroupAverageG3MapValue,
    #          g3_map_key = 'AvMaxTimestreamDerivative',
    #          min_val = 0,
    #          max_val = 3e-8,
    #          per_band = True,
    #          per_wafer = True,
    #          per_squid = True,
    #          flag_reason = 'BadMedianDeriv')
    # pipe.Add(timestreamflagging.flaggingutils.FlagBadG3MapValue,
    #          flag_key = 'Flags',
    #          m_key = 'TimestreamVariance',
    #          flag_reason = 'BadRmsNoise',
    #          min_val = 1e-7,
    #          max_val = 5e0)
    # pipe.Add(timestreamflagging.flaggingutils.FlagBadG3MapValue,
    #          flag_key = 'Flags',
    #          m_key = 'TimestreamMaxDerivative',
    #          flag_reason = 'BadDeriv',
    #          min_val = 1e-11,
    #          max_val = 1e-7)
    if flag_reasons is not None:
        pipe.Add(get_flag_reasons, flag_reasons = flag_reasons)
    # pipe.Add(timestreamflagging.flaggingutils.RemoveFlaggedTimestreams,
    #          input_ts_key = ts_in_key,
    #          input_flag_key = 'Flags',
    #          output_ts_key = 'FlaggedTs')
    pipe.Add(todfilter.dftutils.AddBinnedPsdInfo, 
             ts_map_key = ts_in_key, bins = psd_bins)
    pipe.Add(less_janky_weights, ts_key = ts_in_key,
             freq_low = freq_low * core.G3Units.hz, 
             freq_high = freq_high * core.G3Units.hz)
    pipe.Add(janky_weights, ts_key = ts_in_key, ts_weight_key = 'JankyWeights')
    # pipe.Add(core.Dump)
    outfile = '/poleanalysis/ndhuang/test_frames/focus_{:.02f}_{:.02f}.g3'.format(freq_low, freq_high)
    if not os.path.exists(outfile):
        pipe.Add(core.G3Writer, filename = outfile)
    pipe.Add(maputils.MakeMap, map_in = map_params, map_id = map_id,
             poly_order = 4,
             ts_in_key = ts_in_key,
             boresight_ra_key = bs_x_key, 
             boresight_dec_key = bs_y_key,
             bolo_props_key = bolo_props_key,
             timestream_weight_key = 'JankyWeights',
             **kwargs)
    pipe.Add(maputils.RemoveWeightModule)
    pipe.Add(map_extractor)
    pipe.Run(profile = True)
    return map_extractor.maps[map_id]

# def inspect(frame, type = core.G3FrameType.Scan):
#     if frame.type == type:
#         import IPython
#         IPython.embed()

def less_janky_weights(fr, ts_key = 'CalTimestreams', 
                       freq_low = 11.0 * core.G3Units.hz, 
                       freq_high = 19.0 * core.G3Units.hz,
                       psd_key = "BinAveragedPsds",
                       psd_lower_bound_key = "BinAveragedPsdsLowerBound",
                       psd_upper_bound_key = "BinAveragedPsdsUpperBound",
                       ts_weight_key = "TodWeights"):
    if fr.type != core.G3FrameType.Scan:
        return
    psd_lower = fr[psd_lower_bound_key]
    psd_upper = fr[psd_upper_bound_key]
    psd_inds = np.where((np.array(psd_lower) >= freq_low) & 
                        (np.array(psd_upper) <= freq_high))[0]
    total_bandwidth = 0
    noise = core.G3MapDouble()
    for i in psd_inds:
        band = psd_upper[i] - psd_lower[i]
        total_bandwidth += band
        for bolo in fr[ts_key].keys():
            if bolo not in noise.keys():
                noise[bolo] = band * fr[psd_key][bolo][i]
            else:
                noise[bolo] += band * fr[psd_key][bolo][i]
    weights = core.G3MapDouble()
    for bolo in fr[ts_key].keys():
        if noise[bolo] == 0 or not np.isfinite(noise[bolo]):
            weights[bolo] = 0
        else:
            weights[bolo] = 1 / noise[bolo] * total_bandwidth
    # print np.shape(np.nonzero(weights.values()))
    fr[ts_weight_key] = weights
    return fr

def janky_weights(fr, ts_key = 'CalTimestreams', 
                  ts_weight_key = 'TodWeights'):
    if fr.type != core.G3FrameType.Scan:
        return
    weights = core.G3MapDouble()
    for bolo in fr[ts_key].keys():
        weight = 1 / np.var(fr[ts_key][bolo])
        if not np.isfinite(weight):
            weight = 0
        weights[bolo] = weight
    fr[ts_weight_key] = weights
    
def get_flag_reasons(frame, flag_reasons = None):
    try:
        flag_reasons.append(frame['Flags'])
    except KeyError:
        pass
             
def fp_average(fr, ts_key):
    m = np.zeros_like(fr[ts_key].times())
    for bolo, ts in fr[ts_key]:
        m += ts
    return m / len(fr[ts_key])

def flatten_inhomo(arr):
    lens = [len(a) for a in arr]
    out = np.empty(np.sum(lens))
    i = 0
    for a in arr:
        out[i: i + len(a)] = a
        i += len(a)
    return out

def fix_off_source(off):
    '''
    Translate whacky GCP strings back to bools.  Attempts to make a copy
    '''
    off = off.copy()
    for i, thing in enumerate(off):
        if isinstance(thing, str):
            off[i] = np.zeros(100)
            for j, b in enumerate(thing.encode()):
                off[i][j] = b
        if isinstance(thing, core.G3VectorInt):
            off[i] = np.array(thing)
    return np.concatenate(off)

def split_bolo_list_band(bolos, boloprops, index = True):
    bolos_by_band = {}
    for i, b in enumerate(bolos):
        if index:
            value = i
        else:
            value = b
        try:
            bolos_by_band[int(boloprops[b].band / core.G3Units.GHz)].append(value)
        except KeyError:
            bolos_by_band[int(boloprops[b].band / core.G3Units.GHz)] = [value]
    return bolos_by_band

def get_key_by_band(fr, key, boloprops = None):
    if boloprops is None:
        boloprops = fr['BolometerProperties']
    bolos_by_band = split_bolo_list_band(fr[key].keys(), boloprops, False)
    out = {}
    for band, bolos in bolos_by_band.items():
        out[band] = np.array([fr[key][b] for b in bolos])
    return out

def inject_dt(fr):
    if fr.type != core.G3FrameType.Scan:
        return
    t_ts = fr['RawTimestreams_I'].times()
    dt = core.G3VectorDouble([t1.time - t2.time for 
                              t1, t2 in zip(t_ts, fr['DetectorSampleTimes'])])
    fr['TimeDiff'] = dt

def check_time_diffs(files):
    a = Accumulator('TimeDiff')
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename = files)
    pipe.Add(inject_dt)
    pipe.Add(a)
    pipe.Run()
    out = np.concatenate([np.array(v) for v in a.values])
    out = out.flatten()
    return out

def get_focus_position(obsid, bolodata = '/spt_data/bolodata/fullrate/', 
                       source = '0537-441-pixelraster'):
    f = core.G3File(os.path.join(bolodata, source, str(obsid), '0000.g3'))
    for fr in f:
        if 'BenchCommandedPosition' in fr.keys():
            break
    benchpos = focus.g3bench2xyz(fr['BenchCommandedPosition'])
    opticspos = focus.bench2optical(*benchpos)
    return np.array(opticspos).flatten()

def arcfiles_for_time(start, stop, arcdir = '/spt_data/arc'):
    arc2time = lambda path: core.G3Time(os.path.basename(path).split('.')[0])
    arcfiles = [std_processing.FindARCFiles.ARCForTime(start, arcdir)]
    while arc2time(arcfiles[-1]).time < stop.time:
        arcfiles.append(std_processing.FindARCFiles.NextARCFile(arcfiles[-1]))
    return arcfiles[:-1]

def recent_obsids(dir, dur = 7 * core.G3Units.days):
    obsids = []
    for p in os.listdir(dir):
        try:
            obsids.append(guess_obsid(p))
        except IndexError:
            continue
    obsids = sorted(obsids)
    times = [obsid_to_time(oid) for oid in obsids]
    last = core.G3Time.Now() - dur
    last_i = np.searchsorted(times, last)
    return obsids[last_i:]

def obsid_to_datetime(obsid):
    return (datetime.datetime(2017, 1, 1, 0, 0, 0, 0) + 
            datetime.timedelta(seconds = obsid))
