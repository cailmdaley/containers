import os
import glob
import re
from collections import deque
import numpy as np
from spt3g import core, mapmaker, gcp, timestreamflagging
from spt3g import std_processing, calibration
from spt3g.mapmaker import mapmakerutils as maputils
from spt3g.todfilter import util

class Accumulator(object):
    def __init__(self, keys, callback = None, frametype = None):
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

    def get_ret(self):
        return np.concatenate(self.values)

    @staticmethod
    def getval(fr, keys):
        for k in keys:
            fr_k = fr[k]
        try:
            if isinstance(fr_k, core.G3Timestream):
                fr_k = np.array(fr_k)
            if np.iterable(fr_k):
                fr_k = np.array([v.value for v in fr_k])
            else:
                fr_k = fr_k.value
        except AttributeError:
            pass
        try:
            if np.iterable(fr_k):
                fr_k = np.array([v.time for v in fr_k])
            else:
                fr_k = fr_k.time
        except AttributeError:
            pass
        return fr_k

    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        if self.frametype is not None and fr.type != self.frametype:
            return
        thingy = fr[self.keys[0]]
        for k in self.keys[1:]:
            thingy = thingy[k]
        self.values.append(thingy)
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

def extract_keys(file, keys, frametype = None):
    if isinstance(file, str):
        file = core.G3File(file)
    if isinstance(keys, str):
        keys = [keys]
    out = {}
    for k in keys:
        out[k] = []
    for fr in file:
        if frametype is not None and frametype != fr.type:
            continue
        for k in keys:
            out[k].append(fr[k])
    for k in keys:
        n = sum([len(a) for a in out[k]])
        vec = np.zeros(n)
        i = 0
        for a in out[k]:
            vec[i: i + len(a)] = a
            i += len(a)
        out[k] = vec
    return out

def find_cal(obsid, dir):
    fobs = sorted(glob.glob(os.path.join(dir, '*')))
    calids = []
    for f in fobs:
        match = re.match('(\d+)(\.g3)?', os.path.basename(f))
        calids.append(int(match.group(1)))
    if calids[0] > obsid:
        raise RuntimeError("Can't find calibration before {:d}".format(obsid))
    if calids[-1] < obsid:
        return fobs[-1]
    for i, cal in enumerate(calids):
        if cal > obsid:
            break
    return fobs[i - 1]

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
                  bolos_to_use = None, **kwargs):
    psd_bins = [(0.0*core.G3Units.hz, 1.0*core.G3Units.hz), 
                (1.0*core.G3Units.hz, 2.0*core.G3Units.hz),
                (2.0*core.G3Units.hz, 3.0*core.G3Units.hz), 
                (3.0*core.G3Units.hz, 6.0*core.G3Units.hz), 
                (6.0*core.G3Units.hz, 9.0*core.G3Units.hz),
                (9.0*core.G3Units.hz, 12.0*core.G3Units.hz),
                (12.0*core.G3Units.hz, 15.0*core.G3Units.hz),
                (15.0*core.G3Units.hz, 20.0*core.G3Units.hz),
                (20.0*core.G3Units.hz, 30.0*core.G3Units.hz),
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
    pipe.Add(janky_weights, ts_key = ts_in_key)
    pipe.Add(core.Dump)
    # pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
    pipe.Add(maputils.MakeMap, map_in = map_params, map_id = map_id,
             poly_order = 4,
             ts_in_key = ts_in_key,
             boresight_ra_key = bs_x_key, 
             boresight_dec_key = bs_y_key,
             bolo_props_key = bolo_props_key,
             timestream_weight_key = 'TodWeights',
             **kwargs)
    pipe.Add(maputils.RemoveWeightModule)
    pipe.Add(map_extractor)
    pipe.Run()
    return map_extractor.maps[map_id]

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
    
def delete_mapping_cache(fr):
    for k in ['MapPointing', 'DetectorAlphaPointing', 'DetectorDeltaPointing']:
        # also stokes
        fr.pop(k, None)

def unpack_pointing_from_arc(files):
    az = Accumulator(['antenna0', 'tracker', 'actual', 0])
    el = Accumulator(['antenna0', 'tracker', 'actual', 1])
    azr = Accumulator(['antenna0', 'tracker', 'actual_rates', 0])
    elr = Accumulator(['antenna0', 'tracker', 'actual_rates', 1])
    azerr = Accumulator(['antenna0', 'tracker', 'errors', 0])
    elerr = Accumulator(['antenna0', 'tracker', 'errors', 1])
    cryo = Accumulator(['array', 'cryo', 'temperature', 0, 0])
    accums = [az, el, azr, elr, azerr, elerr, cryo]

    pipe = core.G3Pipeline()
    if isinstance(files, str):
        pipe.Add(gcp.ARCFileReader, filename = files)
    else:
        pipe.Add(gcp.ARCFileReader, filenames = files)
    for acc in accums:
        pipe.Add(acc)
    pipe.Run()
    arrs = [np.array(acc.values).flatten() / core.G3Units.deg for acc in accums]
    arrs[-1] *= core.G3Units.deg
    return arrs

