import os
import datetime
import glob
import itertools
import re
from collections import namedtuple
import numpy as np
from scipy import io
from spt3g import core, mapmaker, gcp, timestreamflagging, todfilter
from spt3g import std_processing, calibration
from spt3g.mapmaker import mapmakerutils as maputils
from spt3g.todfilter import util
from spt3g.pointing import focus
from spt3g.util.extractdata import extract_keys, Accumulator, convert_G3_types, MultiAccumulator

anttrack = ['antenna0', 'tracker']
keys = {'pointing': {'az': anttrack + ['actual', 0],
                     'el': anttrack + ['actual', 1],
                     'azerr': anttrack + ['errors', 0],
                     'elerr': anttrack + ['errors', 1],
                     'azexp': anttrack + ['expected', 0],
                     'elexp': anttrack + ['expected', 1],
                     'time': anttrack + ['utc']},
        'elencoders': {0: anttrack + ['raw_encoder', 0],
                       1: anttrack + ['raw_encoder', 1],
                       'el': anttrack + ['actual', 1],
                       'time': anttrack + ['utc']},
        'scanification': {'scanflags': anttrack + ['scan_flag'],
                          'fbits': ['antenna0', 'frame', 'features'],
                          'offsource': anttrack + ['off_source']}
    }

def extract_arctime(start, stop, keys, arcdir = '/spt/data/arc'):
    if isinstance(start, str):
        start = core.G3Time(start.strip())
    if isinstance(stop, str):
        stop = core.G3Time(stop.strip())
    if isinstance(start, datetime.datetime):
        start = core.G3Time(start.strftime('%Y%m%d_%H%M%S'))
    if isinstance(stop, datetime.datetime):
        stop = core.G3Time(stop.strftime('%Y%m%d_%H%M%S'))
    pipe = core.G3Pipeline()
    macc = MultiAccumulator(keys = keys)
    pipe.Add(std_processing.ARCTimerangeReader, start_time = start, 
             stop_time = stop, basedir = arcdir)
    pipe.Add(gcp.ARCExtract)
    pipe.Add(macc)
    pipe.Run()
    return macc.extract_values()

def guess_obsid(txt):
    if txt.endswith('.g3'):
        return int(re.findall('([0-9]+)', txt)[-2])
    return int(re.findall('([0-9]+)', txt)[-1])

def obsid_to_time(obsid):
    obsid = int(obsid)
    mjd0 = core.G3Time('01-Jan-2017:00:00:00').mjd
    time_out = core.G3Time()
    time_out.mjd = mjd0 + obsid/86400.
    return time_out

def find_cal(obsid, dir = '/spt/user/production/calibration/calibrator/'):
    if not os.path.exists(dir):
        raise ValueError('{} does not exist'.format(dir))
    if isinstance(obsid, str):
        obsid = int(obsid[:-3])
    fobs = sorted(glob.glob(os.path.join(dir, '*.g3')))
    if len(fobs) == 0:
        fobs = os.listdir(dir)
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
    # this return is wrong if we're looking for bolodata dirs
    return os.path.join(dir, str(calids[i - 1]) + '.g3')

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
    return std_processing.FindARCFiles.ARCForTimerange(start, stop, arcdir)

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

def unwrap(a, rad = True):
    if rad:
        split = np.pi
    else:
        split = 180
    a[a > split] -= 2 * split

def wrap_ra(x, new_center = 0):
    x[x > new_center + 180 % 360] -= 360
    return x

def at20_sources_for_field(radec0, extent, buffer = 0, flux_cut = 300,
                           catalogfile = '/home/ndhuang/spt_code/spt_analysis/idlsave/at20g.sav',
                           verbose = False):
    # speaks degrees, not G3 units
    at20g = io.readsav(catalogfile)['at20g']
    ra_extent = [radec0[0] - extent[0] / 2. + buffer,
                 radec0[0] + extent[0] / 2. - buffer]
    dec_extent = [radec0[1] - extent[1] / 2. + buffer,
                  radec0[1] + extent[1] / 2. - buffer]
    if ra_extent[0] < 0:
        at20g.ra = wrap_ra(at20g.ra)
    if verbose:
        print(ra_extent)
        print(dec_extent)
    whinmap = np.where((at20g.ra < ra_extent[1]) & 
                       (at20g.ra > ra_extent[0]) &
                       (at20g.dec < dec_extent[1]) &
                       (at20g.dec > dec_extent[0]) &
                       (at20g.flux > flux_cut))
    return at20g[whinmap]

def split_map_by_band(map, calfile):
    calframe = core.G3File(calfile).next()
    if 'BolometerProperties' in calframe:
        bp = calframe['BolometerProperties']
    else:
        bp = calframe['NominalBolometerProperties']

    byband = {}
    for b in map.keys():
        band = bp[b].band
        if band not in byband.keys():
            byband[band] = [b]
        else:
            byband[band].append(b)
    map_by_band = {}
    for band in byband:
        map_by_band = np.array([map[b] for b in byband[band]])
    return map_by_band, byband

def getbp(fname):
    try:
        return core.G3File(fname).next()['BolometerProperties']
    except KeyError:
        fr = core.G3File(fname).next()
        x = fr['PointingOffsetX']
        y = fr['PointingOffsetY']
        fakeBP = namedtuple('fakeBP', ['x_offset', 'y_offset'])
        bp = {b:fakeBP(x[b], y[b]) for b in x.keys()}
        return bp

def compare_bp(fname1, fname2):
    bp1 = getbp(fname1)
    bp2 = getbp(fname2)
    bolos = np.array(tuple(set(bp1.keys()).intersection(set(bp2.keys()))))
    nbolos = len(bolos)
    dx = np.zeros(nbolos)
    dy = np.zeros(nbolos)
    for i, bolo in enumerate(bolos):
        dx[i] = (bp1[bolo].x_offset - bp2[bolo].x_offset)# / bp1[bolo].x_offset
        dy[i] = (bp1[bolo].y_offset - bp2[bolo].y_offset)# / bp1[bolo].y_offset
    good = np.isfinite(dx)
    return dx[good], dy[good], bolos[good]

class bpcomparator(object):
    def __init__(self):
        self.bp = None
        self.fakebp = None
        
    def __call__(self, fr):
        if fr.type == core.G3FrameType.EndProcessing:
            return
        if 'BolometerProperties' in fr:
            self.bp = fr['BolometerProperties']
        bp = getbp('/spt_data/bolodata/fullrate/PMNJ0210-5101-pixelraster/47078333/offline_calibration.g3')
        for b in bp.keys():
            if np.isfinite(bp[b].x_offset):
                try:
                    assert((bp[b].x_offset == self.bp[b].x_offset) &
                           (bp[b].y_offset == self.bp[b].y_offset))
                except AssertionError:
                    print(bp[b].x_offset, self.bp[b].x_offset, '\n', 
                          bp[b].y_offset, self.bp[b].y_offset)
        '''
        if 'PointingOffsetX' in fr:
            x = fr['PointingOffsetX']
            y = fr['PointingOffsetY']
            fakeBP = namedtuple('fakeBP', ['x_offset', 'y_offset'])
            self.fakebp = {b:fakeBP(x[b], y[b]) for b in x.keys()}
        if self.bp is not None and self.fakebp is not None:
            print('doing')
            bp1 = self.bp
            bp2 = self.fakebp
            bolos = np.array(tuple(set(bp1.keys()).intersection(set(bp2.keys()))))
            nbolos = len(bolos)
            dx = np.zeros(nbolos)
            dy = np.zeros(nbolos)
            for i, bolo in enumerate(bolos):
                dx[i] = (bp1[bolo].x_offset - bp2[bolo].x_offset) / bp1[bolo].x_offset
                dy[i] = (bp1[bolo].y_offset - bp2[bolo].y_offset) / bp1[bolo].y_offset
            good = np.isfinite(dx)
            print('making file')
            np.savez('/poleanalysis/ndhuang/bp_check_wtf.npz', dx = dx[good], dy = dy[good], bolos = bolos[good])
        '''
def comparepixpointing(fname1, fname2):
    out = {}
    f2 = core.G3File(fname2)
    for fr in core.G3File(fname1):
        frs = (fr, f2.next())
        pps = [fr['PixelPointing'] for fr in frs]
        for bolo in pps[0].keys():
            if bolo not in pps[1]:
                continue
            if bolo not in out:
                out[bolo] = [np.mean(np.array(pps[0][bolo]) - 
                                     np.array(pps[1][bolo]))]
            else:
                out[bolo].append(np.mean(np.array(pps[0][bolo]) - 
                                         np.array(pps[1][bolo])))
    return out

def reconstruct_el(enc1, enc2, off1 = 168.884, off2 = 94.6230):
    # Someone has done something evil with enc2.
    # Based on the readout on the ACU, el2 = off2 - enc2
    # However, in the gcp register, el2 = off2 - (360 - enc2)
    el = np.array([enc1 - off1, off2 - (360 - enc2)])
    el[np.where(el < 0)] += 360
    return np.mean(el, axis = 0)

def amy_check():
    files = sorted(glob.glob('/poleanalysis/sptdaq/calresult/calibration/calibrator/62*.g3'))
    boloids = ['2019.hv2', '2019.a8g']
    for f in files:
        oid = guess_obsid(f)
        if oid < 62242600:
            continue
        calsn = core.G3File(f).next()['CalibratorResponseSN']
        print(f)
        for b in boloids:
            if b in calsn:
                print(b, calsn[b])
    
# def make_bitarr(intarr, nbits = 64):
#     def uint_to_bin(uint, nbits):
#         out = np.zeros(
#     # Currently assumes uints
#     bitarr = np.zeros((len(intarr), nbits), dtype = np.bool)
    
        
def fbits_to_flags(fbits):
    fbits = fbits.astype(np.uint32)
    features = {0: 'analyze', 1: 'on_source', 2: 'shutter', 3: 'elnod', 
                4: 'cross', 5: 'calibrator', 6: 'raster', 7: 'skydip', 
                8: 'optical', 9: 'noise', 10: 'leadtrail', 11: 'elscan',
                19: 'debug'}
    out = {feat: fbits & 2**bit for bit, feat in features.items()}
    return out

def fbits_to_whatdoing(fbits):
    fbits = fbits.astype(np.uint32)
    features = {0: 'nothing', 1: 'analyze', 8: 'elnod', 32: 'calibrator',
                512: 'noise'}
    keys = np.array(list(features.keys()), dtype = np.uint)
    keyfunc = lambda k: 0 if not k else 2**int(np.log2(k))
    out = list(itertools.groupby(fbits, keyfunc))
    out = np.array([o[0] for o in out])
    out[out == 0] = -1
    out[out > 0] = np.log2(out[out > 0])
    out = out.astype(int)
    return out

