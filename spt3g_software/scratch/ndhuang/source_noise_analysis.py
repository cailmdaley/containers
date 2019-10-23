import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from spt3g import core, mapmaker
import nhutils as nhu

def unpack_lr_time(filename, bolos = None):
    f = core.G3File(filename)
    fr = f.next()
    while 'TimestreamWeights' not in fr:
        try:
            fr = f.next()
        except StopIteration as err:
            raise RuntimeError("Didn't find any frames with {}".format('TimestreamWeights'))
    bolokeys = {}
    if bolos is None:
        for bolo in fr['TimestreamWeights'].keys():
            bolokeys[bolo] = ['TimestreamWeights', bolo]
    else:
        for bolo in bolos:
            bolokeys[bolo] = ['TimestreamWeights', bolo]
    for fr in f:
        pass
    return nhu.extract_keys(filename, keys = bolokeys)

def unpack_lr_map(lmap_file, rmap_file, bolos = None):
    lfile = core.G3File(lmap_file)
    rfile = core.G3File(rmap_file)
    ts_weights = {}
    for lfr, rfr in zip(lfile, rfile):
        if lfr.type != core.G3FrameType.Map and rfr.type != core.G3FrameType.Map:
            continue
        if lfr['Id'] != rfr['Id']:
            raise RuntimeError("Bolo ids don't match!  Are you sure these are the same observation?")
        if lfr['Id'] == 'bsmap':
            continue
        if 'Wunpol' in lfr:
            rweight = np.array(rfr['Wunpol'].TT)
            lweight = np.array(lfr['Wunpol'].TT)
        bolo = lfr['Id']
        if bolos is not None and bolo not in bolos:
            continue
        bad_weight = np.where((lweight == 0) | (rweight == 0))
        dmap = np.array(lfr['T'] / lweight - rfr['T'] / rweight)
        dmap[bad_weight] = 0
        good = np.nonzero(dmap) # only use pixels where both maps are nonzero
                                # This is a small change in the number of pixels used,
                                # so it shouldn't make a difference
        if np.product(np.shape(good)) == 0:
            ts_weights[lfr['Id']] = 0
        else:
            ts_weights[lfr['Id']] = 1 / np.var(dmap[good])
    return ts_weights

def unpack_cal(obsid, bolos = None):
    '''obsid of source'''
    calfile = nhu.find_cal(obsid)
    calframe = core.G3File(calfile).next()
    if bolos is None:
        return calframe['CalibratorResponseSN']
    return {bolo: calframe['CalibratorResponseSN'][bolo] for bolo in bolos}


def unpack_elnod(obsid, bolos = None):
    elnodfile = nhu.find_cal(obsid, '/spt/user/production/calibration/elnod')
    elnodframe = core.G3File(elnodfile).next()
    if bolos is None:
        return elnodframe['ElnodSNSlopes']
    return {bolo: elnodframe['ElnodSNSlopes'][bolo] for bolo in bolos}

def unpack_source_signal(filename):
    sourceframe = core.G3File(filename).next()
    return sourceframe['RCW38FluxCalibration']

def unpack_off_source(filename):
    return np.load(filename, encoding = 'latin1').item()

def reduce_time(src_td_weights):
    mean_out = {}
    std_out = {}
    for bolo in src_td_weights.keys():
        mean_out[bolo] = np.mean(1 / np.sqrt(src_td_weights[bolo]))
        std_out[bolo] = np.std(1 / np.sqrt(src_td_weights[bolo]))
    return mean_out, std_out

def reduce_map(src_map_weights):
    noise_out = {}
    for bolo in src_map_weights.keys():
        noise_out[bolo] = 1 / np.sqrt(src_map_weights[bolo])
    return noise_out

def get_all_sn(obsid):
    # obsid = 12078921
    datadir = os.path.join('/spt/user/ndhuang/source_noise/', str(obsid))
    datafiles = ['maps_l.g3', 'maps_r.g3', 'lrnoise_minimal_bins.g3']
    for df in datafiles:
        if not os.path.exists(os.path.join(datadir, df)):
            raise RuntimeError('Missing ' + os.path.join(datadir, df))        
    src_signal = unpack_source_signal('/spt/user/production/calibration/RCW38-pixelraster/{:d}.g3'.format(obsid))
    bolos = src_signal.keys()
    en = unpack_elnod(obsid, bolos)
    cal = unpack_cal(obsid, bolos)
    src_ms = reduce_map(unpack_lr_map(os.path.join(datadir, 'maps_l.g3'),
                                      os.path.join(datadir, 'maps_r.g3'), bolos))
    src_td, junk = reduce_time(unpack_lr_time(os.path.join(datadir, 'lrnoise_minimal_bins.g3'),
                                              bolos))
    src_off = unpack_off_source(os.path.join('/spt/user/ndhuang/source_noise', 
                                             '{:d}MapNoise.npy'.format(obsid)))
    for bolo in src_ms.keys():
        try:
            src_td[bolo] = src_signal[bolo] / src_td[bolo]
            src_ms[bolo] = src_signal[bolo] / src_ms[bolo]
        except KeyError:
            src_td[bolo] = 0.
            src_ms[bolo] = 0.
        if bolo in src_td.keys() and src_td[bolo] == 0:
            print(src_signal[bolo])
            
    all_types = [en, cal, src_ms, src_td]
    # all_types = [en, cal, src_ms]
    all_sn = {}
    for bolo in en.keys():
        try:
            all_sn[bolo] = [v[bolo] for v in all_types]
            all_sn[bolo].append(src_off[bolo.encode()])
        except KeyError:
            pass
    # import IPython
    # IPython.embed()
    return all_sn

def proper_unpack(allsn):
    bolos = [b for b in allsn.keys()]
    cal = np.zeros(len(bolos))
    en = np.zeros(len(bolos))
    src_ms = np.zeros(len(bolos))
    src_td = np.zeros(len(bolos))
    src_off = np.zeros(len(bolos))
    for i, b in enumerate(bolos):
        en[i], cal[i], src_ms[i], src_td[i], src_off = allsn[b]
    return bolos, cal, en, src_ms, src_td, src_off

def iterative_drop(a, thresh = 5, return_inds = False, return_both = False):
    inds_used = np.arange(len(a), dtype = int)
    inds = np.where(np.isfinite(a))
    a = a[inds]
    inds_used = inds_used[inds]
    inds = range(len(a) + 1)
    while len(inds) != len(a) and len(a) > 0:
        if len(a) >= len(inds):
            a = a[inds]
            inds_used = inds_used[inds]
        std = np.std(a)
        mean = np.mean(a)
        sig = np.abs(a - mean) / std
        inds = np.where(sig < thresh)[0]
    if return_both:
        return inds_used, a
    if return_inds:
        return inds_used
    return a

def get_bolo_props(obsid):
    datadir = os.path.join('/spt/user/ndhuang/source_noise/', str(obsid))
    f = core.G3File(os.path.join(datadir, 'maps_l.g3'))
    fr = f.next()
    return fr['BolometerProperties']

def make_plots(obsid, directory = '/home/ndhuang/plots/source_noise'):
    bandspec = {90: 'bx',
                150: 'go',
                220: 'kv'}
    try:
        allsn = get_all_sn(obsid)
    except Exception as err:
        print(err)
        return
    stuff = proper_unpack(allsn)
    bolos = np.array(stuff[0])
    # all = np.array(stuff[1:])
    all = np.array(stuff[3:])
    boloprop = get_bolo_props(obsid)
    labels = ['Map space', 'Time domain', 'Off source']
    lims = [(0, 0) for i in labels]
    # for i in range(len(all)):
    #     thing = all[i]
    #     inds = iterative_drop(thing, thresh = 5, return_inds = True)
    #     if len(inds) == 0:
    #         raise RuntimeError('Killed all the bolos on {}'.format(labels[i]))
    #     all = all[:, inds]
    #     bolos = bolos[inds]
    for i, thing in enumerate(all):
        if np.mean(thing) < 0:
            all[i] = -thing
    # lims = [(20, 200), (0, 700), (0, 50), (0, 50)]
    lims = [(0, 50), (0, 50), (0, 50)]
    fig, axes = pl.subplots(2, 2, sharex = 'col', sharey = 'row')
    actually_live_bolos = sum(np.isfinite(np.sum(all, 0)))
    bands = np.array([boloprop[bolo].band / core.G3Units.GHz for bolo in bolos])
    band_inds = {90.: np.where(bands == 90.), 
                 150.: np.where(bands == 150.),
                 220.: np.where(bands == 220.)}
    # import IPython
    # IPython.embed()
    for i in range(2):
        for j in range(i + 1, 3):
            ax = axes[1 - i, 2 - j]
            for band, bi in band_inds.items():
                ax.plot(all[j][bi], all[i][bi], bandspec[band], label = '{:0.0f} GHz'.format(band))
            pl.legend()
            if i == 0:
                ax.set_xlabel(labels[j])
            if j == 3:
                ax.set_ylabel(labels[i])
            xlim = ax.get_xlim()
            ax.set_xlim(0, xlim[1])
            ylim = ax.get_ylim()
            ax.set_ylim(0, ylim[1])
            # if i == 1 and j == 2:
            #     ax.set_title('# Reasonable bolos: {}'.format(actually_live_bolos))
            # if i == 0 and j == 1:
            #     ax.set_title('# Additional bolos cut: {}'.format(len(bolos) - np.shape(all)[1]))
    # for i, j in [[0, 1], [0, 2], [1, 2]]:
        # fig.delaxes(axes[i, j])
    fig.delaxes(axes[1, 0])
    fig.set_figheight(10)
    fig.set_figwidth(12)
    fig.tight_layout()
    fig.savefig(os.path.join(directory, str(obsid) + '_triangle.png'))
    pl.close(fig)

def make_other_plots(obsid):
    datadir = os.path.join('/spt/user/ndhuang/source_noise/', str(obsid))
    datafiles = ['maps_l.g3', 'maps_r.g3', 'lrnoise_minimal_bins.g3']
    for df in datafiles:
        if not os.path.exists(os.path.join(datadir, df)):
            raise RuntimeError('Missing ' + os.path.join(datadir, df))        
    src_signal = unpack_source_signal('/spt/user/production/calibration/RCW38-pixelraster/{:d}.g3'.format(obsid))
    bolos = src_signal.keys()
    src_ms = reduce_map(unpack_lr_map(os.path.join(datadir, 'maps_l.g3'),
                                      os.path.join(datadir, 'maps_r.g3'), bolos))
    src_td, junk = reduce_time(unpack_lr_time(os.path.join(datadir, 'lrnoise_minimal_bins.g3'),
                                              bolos))
    src_off = unpack_off_source(os.path.join('/spt/user/ndhuang/source_noise', 
                                             '{:d}MapNoise.npy'.format(obsid)))
    bolos = src_ms.keys()
    tdom = np.zeros(len(bolos))
    mapsp = np.zeros(len(bolos))
    offsrc = np.zeros(len(bolos))
    for i, bolo in enumerate(bolos):
        try:
            tdom[i] = src_signal[bolo] / src_td[bolo]
            mapsp[i] = src_signal[bolo] / src_ms[bolo]
            offsrc[i] = src_signal[bolo] / src_off[bolo]
        except KeyError:
            tdom[i] = 0.
            mapsp[i] = 0.
            offsrc[i] = 0.
    boloprop = get_bolo_props(obsid)
    bands = np.array([boloprop[bolo].band / core.G3Units.GHz for bolo in bolos])
    band_inds = {90.: np.where(bands == 90.), 
                 150.: np.where(bands == 150.),
                 220.: np.where(bands == 220.)}
    bandspec = {90: 'bx',
                150: 'go',
                220: 'kv'}
    f = pl.figure()
    f.set_size_inches(12.8, 9.6)
    pl.subplot(221)
    for band, bi in band_inds.items():
        pl.plot(mapsp[bi], tdom[bi], bandspec[band], label = '{:0.0f} GHz'.format(band))
    pl.ylabel('Time domain')

    pl.subplot(223)
    for band, bi in band_inds.items():
        pl.plot(mapsp[bi], offsrc[bi], bandspec[band], label = '{:0.0f} GHz'.format(band))
    pl.ylabel('Off source')
    pl.xlabel('Map space')

    pl.subplot(224)
    for band, bi in band_inds.items():
        pl.plot(tdom[bi], offsrc[bi], bandspec[band], label = '{:0.0f} GHz'.format(band))
    pl.xlabel('Time domain')
    
    
