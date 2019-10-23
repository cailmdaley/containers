import os, glob
import re
try:
    import cPickle as pickle
except ImportError:
    import pickle
import numpy as np
import matplotlib.cm
from matplotlib import pyplot as pl
from spt3g import core, calibration
import nhutils as nhu

def net(fname, plotdir = None):
    obsid = None
    for junks in os.path.basename(fname).split('_'):
        for junk in junks.split('.'):
            try:
                obsid = int(junk)
            except ValueError:
                continue
    if obsid is None:
        print (os.path.basename(fname).split('_'))
        raise RuntimeError("Can't detect obsid")
    with open(fname, 'rb') as f:
        nets_by_band = pickle.load(f)
    pl.figure()
    bins = np.linspace(0, 4000, num = 101)
    bin_centers = bins[1:] - np.diff(bins[:2])
    modes = {}
    for band in sorted(nets_by_band.keys()):
        nets, bolos = nets_by_band[band]
        inds = np.where(np.isfinite(nets))
        try:
            n, junk, morejunk = pl.hist(nets[inds], bins = bins, alpha = .6,
                                        label = '{:d} GHz'.format(band),
                                        histtype = 'stepfilled')
        except UnboundLocalError:
            print('{:d} GHz band has no detectors below {:.00f} uk-st(s)'.format(band, max(bins)))
            continue
        modes[band] = bin_centers[np.argmax(n)]
    pl.legend()
    pl.title('NETs by band (obsid: {:d})'.format(obsid))
    pl.xlabel('NET ($\mu K \sqrt{s}$)')
    pl.ylabel('Number of bolometers')
    if plotdir is not None:
        pl.savefig(os.path.join(plotdir, str(obsid) + '.png'))
    return modes

def focal_plane(fname):
    fr = core.G3File(fname).next()
    if 'BolometerProperties' not in fr:
        obsid = int(nhu.guess_obsid(fname))
        bprop = core.G3File(nhu.find_cal(obsid, '/poleanalysis/sptdaq/calresult/calibration/boloproperties/')).next()['BolometerProperties']
        bands = nhu.split_bolo_list_band(fr['PointingOffsetX'].keys(), bprop, False)
    else:
        bprop = fr['BolometerProperties']
        bands = calibration.template_groups.get_template_groups(bprop, 
                                                                include_keys = True)
    
    for band, bolos in bands.items():
        pl.figure()
        if 'BolometerProperties' not in fr:
            pl.scatter([fr['PointingOffsetX'][bolo] for bolo in bolos],
                       [fr['PointingOffsetY'][bolo] for bolo in bolos],
                       alpha = .5, edgecolors = 'none')
        else:
            pl.scatter([bprop[b].x_offset for b in bolos],
                       [bprop[b].y_offset for b in bolos],
                       alpha = .5, edgecolors = 'none')
        pl.title(band)
        pl.xlim(-.02, .02)
        pl.ylim(-.017, .017)
        ax = pl.gca()
        ax.set_aspect('equal')
        
def singlebolomaps(fname, plotdir = None):
    for fr in core.G3File(fname):
        if fr.type != core.G3FrameType.Map or fr['Id'] == 'bsmap':
            continue
        if 'Wunpol' in fr:
            weight = fr['Wunpol'].TT
        f = pl.figure()
        pl.imshow(fr['T'] / weight, interpolation = 'None')
        ax = f.axes[0]
        ax.axis('off')
        pl.colorbar()
        pl.title(fr['Id'])
        if plotdir is not None:
            pl.savefig(os.path.join(plotdir, fr['Id'] + '.png'))
            pl.close()

def net_over_time(net_dir):
    predicted = {90: 731,
                 150: 442,
                 220: 1173}
    nets_by_oid = {}
    for fname in glob.glob(os.path.join(net_dir, '*.pkl')):
        oid = nhu.guess_obsid(fname)
        with open(fname, 'rb') as f:
            nets_by_oid[oid] = pickle.load(f)
    obsids = sorted(nets_by_oid.keys())
    net_modes = {}
    net_medians = {}
    for b in [90, 150, 220]:
        net_modes[b] = []
        net_medians[b] = []
        for oid in obsids:
            try:
                nets, bolos = nets_by_oid[oid][b]
            except KeyError:
                net_modes[b].append(np.nan)
                net_medians[b].append(np.nan)
                continue
            nets = nets[np.isfinite(nets)]
            range = np.round(np.max(nets), -1) - np.round(np.min(nets), -1)
            bins = np.linspace(np.round(np.min(nets), -1) - 5, 
                               np.round(np.max(nets), -1) + 5,
                               num = range // 10 + 1)
            bin_centers = bins[:-1] + 5
            hist, derp = np.histogram(nets, bins = bins)
            net_modes[b].append(bin_centers[np.argmax(hist)])
            net_medians[b].append(np.median(nets))
    times = [nhu.obsid_to_datetime(oid) for oid in obsids]
    fig = pl.figure(figsize = (10, 10))
    for b in [90, 150, 220]:
        pl.subplot(211)
        pl.title('Mode')
        l = pl.plot(times, net_modes[b], '.-', label = '{} GHz'.format(b),
                    linewidth = 3)[0]
        xl = pl.xlim()
        pl.hlines(predicted[b], *xl, linestyle = '--', color = l.get_color())
        pl.xlim(*xl)
        pl.subplot(212)
        pl.title('Median')
        pl.plot(times, net_medians[b], '.-', label = '{} GHz'.format(b),
                linewidth = 3)
    pl.subplot(211)
    pl.title('Mode')
    pl.ylabel('NET ($\mu K \sqrt{Hz}$)')
    pl.legend()
    pl.subplot(212)
    pl.title('Median')
    pl.suptitle('NET over time')
    pl.ylabel('NET ($\mu K \sqrt{Hz}$)')
    pl.legend()
    fig.autofmt_xdate()
    pl.xlim()
    pl.show()
    return net_medians, net_modes

def field_maps(fname, pol = 'T', plotdir = None):
    maps = list(core.G3File(fname))[1:]
    pl.figure(figsize = (6, 10))
    obsid = nhu.guess_obsid(fname)
    t = nhu.obsid_to_time(obsid)
    pl.suptitle(t.GetFileFormatString())
    for i, m in enumerate(maps):
        pl.subplot(3, 1, i + 1)
        weight = getattr(m['Wpol'], pol * 2)
        pl.imshow(m[pol] / weight / core.G3Units.K, vmin = -.01, vmax = .01,
                  interpolation = 'None', cmap = matplotlib.cm.viridis)
        pl.title(m['Id'])
        pl.axis('off')
    # pl.tight_layout()
    if plotdir is not None:
        pl.savefig(os.path.join(plotdir, str(obsid) + '.png'))
