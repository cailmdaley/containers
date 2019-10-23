import os
import numpy as np
from matplotlib import pyplot as pl
from scipy import ndimage, optimize
from spt3g import core
import focus

def profile_gauss((A, sig_x, sig_y, mu_x, mu_y, offset), xy):
    sig = (sig_x, sig_y)
    mu = (mu_x, mu_y)
    if abs(sig_x) > 30 or abs(sig_y) > 30:
        return np.ones(len(xy[0])) * (abs(sig_x) + abs(sig_y))
    return offset + A * np.exp(-sum([(mu[i] - xy[i])**2 / (2 * sig[i]**2)
                                     for i in range(2)]))

def fit_focus(m, res = .25, verbose = False):
    xlen, ylen = np.shape(m)
    mu = np.array(np.shape(m)) / 2.
    # sig = np.array([1.5, 1.5]) / res
    sig = 1.5 / res
    xy = np.meshgrid(*np.shape(m))
    pars = [.5, sig, sig, mu[0], mu[1], 0]
    # p, cov = optimize.curve_fit(lambda x, *pars: profile_gauss(pars, x), xy, 
    #                        m[xy[0], xy[1]], p0 = pars)
    out = opt.leastsq(lambda pars, xy: profile_gauss(pars, xy) - m[xy[0], xy[1]],
                      pars, args = (xy,), full_output = True)
    p = out[0]
    if verbose:
        fwhm = p[1:3] * np.sqrt(8 * np.log(2)) * res / core.G3Units.arcmin
        print('FWHM (arcmin) x  \ty\n            {:.02f}\t{:.02f}'.format(*abs(fwhm)))
        # print(out[1])
        # print(out[2])
        # print(out[3])
        # fwhm = p[1] * np.sqrt(8 * np.log(2)) * res / core.G3Units.arcmin
        # print('FWHM: {:.02f} arcmin'.format(fwhm))
    return p
    
# def make_xy_flat(x, y):
#     try:
#         grid = np.zeros(len(x), len(y), 2, dtype = int)
#     except TypeError:
#         grid = np.zeros((x, y, 2), dtype = int)
#         x = range(x)
#         y = range(y)
#     for i, _x in enumerate(x):
#         for j, _y in enumerate(y):
#             grid[i, j] = (_x, _y)
#     return grid.reshape(len(x) * len(y), 2).T

def get_bench_pos(obsid, bolodata = '/spt_data/bolodata/fullrate/', 
                  source = '0537-441-pixelraster'):
    f = core.G3File(os.path.join(bolodata, source, str(obsid), '0000.g3'))
    for fr in f:
        if 'BenchCommandedPosition' in fr.keys():
            return fr['BenchCommandedPosition']

def reconstruct_map(pars, shape):
    xy = make_xy_flat(*shape)
    return profile_gauss(pars, xy).reshape(*shape)

def fit_one_ob(obsid, caldir = '/poleanalysis/sptdaq/calresult/calibration/',
               source = '0537-441-pixelraster'):
    maps = core.G3File(os.path.join(caldir, source, '{}.g3'.format(obsid))) 
    results = {}
    smoothed_maps = {}
    fits = {}
    for mfr in maps:
        band = int(mfr['Id'].split('-')[-1].strip('GHz')) * core.G3Units.GHz
        m = mfr['T']
        # get the central square degree of the map
        xcenter = m.shape[0] / 2
        ycenter = m.shape[1] / 2
        npix = int(.25 / (m.res / core.G3Units.deg))
        m = np.array(mfr['T'] / mfr['Wunpol'].TT)[xcenter - npix:xcenter + npix,
                                                  ycenter - npix:ycenter + npix]
        # smooth, hopefully not too aggressively
        m = ndimage.gaussian_filter(m, 1)
        smoothed_maps[band] = m
        if np.sum(np.isfinite(m)) == 0:
            print('{} is all NaNs.'.format(band))
            continue
        pars = fit_focus(m, res = mfr['T'].res, verbose = False)
        results[band] = {'fwhm': abs(pars[1:3]) * mfr['T'].res 
                         / core.G3Units.arcmin * np.sqrt(8 * np.log(2)),
                         'loc': pars[3:5], 'full': pars}
        fits[band] = reconstruct_map(pars, np.shape(m))
    return results, smoothed_maps, fits

def oplot_fits(band, maps, fits):
    pl.figure(figsize = (10, 10))
    pl.imshow(maps[band])
    pl.imshow(fits[band], alpha = .4)

def fit_many(obsids, caldir = '/poleanalysis/sptdaq/calresult/calibration/',
             source = '0537-441-pixelraster'):
    out = {}
    for obsid in obsids:
        results, smooth_maps, fit_maps = fit_one_ob(obsid, caldir, source)
        bench = get_bench_pos(obsid, source = source)
        optical = focus.bench2optical(*focus.g3bench2xyz(bench))
        out[obsid] = {'focus': optical, 'fit_results': results,
                      'smooth_maps': smooth_maps, 'fit_maps': fit_maps} 
    return out

def unpack_results(res, band):
    obsids = sorted(res.keys())
    fwhms = np.array([np.sqrt(sum(res[obsid]['fit_results'][band]['fwhm']**2)) 
                      for obsid in obsids])
    benches = np.array([np.array(res[obsid]['focus']).flatten()
                      for obsid in obsids])
    return fwhms, benches
