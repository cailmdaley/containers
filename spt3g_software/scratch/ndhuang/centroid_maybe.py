import numpy as np
from scipy import optimize as opt
from scipy import ndimage
from spt3g import core, mapmaker

def profile_gauss(xy, A, sigma, mu_x, mu_y, offset):
    mu = (mu_x, mu_y)
    # import IPython
    # IPython.embed()
    return offset + A * np.exp( -np.sum([(xy[i] - mu[i])**2 / 
                                         (2 * sigma**2)
                                         for i in (0, 1)], axis = 0))
def minimer(pars, xy, m):
    resid = profile_gauss(xy, *pars) - m
    chi2 = np.sum(resid**2)
    # Apply some penalties, 'cause I'm a dyed-in-the wool Bayesian, or something
    penalty = 1
    return chi2 * penalty
    # A, sigma, mu_x, mu_y, offset = pars
    # xpix = xy[0][0, -1]
    # ypix = xy[1][-1, -1]
    # if mu_x < 0 or mu_x > xpix:

def fit_map(m, reso_arcmin):
    xpix, ypix = np.shape(m)
    xy = np.meshgrid(*[np.arange(pix) for pix in (ypix, xpix)])
    mfilt = ndimage.gaussian_filter(m, 5)
    xcen, ycen = np.unravel_index(np.nanargmax(mfilt), np.shape(m))
    pars = [np.max(m), .6 / reso_arcmin, 
            xcen, ycen,
            np.mean(m[:10, :10])]
    # p, cov = opt.curve_fit(lambda x, *pars: profile_gauss(x, *pars).ravel(), 
    #                        xy, m.ravel(), p0 = pars)

    # minimer = lambda pars, xy: profile_gauss(xy, *pars).ravel() - m.ravel()
    # out = opt.leastsq(minimer, pars, args = (xy,),
    #                   full_output = True)
    # return out
    
    bounds = ((0, None), # A
              (.5, 5 / reso_arcmin), # sigma
              (0, xpix), # mu_x
              (0, ypix), # mu_y
              (None, None)) # offset
    out = opt.minimize(minimer, pars, args = (xy, m),
                       bounds = bounds)
    p = out.x, out.fun
    # do some translating to sensible units
    # p[1:3] *= reso_arcmin
    return p

def get_ptsrc_centers(filename):
    f = core.G3File(filename)
    results = {}
    for fr in f:
        source = fr['Id'].split('_')[0]
        band = int(fr['Id'].split('_')[-1].split('GHz')[0])
        direction = fr['Id'].split('GHz')[-1]
        if band not in results.keys():
            results[band] = {}
        res = fr['T'].res / core.G3Units.arcmin
        m = np.array(fr['T'] / fr['Wunpol'].TT)
        xpix, ypix = np.shape(m)
        halfpix = 5 / res
        m = m[xpix // 2 - halfpix:xpix // 2 + halfpix,
              ypix // 2 - halfpix:ypix // 2 + halfpix]
        try:
            p = fit_map(m, res)
        except RuntimeError:
            continue
        i = 0 if direction == 'Left' else 1
        if source not in results[band].keys():
            results[band][source] = np.zeros((2, len(p)))
        results[band][source][i] = p
    # for band in results.keys():
    #     for source, fits in results[band].items():
    #         radiff = np.diff(fits[:, 2]) * res
    #         results[band][source] = (radiff, 
    return results

def collect_src_fits(filename):
    f = core.G3File(filename)
    results = {}
    for fr in f:
        res = fr['T'].res / core.G3Units.arcmin
        m = np.array(fr['T'] / fr['Wunpol'].TT)
        xpix, ypix = np.shape(m)
        halfpix = 5 // res
        s = slice(int(xpix // 2 - halfpix), int(xpix // 2 + halfpix))
        m = m[s, s] # assuming square map
        # results[fr['Id']] = (fit_map(m, res), m)
        results[fr['Id']] = fit_map(m, res)
    return results

def analyze_files(files):
    results = {90: [], 150: [], 220: []}
    for f in files:
        res = collect_src_fits(f)
        splitter = {}
        for src in res:
            # Do cuts, collect Ra offsets
            resid = res[src][1]
            sigma = res[src][0][1]
            if sigma < 1 or sigma > 2 or resid > 3000:
                continue
            source = src.split('_')[0]
            band = int(src.split('_')[-1].split('GHz')[0])
            direction = src.split('GHz')[-1]
            splitter[source] = {band: {direction: res[src][0][2] * .4}}
        for source in splitter:
            # collect delta Ra
            for band in splitter[source]:
                if len(splitter[source][band].keys()) != 2:
                    # left or right was cut, so skip
                    continue
                results[band].append(splitter[source][band]['Left'] - 
                                     splitter[source][band]['Right'])
    return results
            


# Cuts:
#    residual < 3000
#    1 < sigma < 2
