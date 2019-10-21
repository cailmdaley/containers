"""
Lots of stuff that I (TC) use a lot
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import pylab
from scipy import ndimage
from spt3g import core, calibration
import glob
import pickle
#from sptpol_software import util
#from sptpol_software.util import math
#from sptpol_software.util import tools

#################################
def shift(array, shift):
    """
    A wrapper for np.roll. You may specify a shift for each axis, and this
    function will perform each shift at once. This method of specifying the shift is
    sometimes more convenient.

    INPUT
       array: (ndarray) The array which you wish to shift.

       shift: (iterable) The shift of each axis of the input array. "shift" must have
          a number of elements equal to the number of dimensions of "array"!

    OUTPUT
       A shifted version of the input ndarray.
    """
    # Check the input array sizes.
    if len(array.shape) != len(shift):
        raise ValueError("You must specify a shift for each axis of the array!")
       
    output = array
    for axis, this_shift in enumerate(shift):
        if this_shift:
            # Only make this shift if it's non-zero.
            output = np.roll(output, this_shift, axis=axis)
        
    return output

#################################
def makeFFTGrid(shape, resolution=1.):
    """
    Creates a grid in which each element is proportional to the
    distance from the point to the origin (with wrapping).
    The distance between neighboring grid squares is 1/(n*resolution).
    
    >>> utils.makeFFTGrid(10, 2.)
    array([ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25, -0.2 , -0.15, -0.1 , -0.05])
    
    >>> utils.makeFFTGrid([5,5], 1.)
    array([[ 0.        ,  0.2       ,  0.4       ,  0.4       ,  0.2       ],
           [ 0.2       ,  0.28284271,  0.4472136 ,  0.4472136 ,  0.28284271],
           [ 0.4       ,  0.4472136 ,  0.56568542,  0.56568542,  0.4472136 ],
           [ 0.4       ,  0.4472136 ,  0.56568542,  0.56568542,  0.4472136 ],
           [ 0.2       ,  0.28284271,  0.4472136 ,  0.4472136 ,  0.28284271]])
    """
    # Make sure the input "shape" is an iterable.
    try:
        len(shape)
    except TypeError:
        shape = [shape]

    # Create a set of coordinates. We can then find the "distance" of
    # each point from the origin.
    slices = [slice( -np.floor((dim-1)/2.), np.ceil((dim-1)/2.)+1 ) for dim in shape]
    coordinates = np.mgrid[slices]

    # Get the correct normalization - divide each dimension by
    # its size and by the overall resolution.
    for index, dim_size in enumerate(shape):
        coordinates[index] /= (dim_size*resolution)

    # If we want a 1D array, then we already have what we want.
    # Don't apply the square and square-root, so that we can
    # keep the negative signs.
    if len(shape)!=1:
        grid = np.sqrt( np.sum( coordinates**2, axis=0 ) )
    else:
        grid = coordinates[0]
    
    for index, dim_size in enumerate(shape):
        grid = np.roll(grid, int(np.ceil((dim_size-1)/2.)+1), axis=index)

    return grid

#################################
def makeEllGrid(shape, resolution=1.):
    """
    OUTPUT 
       2.*np.pi*makeFFTGrid(shape=shape, resolution=resolution)
    """
    return 2.*np.pi*makeFFTGrid(shape=shape, resolution=resolution)


def makeEllGrids(shape, resolution_radians):
    '''
    This returns the 2d ELL grids (ELL, ELL_X, ELL_Y) given a map 
    shape and resolution.  It's like makeEllGrid, but this also 
    returns ELL_X and ELL_Y.

    INPUTS
      shape - the shape of the map you want the ell grids for.
      resolution_radians - the map resolution in radians.
    OUTPUTS
      A tuple of three 2d numpy arrays, (ELL, ELL_X, ELL_Y).

    created 12 Aug 2014, RK
    '''
    ell_x, ell_y = np.meshgrid(2*np.pi*np.fft.fftfreq(shape[1], resolution_radians),
                               2*np.pi*np.fft.fftfreq(shape[0], resolution_radians))
    ell = np.sqrt(ell_x**2 + ell_y**2)
    return ell, ell_x, ell_y

#################################
def smooth_flatsky(grid_in, reso_arcmin, fwhm_arcmin):
    """
    2d smooth with gaussian kernel
    """

    sig_smooth_arcmin = fwhm_arcmin/np.sqrt(8.*np.log(2.))
    sig_smooth_pix = sig_smooth_arcmin/reso_arcmin

    result = scipy.ndimage.gaussian_filter(grid_in,sig_smooth_pix)

    return result


#################################
def gaussbeam2d(ngrid, reso_arcmin, fwhm_arcmin, realspace=False):
    """
    Returns 2d Fourier-space Gaussian beam. If realspace is set to
    True, returns real-space beam instead, centered at
    [ngrid/2,ngrid/2].
    """

    ellg = makeEllGrid([ngrid,ngrid],reso_arcmin/60.*np.pi/180.)
    sig_smooth = fwhm_arcmin/np.sqrt(8.*np.log(2.))/(60.*180.)*np.pi
    bl2d = np.exp(-ellg*(ellg+1.)*sig_smooth**2/2.)

    if realspace:
        bl2d = np.real(np.fft.fft2(bl2d))
        bl2d /= np.max(bl2d)
        bl2d = shift(bl2d,[ngrid/2,ngrid/2])

    return bl2d


#################################
def bin_into_pixels(data, pointing, dweight=1, 
                    use_histogram=True, return_weighted_map=False, return_weights=True):
    """
    Make simple map, using simplest or next-to-simplest algorithm.
    """

    which_alg = 2
    if use_histogram: which_alg = 1

    indp = pointing.astype(int)
    npixels = np.max(indp) + 1
    map = np.zeros(npixels)
    weights = np.zeros(npixels)

# check if dweight is a vector the same length as data
    if hasattr(dweight,"__len__"):
        if len(dweight) == len(data):
            dw2use = dweight
        else:
            dw2use = np.zeros(len(data)) + dweight[0]
    else:
        dw2use = np.zeros(len(data)) + dweight

    if which_alg == 1:
        sp = np.argsort(indp)
        sindp = indp[sp]
        sdata = data[sp]
        [upix,invind] = np.unique(sindp,return_inverse=True)
        histind = (np.histogram(invind,bins=np.max(invind)+1,range=(0,max(invind)+1)))[0]
        ntot = 0
        for i in np.arange(len(histind)):
            h = histind[i]
            if h > 0:
                theseinds = np.arange(h) + ntot
                thispix = upix[i]
                map[thispix] = np.sum(sdata[theseinds]*dw2use[theseinds])
                weights[thispix] = np.sum(dw2use[theseinds])
                ntot += h

    if which_alg == 2:
        for i in np.arange(len(data)):
            map[indp[i]] += data[i]*dw2use[i]
            weights[indp[i]] += dw2use[i]
    
    if return_weighted_map == False:
        whn0 = (np.where(weights > 0))[0]
        map[whn0] /= weights[whn0]

    if return_weights:
        return [map,weights]
    else:
        return map


#################################
def put_point_sources_on_grid(fluxvec, dnds, ngrid, reso_arcmin,
                              mins=1e-6, maxs=100., 
                              bands=150., nu0=150., spectral_index=0., 
                              sigspec=0., nomax=False, flux_cut=1e6, 
                              fcband=-100., return_map=True):
    """

    """

# bookkeeping
    reso_rad = reso_arcmin/60.*np.pi/180.
    npix = ngrid**2
    nsr = np.float(npix)*reso_rad**2
    ndeg = nsr*(180./np.pi)**2
    nstot = 0

# flux bins
    s = fluxvec
    maxs2use = np.min([maxs,np.max(s)])
    mins2use = np.max([mins,np.min(s)])
    dlogs = 0.1
    nbins_f = np.log10(maxs2use/mins2use)/dlogs + 1.
    nbins = np.int(np.floor(nbins_f))
    lmins = np.log10(mins2use)
    lmaxs = np.log10(maxs2use)
    logs = np.log10(s)

# interpolate onto grid to sample
    logs2int = lmins + np.array(range(nbins+1))*dlogs
    wh0 = (np.where(dnds <= 0.))[0]
    dnds2int = dnds
    dnds2int[wh0] = 1e-24
    ldnds_int = np.interp(logs2int,logs,np.log10(dnds2int))

# check expected # of output sources
    wh2int = np.where(np.logical_and(s > mins2use,s < maxs2use))
    nsource_expected = np.trapz(dnds[wh2int],s[wh2int])*ndeg
    if (nsource_expected > 1e8 and nomax == False):
        print('Do you really want >10,000,000 sources \n')
        print('actually, ' + repr(nsource_expected) + ', give or take a few)?\n')
        print('If so, set nomax = True.')
        flux = -1.
        return

#    print('nsource_exp: ' + repr(nsource_expected) + '.\n')

# initialize output arrays
    nsource_safe = nsource_expected*2
    if nsource_safe < 1000: nsource_safe=1000
    flux = np.zeros(nsource_safe)
    pixno_f = np.zeros(nsource_safe)

# step through flux bins, creating poisson realization of sources in
# each
    for i in np.arange(nbins):
        lsmin = logs2int[i]
        lsmax = lsmin + dlogs
        ldnds_mean = (ldnds_int[i+1] + ldnds_int[i])/2
        nsource_mean = 10**(ldnds_mean+(lsmax+lsmin)/2)*dlogs*np.log(10)*ndeg
        nsource = np.random.poisson(nsource_mean)
        if nsource > 0:
            lflux = np.random.uniform(lsmin,lsmax,nsource)
            flux[nstot:nstot+nsource] = 10**lflux
            pixno_f[nstot:nstot+nsource] = np.floor(np.random.uniform(0,npix,nsource))
            nstot += nsource

# get rid of empty array elements
    flux = flux[0:nstot]
    pixno_f = pixno_f[0:nstot]
    pixno = pixno_f.astype(int)

    if return_map:
# put sources on map
        psmap = np.zeros(npix)
        if nstot > 0:
            [psmtemp,nsperpix] = bin_into_pixels(flux,pixno,return_weighted_map=True)
            psmap[0:len(psmtemp)] = psmtemp
        return [flux,pixno,psmap]
    else:
        return [flux,pixno]
    

#################################
def counts_bpl(fluxvec, params):
    """
    Returns broken-power-law model dN/dS for supplied flux vector (in Jy) and
    model params (amplitude at 1 Jy, bright-end power law, break flux,
    faint-end power law).
    """
    
    flux2use = np.copy(fluxvec)
    flux2use[(np.where(fluxvec <= 0.))[0]] = 1e-24
    ampl_1jy = params[0]
    alpha_b = params[1]
    f_break = params[2]
    fbg1 = f_break > 1.
    alpha_f = params[3]
    whb = (np.where(fluxvec > f_break))[0]
    whf = (np.where(fluxvec < f_break))[0]
    ls = np.log(fluxvec)
    lsb = np.log(f_break)
    if fbg1:
        ldnds = alpha_f*ls
        ldnds[whb] = alpha_b*ls[whb] + alpha_f*lsb - alpha_b*lsb
    else:
        ldnds = alpha_b*ls
        ldnds[whf] = alpha_f*ls[whf] + alpha_b*lsb - alpha_f*lsb
    dnds = np.exp(ldnds)*ampl_1jy

    return dnds


#################################
def ud_grade_grid(grid_in, fac_rebin, do_up=False, do_sum=False):
    """
    Re-bin a 2-dimensional grid, either to a larger grid by resampling
    or to a smaller grid by averaging. If do_up is set to True, the
    output grid is larger than the input by a factor of fac_rebin on a
    side; otherwise it is smaller by a factor of fac_rebin on a side.
    """
    
    if np.float(fac_rebin) <= 1.:
        if np.float(fac_rebin) < 1.:
            print("UD_GRADE_GRID: TO GET AN OUTPUT GRID THAT IS SMALLER THAN THE INPUT,")
            print("SET FAC_REBIN TO THE FACTOR (> 1) BY WHICH YOU WANT TO REDUCE THE SIZE.\n")
        grid_out = np.copy(grid_in)
        return grid_out

    shape_orig = np.shape(grid_in)

    fac2use = np.int(np.floor(fac_rebin))
    if do_up:
        grid_out = np.zeros((shape_orig[0]*fac2use,shape_orig[1]*fac2use))
        for i in np.arange(shape_orig[0]):
            for j in np.arange(shape_orig[1]):
                grid_out[fac2use*i:fac2use*(i+1),fac2use*j:fac2use*(j+1)] = grid_in[i,j]
    else:
        nx_new = np.int(np.floor(np.float(shape_orig[0])/np.float(fac2use)))
        ny_new = np.int(np.floor(np.float(shape_orig[0])/np.float(fac2use)))
        grid_out = np.zeros((nx_new,ny_new))
        for i in np.arange(nx_new):
            for j in np.arange(ny_new):
                grid_out[i,j] = np.mean(grid_in[fac2use*i:fac2use*(i+1),fac2use*j:fac2use*(j+1)])
        

    return grid_out


#################################
def gaussfit_hist(values, nbins, minval, maxval, 
                  do_plot=True, do_oplot=False, do_log=False):
    """
    histogram some values, fit the histogram to a gaussian, and return
    the best-fit parameters.
    """
    
    binsize = (maxval-minval)/(np.float(nbins-1))
    ht = np.histogram(values*1.0,bins=nbins,range=(minval,maxval))
    dx = ht[1][1] - ht[1][0]
    x_out = ht[1][0:nbins] + dx/2.
    hist_out = ht[0]

# throw out some outliers before fittin
    (mutemp,sigmatemp) = scipy.stats.norm.fit(values)
    whgood = np.where(np.abs(values-mutemp)/sigmatemp < 5.)
    (mu,sigma) = scipy.stats.norm.fit(values[whgood])
    yfit_unnorm = np.exp(-(x_out-mu)**2/2./sigma**2)
    ampl = np.sum(yfit_unnorm*hist_out)/np.sum(yfit_unnorm**2)
    yfit = ampl*yfit_unnorm
    hfit_out = yfit

    params = [ampl,mu,sigma]

    if do_plot:
        plt.hist(values, bins=nbins, range=(minval,maxval), color='g')
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
#        p = scipy.stats.norm.pdf(x, mu, sigma)
        p = ampl*np.exp(-(x-mu)**2/2./sigma**2)
        plt.plot(x, p, 'k', linewidth=2)
        title = "Fit results: mu = %.4g,  std = %.4g" % (mu, sigma)
        plt.title(title)

#        if do_log:
#            loglog(x_out,hist_out,'k',linewidth=2)

    dict_out = {'params':params,'x_out':x_out,'hist_out':hist_out,'hfit_out':hfit_out}


    return dict_out


#################################
def cov2corr(cov):

    """
    convert covariance matrix to correlation matrix
    """
    
    corr = np.copy(cov)
    corr[:,:] = 0.
    nrow = np.size(cov[:,0])

    for i in np.arange(nrow):
        idiag = cov[i,i]
        if np.abs(idiag) > 0:
            for j in np.arange(nrow):
                jdiag = cov[j,j]
                if np.abs(idiag) > 0:
                    corr[i,j] = cov[i,j]/np.sqrt(idiag*jdiag)

    return corr


#################################
def griddist(ysize, xsize=None):

    """
    grid of distances, with 0 at [0,0] and wrapped. equivalent to IDL
    dist.pro.
    """

    if xsize is None:
        xsize = ysize

    grids = np.mgrid[0:ysize,0:xsize]
    gdist = np.roll(np.roll(np.sqrt((grids[0,:,:]-ysize/2.)**2 + (grids[1,:,:]-xsize/2.)**2),ysize/2,axis=0),xsize/2,axis=1)

    return gdist

#################################
def int_src(map_in, r_cutoff, xcenter=None, ycenter=None, resolution=1.):

    """
    integrate around a point in a map. answer returned in units of map
    intensity units times map area units (so if you give it a map in K
    and resolution in radians, the answer is in K-sr).
    """

    amap = np.asarray(map_in)
    mapshape = np.shape(amap)
    
    if xcenter is None or ycenter is None:
        if xcenter is None and ycenter is None:
            amap_sm = ndimage.gaussian_filter(amap,np.int(2./resolution))
            ycenter, xcenter = np.unravel_index(np.argmax(np.abs(amap_sm)),[mapshape[0],mapshape[1]])
        else:
            print('You have supplied an x coordinate to integrate around but not a y coordinate (or vice-versa).')
            return 0
        
    dist2d = shift(griddist(mapshape[0],mapshape[1]),[np.round(ycenter),np.round(xcenter)])*resolution
    whint = np.where(dist2d <= r_cutoff)
    flux = np.sum(amap[whint])*resolution**2

    return flux

#################################
def quick_net_from_rcw38(rcw38_map, data1d, band=150, xcenter=None, ycenter=None, r_cutoff=5., reso_arcmin=0.5, rate=152.58789, frange=[10.,20.], data_is_wn=False, invert_map=False):

    """
    calculate NET from single-bolo map of rcw38 and some noise data
    (timestream). both should be in the same units (assuming ADC
    counts for now, but it doesn't matter as long as they're the
    same). result is NET in K_RJ-sqrt(s).
    """

    bands = [90, 150, 220]
#    rcw38_integrated_k_sr = np.array([6.4, 2.0, 2.3])*1e-7 # from https://anal.spt/trac/wiki/QuickNets06Feb08
    rcw38_integrated_k_sr = np.array([4.7, 2.0, 1.3])*1e-7 # update with 5'-radius integration of Planck-recalibrated SPT-SZ maps
    rcw38_band = (rcw38_integrated_k_sr[np.where(np.asarray(bands) == np.int(band))])[0]
    rcw38_bolo = int_src(rcw38_map, r_cutoff, xcenter=xcenter, ycenter=ycenter, resolution=reso_arcmin)*core.G3Units.arcmin**2 # in input units-sr
    if invert_map:
        rcw38_bolo *= -1.
    cts_per_k = rcw38_bolo / rcw38_band # input-units / K
    if data_is_wn:
        whitelevel = data1d
    else:
        npts_psd = np.min([1024,((len(data1d)-10)/4)*2])
        psdict = quick_pspec(data1d,rate=rate,npts_psd=npts_psd) # in input units / sqrt(Hz)
        whint = np.where(np.logical_and(psdict['freqs'] > frange[0],psdict['freqs'] < frange[1]))
        whitelevel = np.sqrt(np.mean(psdict['psd'][whint]**2))
    net = whitelevel / cts_per_k / np.sqrt(2.) # in K-sqrt(s)

    return net, cts_per_k, whitelevel, rcw38_bolo

#################################                                                          
def quick_net_from_planet(planet_map, planet_name, mjd, data1d, band=150, xcenter=None, ycenter=None, r_cutoff=5., reso_arcmin=0.5, rate=152.58789, frange=[10.,20.], data_is_wn=False, invert_map=False):

    from spt3g.calibration import planets

    """                                                                                     
    calculate NET from single-bolo map of a planet and some noise data                      
    (timestream). both should be in the same units (assuming ADC                            
    counts for now, but it doesn't matter as long as they're the                            
    same). result is NET in K_RJ-sqrt(s).                                                   
    """

    integrated_k_sr_band = planets.temperature_solid_angle_product(planet_name,mjd,band,TCMB=False)
    integrated_k_sr_bolo = int_src(planet_map, r_cutoff, xcenter=xcenter, ycenter=ycenter, resolution=reso_arcmin)*core.G3Units.arcmin**2 # in input units-sr                          
    if invert_map:
        integrated_k_sr_bolo *= -1.
    cts_per_k = integrated_k_sr_bolo / integrated_k_sr_band # input-units / K               
    if data_is_wn:
        whitelevel = data1d
    else:
        npts_psd = np.min([1024,((len(data1d)-10)/4)*2])
        psdict = quick_pspec(data1d,rate=rate,npts_psd=npts_psd) # in input units / sqrt(Hz)
        whint = np.where(np.logical_and(psdict['freqs'] > frange[0],psdict['freqs'] < frange[1]))
        whitelevel = np.sqrt(np.mean(psdict['psd'][whint]**2))
    net = whitelevel / cts_per_k / np.sqrt(2.) # in K-sqrt(s)                               
    return net, cts_per_k, whitelevel, integrated_k_sr_bolo

#################################                                                          
def quick_nets_from_planet(mapdict, noisedict, banddict, planet_name, mjd, xcenter=None, ycenter=None, r_cutoff=5., reso_arcmin=0.5, rate=152.58789, frange=[10.,20.], data_is_wn=False, invert_map=False):

    from spt3g.calibration import planets

    """                                                                                     
    calculate NET from single-bolo map of a planet and some noise data                      
    (timestream). both should be in the same units (assuming ADC                            
    counts for now, but it doesn't matter as long as they're the                            
    same). result is NET in K_RJ-sqrt(s). update to use many maps at once.
    """
    
    ubands = np.unique(np.asarray(banddict.values()))
    integrated_k_sr_band = {}
    for uband in ubands:
        integrated_k_sr_band[str(uband)] = planets.temperature_solid_angle_product(planet_name,mjd,uband,TCMB=False)

    integrated_k_sr_bolo = {}
    cts_per_k = {}
    whitelevel = {}
    net = {}

    for name in mapdict.keys():
        integrated_k_sr_bolo[name] = int_src(mapdict[name], r_cutoff, xcenter=xcenter, ycenter=ycenter, resolution=reso_arcmin)*core.G3Units.arcmin**2 # in input units-sr                          
        if invert_map:
            integrated_k_sr_bolo[name] *= -1.
        cts_per_k[name] = integrated_k_sr_bolo[name] / integrated_k_sr_band[str(banddict[name])] # input-units / K               
        data1d = noisedict[name]
        if data_is_wn:
            whitelevel[name] = data1d
        else:
            npts_psd = np.min([1024,((len(data1d)-10)/4)*2])
            psdict = quick_pspec(data1d,rate=rate,npts_psd=npts_psd) # in input units / sqrt(Hz)
            whint = np.where(np.logical_and(psdict['freqs'] > frange[0],psdict['freqs'] < frange[1]))
            whitelevel[name] = np.sqrt(np.mean(psdict['psd'][whint]**2))
        net[name] = whitelevel[name] / cts_per_k[name] / np.sqrt(2.) # in K-sqrt(s)                               
    return net, cts_per_k, whitelevel, integrated_k_sr_bolo

#################################
def quick_pspec(data_1d, rate = 152.58789, npts_psd = 1024):

    npts = len(data_1d)
    npts_fft = npts_psd*2
    nchunk = npts/npts_fft
    win = np.hanning(npts_fft)
    winfac = np.mean(win**2)
    freqs = (makeFFTGrid(npts_fft,1./rate))[0:npts_psd]
    psd = np.zeros(npts_fft)
    for i in np.arange(nchunk):
        psd += (np.abs(np.fft.fft(data_1d[npts_fft*i:npts_fft*(i+1)]*win)))**2
    df = rate/2./np.float(npts_psd)
    psd = np.sqrt(psd/np.float(npts_fft)/np.float(npts_psd)/np.float(nchunk)/winfac/df)
    psd = psd[0:npts_psd]
    outdict = {}
    outdict['freqs'] = freqs
    outdict['psd'] = psd

    return outdict

#################################
def obsid_to_g3time(obsid_in, verbose=True):
    
    obsid = np.int(obsid_in)
    mjd0 = core.G3Time('01-Jan-2017:00:00:00').mjd
    time_out = core.G3Time()
    time_out.mjd = mjd0 + obsid/86400.

    if verbose:
        print(time_out)
    
    return time_out

#################################
def time_to_obsid(time_in):

    if type(time_in) is str:
        time1 = core.G3Time(time_in)
    mjd0 = core.G3Time('01-Jan-2017:00:00:00').mjd
    obsid_sec = np.int(np.round((time1.mjd - mjd0)*86400.))
    obsid_out = core.G3Int(obsid_sec)

    return obsid_out

#################################
def list_to_array(list_in):

    npts_tot = 0
    for item in list_in:
        npts_tot += len(item)
    array_out = np.zeros(npts_tot)
    npts = 0
    for item in list_in:
        array_out[npts:npts+len(item)] = item
        npts += len(item)

    return array_out

#################################
def list2d_to_array(list_in):

    nframes = len(list_in)
    nbolos = len(list_in[0])
    npts_tot = 0
    for item in list_in:
        npts_tot += len(item)
    array_out = np.zeros([nbolos,npts_tot])
    npts = 0
    for item in list_in:
        tnpts = len(item[0])
        for i in np.arange(nbolos):
            array_out[i,npts:npts+tnpts] = item[i]
        npts += tnpts

    return array_out

#################################
def ten_tc(sexagesimal_string_in, h2deg=False):

    ssi = sexagesimal_string_in
# check whether it really looks sexagesimal; if not, 
# just return decimal or decimal*15 (if h2deg)
    if ':' not in ssi:
        try:
            decimal_out = np.float(ssi)
            if h2deg:
                decimal_out *= 15.
        except:
            core.log_warn('Input string does not look like sexagesimal or decimal, giving up.')
            decimal_out = -1e32
        return decimal_out
    ssi2 = ssi.split(':')
    decimal_out = 0.
    gains = np.array([1., 1./60., 1./3600.])
    if ssi2[0][0] == '-':
        ssi2[0] = ssi2[0].replace('-','')
        gains *= -1.
    for num, gain in zip(ssi2,gains):
        decimal_out += np.float(num)*gain
    if h2deg:
        decimal_out *= 15.

    return decimal_out

#################################
def sixty_tc(decimal_in, deg2h=False, floatsec=False):

    dec2use = decimal_in
    if deg2h:
        dec2use = decimal_in/15.

    dd = np.abs(dec2use) 
    mm = np.abs(60.0*dec2use) 
    ss = np.abs(3600.0*dec2use)

    idd = np.int(dd)
    imm = np.int(mm-60.*idd)
    fss = ss-3600.*idd-60.*imm

    dds = '%02d' % idd
    mms = '%02d' % imm
    if floatsec:
        sss = '%06.3f' % fss
    else:
        sss = '%02d' % fss

    if decimal_in < 0.:
        dds = '-'+dds

    string_out = dds+':'+mms+':'+sss

    return string_out

#################################
def gcirc(ra1,dec1,ra2,dec2):

    ''' 
    great-circle distance between two points 
    on a sphere. assumes ra/dec in degrees, 
    outputs result in degrees. copied from 
    gcirc.pro in idlastro.
    '''
    
    d2r = np.pi/180.

# if input appears to be strings, assume 
# ra/dec in hhmmss / ddmmss, and convert
    racheck = ra1
    if np.size(racheck) > 1: 
        racheck = ra1[0]
    if type(racheck) is str:
        if np.size(ra1) > 1:
            for ratemp, dectemp in zip(ra1,dec1):
                rarad1 = ten_tc(ratemp,h2deg=True)*d2r
                dcrad1 = ten_tc(dectemp)*d2r
        else:
            rarad1 = ten_tc(ra1,h2deg=True)*d2r
            dcrad1 = ten_tc(dec1)*d2r
        if np.size(ra2) > 1:
            for ratemp, dectemp in zip(ra2,dec2):
                rarad2 = ten_tc(ratemp,h2deg=True)*d2r
                dcrad2 = ten_tc(dectemp)*d2r
        else:
            rarad2 = ten_tc(ra2,h2deg=True)*d2r
            dcrad2 = ten_tc(dec2)*d2r
    else:
        rarad1 = ra1*d2r
        rarad2 = ra2*d2r
        dcrad1 = dec1*d2r
        dcrad2 = dec2*d2r
    deldec2 = (dcrad2-dcrad1)/2.0
    delra2 =  (rarad2-rarad1)/2.0
    sindis = np.sqrt( np.sin(deldec2)*np.sin(deldec2) + 
                      np.cos(dcrad1)*np.cos(dcrad2)*
                      np.sin(delra2)*np.sin(delra2) )
    distance = 2.0*np.arcsin(sindis)/d2r

    return distance

#################################
def annotate_focal_plane():

    aptemp = {}
    aptemp['width'] = 2
    aptemp['color'] = 'k'

    xpos2 = [-37,-6,24,-73,-67,56,73,-23,6,37]
    ypos2 = [40,47,55,-28,14,38,0,-56,-48,-42]
   
    xpos2_text = (np.asarray(xpos2)).copy()
    ypos2_text = (np.asarray(ypos2)).copy()

    xpos2[4] = -31
    ypos2[4] = 1

    xpos2[5] = 26
    ypos2[5] = 15

#    print(xpos2)
#    print(ypos2)

#    wnames = ['W142','W148','W152','W139','W157','W147','W153','W136','W158','W162']
    wnames = ['W181','W180','W206','W204','W172','W177','W176','W203','W174','W188']

    for wname,xp,yp,xpt,ypt in zip(wnames,xpos2,ypos2,xpos2_text,ypos2_text):
#        if wname == 'W157' or wname == 'W147':
        if wname == 'W172' or wname == 'W177':
            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k',arrowprops=aptemp)
        else:
            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k')

    pylab.draw()
              
#def annotate_focal_plane():
#
#    aptemp = {}
#    aptemp['width'] = 2
#    aptemp['color'] = 'k'
#
#    xpos2 = [-37,-6,24,-66,-67,56,67,-23,6,37]
#    ypos2 = [40,47,55,-22,14,38,0,-50,-45,-37]
#   
#    xpos2_text = (np.asarray(xpos2)).copy()
#    ypos2_text = (np.asarray(ypos2)).copy()
#
#    xpos2[4] = -31
#    ypos2[4] = 1
#
#    xpos2[5] = 26
#    ypos2[5] = 15
#
#    print(xpos2)
#    print(ypos2)
#
#    wnames = ['W142','W148','W152','W139','W157','W147','W153','W136','W158','W162']
#
#    for wname,xp,yp,xpt,ypt in zip(wnames,xpos2,ypos2,xpos2_text,ypos2_text):
#        if wname == 'W157' or wname == 'W147':
#            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k',arrowprops=aptemp)
#        else:
#            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k')
#
#    pylab.draw()
              
#################################
def robust_sigma(array, nsthresh=3.):

    '''
    Calculate standard deviation, iteratively throwing out > N*sigma
    outliers. Default is N = 3.
    '''

    ftol = 1e-4
    npts = np.size(array)
    fnpts = np.float(npts)
    med0 = np.median(array)
    sd0 = np.sqrt(np.sum((array-med0)**2)/(fnpts-1.))
    sd = sd0
    md = med0
    array_cut = array.copy()
    tol = 1e6
    while tol > ftol:
        whgood = np.where(np.abs(array_cut-md) <= nsthresh*sd)
        array_cut = array_cut[whgood]
        if len(array_cut) == 0:
            break
        newmd = np.mean(array_cut)
        newsd = np.std(array_cut)
        tol = np.abs(newsd - sd)/sd
        sd = newsd
        md = newmd

    return sd

#################################
def average_rawdump_onedir_onecrate(dir1, crate1):

    '''
    average rawdump from all squids on a crate in one directory
    '''

    files = glob.glob(dir1 + 'IceCrate*' + str(np.int(crate1)) + '*IceBoard*.pkl')

    d1 = pickle.load(open(files[0]))
    freqs = d1['freq_domain']['x']
    rdave = np.zeros(len(freqs))

    for file1 in files:
        dtemp = pickle.load(open(file1))
        rdave += d1['freq_domain']['y']
        

    rdave /= np.float(len(files))

    return freqs, rdave

#################################
def get_noise_v_freq_two_files(file1, file2):

    '''
    get noise v. freq on the same set of bolos in two files. (if there
    are no overlapping bolos, this won't work.)
    '''

    d1 = pickle.load(open(file1))
    d2 = pickle.load(open(file2))

    freqs = []
    noise1 = []
    noise2 = []

    for key in d1.keys():
        if key in d2.keys():
            freqs.append(d1[key]['frequency'])
            noise1.append(d1[key]['noise']['median_noise'])
            noise2.append(d2[key]['noise']['median_noise'])

    freqs = np.asarray(freqs)
    noise1 = np.asarray(noise1)
    noise2 = np.asarray(noise2)

    return freqs, noise1, noise2

#################################
def trj_to_tcmb(freq):

    '''
    conversion factor between Rayleigh-Jeans temperature
    and CMB fluctuation temperature. input frequency
    should be in native spt3g_software units.
    '''

    tcmb = 2.728
    fghz = freq/core.G3Units.hz/1e9

    h_times_1e9_over_kb = (6.62607/1.38065)*1e-2
    x = h_times_1e9_over_kb*fghz/tcmb
    factor = (np.exp(x)-1.)**2/(x**2*np.exp(x))

    return factor

#################################
def camb2cl(cambfile = None):

    '''
    read CAMB power spectra from file into dictionary suitable for
    handing to quick_flatsky_routines.cmb_flatsky. assumes CAMB output
    is in D_l (l*(l+1)/2pi*C_l) and in units of uK^2.
    '''

    if cambfile is None:
        cambfile = '/home/tcrawfor/code/spt3g_software/simulations/camb/base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat'

    dl_all = np.genfromtxt(cambfile, unpack=True)
    cldict = {}
    ell_orig = dl_all[0,:]
    cl2dl = ell_orig*(ell_orig+1.)/2./np.pi
    ellmax = np.int(np.max(ell_orig))
    cldict['ell'] = np.arange(ellmax+1)
    cldict['cl'] = {}
    cltags = ['TT','EE','BB','TE']
    ellmin = np.int(np.min(ell_orig))
    for i in np.arange(4):
        tag = cltags[i]
        cldict['cl'][tag] = np.zeros(ellmax+1)
        cldict['cl'][tag][ellmin:] = dl_all[i+1,:]/cl2dl/1e6 # return Cls in native 3g units squared

    return cldict
