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
def annotate_focal_plane():

    aptemp = {}
    aptemp['width'] = 2
    aptemp['color'] = 'k'

    xpos2 = [-37,-6,24,-66,-67,56,67,-23,6,37]
    ypos2 = [40,47,55,-22,14,38,0,-50,-45,-37]
   
    xpos2_text = (np.asarray(xpos2)).copy()
    ypos2_text = (np.asarray(ypos2)).copy()

    xpos2[4] = -31
    ypos2[4] = 1

    xpos2[5] = 26
    ypos2[5] = 15

    print(xpos2)
    print(ypos2)

    wnames = ['W142','W148','W152','W139','W157','W147','W153','W136','W158','W162']

    for wname,xp,yp,xpt,ypt in zip(wnames,xpos2,ypos2,xpos2_text,ypos2_text):
        if wname == 'W157' or wname == 'W147':
            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k',arrowprops=aptemp)
        else:
            pylab.annotate(wname,(xp,yp),xytext=(xpt-10,ypt-4),xycoords='data',size='x-large',color='k')

    pylab.draw()
              
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

