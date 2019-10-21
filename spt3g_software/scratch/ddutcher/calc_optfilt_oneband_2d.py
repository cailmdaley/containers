# This is a translation of spt_analysis/sources/calc_optband_2d.pro
#  which was made to:
#      "Calculate optimal filter for source detection in single map
#       with many different assumed source profiles."

import os
import numpy as np
from healpy import gauss_beam
from spt3g import core, mapmaker
from spt3g.mapspectra import map_analysis, apodmask, basicmaputils
from spt3g.util import files

spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']

def calc_optfilt_oneband_2d(profiles, skysigvar, noisepsd, reso_arcmin,
                            sigspect = 1,
                            fwhm_beam_arcmin = 1e-6,# 1E-6 is an infinitely small beam
                            ell_highpass = None,
                            ell_lowpass = None,
                            f_highpass = None, 
                            f_lowpass = None,
                            scanspeed = None, 
                            trans_func_grid= None, 
                            beam_grid = None,
                            keepdc = False, 
                            arcmin_sigma = None,
                            area_eff = None,
                            whitenoise1d = None,
                            pointsource = False,
                            norm_unfilt = None,
                            external_profile_norm = None,
                            bfilt = None,
                            debug1d = False,
                            debug = False,
                            ):
    '''
    This is a translation of spt_analysis/sources/calc_optband_2d.pro
    which was made to:

    "Calculate optimal filter for source detection in single map
       with many different assumed source profiles."

    Parameters:
    -----------
    profiles: (array) 1d ell-space profiles of objects you want to detect.
        This input is not used when pointsource = True. 

    skysigvar: (array) Assumed 1D variance as a function of ell (isotropic)
        of all sky signals other than the signal of
        interest.  (E.g., CMB + point sources for an SZ
        measurement.)

    noisepsd: (ndarray): N_pts-by-N_pts array of 2d Fourier
        representation of noise in map.
        Should be normalized such that
        sqrt(int_tabulated(uv_arr,noisepsd[*,*,j]^2)) =
        noise rms in the map.
        **** Need to confirm this normalization ****

    reso_arcmin: (float) Pixel size in arcminutes of maps to be filtered.

    sigspect[1] - sensitivity to signal of interest in supplied map.
          For clusters, this should be equal to fxsz(band).
          Leave this equal to 1 for single band filtering. 

    fwhm_beam_arcmin[1E-6]: (float) Size (in full-width at half-max in units of
          arcminutes) of beams in each band.  Beams
          assumed to be Gaussian.  Default is 10^-6
          (infinitely small beams).  Ignored if
          TRANS_FUNC_GRID is supplied.

    beam_grid - set this to a 2d array with the 2d beam response
          (if you want to specify a non-Gaussian or
          asymmetric beam but use default filters).

    bfilt[None]: (array) Hand the program a 1-d beam window function, 
          which will be used in place of the default assumed
          Gaussian beam for spatial filtering

    ell_highpass[None]: (float) Ell at which timestreams were high-passed
          before maps were made.  Default is no filtering.
          Ignored if TRANS_FUNC_GRID is supplied.
          *Do not supply both an ell_highpass and f_highpass.*

    ell_lowpass[None]: (float) Ell at which timestreams were low-passed
          before maps were made.  Default is no filtering.
          Ignored if TRANS_FUNC_GRID is supplied.
          *Do not supply both an ell_lowpass and f_lowpass.*

    f_highpass[None]: (float) Frequency at which timestreams were high-passed
          before maps were made.  Default is no filtering.
          Ignored if TRANS_FUNC_GRID is supplied.
          *Do not supply both an ell_highpass and f_highpass.*

    f_lowpass[None]: (float) Frequency at which timestreams were low-passed
          before maps were made.  Default is no filtering.
          Ignored if TRANS_FUNC_GRID is supplied.
          *Do not supply both an ell_lowpass and f_lowpass.*

    scanspeed[None]: (float) what speed were we scanning at when the maps were
          made?  (Need to know this to translate high- and
          low-pass temporal frequencies into spatial
          frequencies.)  Default is 0.25 (deg/s on the
          sky). 
          *This must be supplied to use f_highpass or f_lowpass*

    trans_func_grid[None]: (ndarray) Array the same size as MAP containing the
          effective delta-function response (in 2d
          fourier space) of all the filtering we've
          done to the signal in the map
          (INCLUDING THE BEAM).  If TRANS_FUNC_GRID
          is set, other filtering keywords (including
          FWHM_BEAM_ARCMIN) are ignored.

    keepdc[False]: (bool) Set this to keep the DC information along the scan
          direction.  Otherwise, the u=0 modes
          are set to zero even if F_HIGHPASS isn't
          set. 

    arcmin_sigma - SIGMA scaled to what it would be for a
          1-arcminute FWHM beam.  (SIGMA is reported for
          the smallest beam used.)

    area_eff - effective area of sources with input profile and
          unit amplitude, convolved with smallest beam.

    whitenoise1d[False]: (bool)Set this to ignore all anisotropic noise and
          filtering (mostly for debugging). 
          ****** This is currently not working*******

    pointsource[False]: (bool)Set this to True to tell the routine that you want the
          optimal filter for a point source. This effectively makes the
          profile the beam scale. Setting pointsource=True is equivalent
          to setting profiles = np.ones(ellmax)

    norm_unfilt - set this to normalize the output such that 
          sources of the shape of one of the input
          profiles will have amplitude in the output map
          equal to the amplitude of the input source
          before filtering or beam-smoothing.

    external_profile_norm - set this to the known area (in sr)
          under the real-space version of your
          source profiles, if you don't want the
          area calculated in ell space.  Only
          appropriate if NORM_UNFILT is set.

    Returns:
    --------
    2D Fourirer-space optimal filters, one for each input profile

    HISTORY:
       Big Bang + 300,000 years: Tom Crawford writes the initial IDL function. 
       Dec 15th, 2014: Tyler Natoli translates the function to python. 
    '''

    ellmax = len(skysigvar)-1 # 
    ell = np.arange(0,ellmax+1)

    ngridy, ngridx = np.shape(noisepsd)

    if pointsource:
        # This makes the optimal filtering for a point source
        #   which just makes the profile an array of 1s for now.
        #   later the beam will be multiplied into the profile, 
        #   which means the profile will effectively just be the beam sky
        profiles = np.ones(ellmax +1) 

    if external_profile_norm:
        ##################### I dont really understand what this does at all
        do_ext_norm = True

    if f_highpass or f_lowpass:
        if not scanspeed:
            print('You must supply a scanspeed to use f_highpass or f_lowpass')
            print('....but you did not, so I am returning without doing a damn thing.')
            return
    
    # just making sure this not an integer
    reso_arcmin = float(reso_arcmin)
    reso_rad = reso_arcmin * core.G3Units.arcmin/core.G3Units.rad

    # get a 2d grid of values appropriate fo 2d fft map
    ellg = basicmaputils.make_ellgrid(
        reso_rad*core.G3Units.rad, np.shape(noisepsd))
    
    # create a flat ellg with values above ellmax set to ellmax
    #     this will be used later to make a 2d things in Fourier space
    flat_ellg = np.round(ellg.flatten()).astype(int)
    flat_ellg[flat_ellg>=ellmax] = ellmax -1
    
    # array of residual noise to use in optimal filter &
    #   variance calculation
    ivarr = np.zeros(np.shape(noisepsd))+1e20
    
    # make a gaussian beam to use if a beam is not supplied
    if bfilt == None:
        bfilt = gauss_beam(fwhm_beam_arcmin*core.G3Units.arcmin/core.G3Units.rad,
                           lmax = ellmax)

    # get 2d Fourier representations of timestream filters and beams
    # (if TRANS_FUNC_GRID is set, these will be ignored except for 1d
    #     debugging case.)
    #
    # trans_func_grid is the 2D transfer function (beam+filtering+etc..) in ell-space
    if trans_func_grid == None:
        ellg_to_use = ellg[0,:] # original
        #ellg_to_use = ellg[0,:]
 
        trans_func_grid = np.zeros(np.shape(noisepsd))
        hplp1d = np.ones(len(ellg[0,:]))# high-pass low-pass 1 dimmensional
        
        if f_highpass or f_lowpass:
            # the effective frequency is used to turn the high & low pass
            # frequency numbers into ell ranges
            f_eff = ellg[0,:]*(scanspeed*(core.G3Units.deg/core.G3Units.rad)/
                               2./np.pi) # f_eff is the effective frequency 
            if f_highpass:
                hplp1d *= np.exp(-1.*(f_highpass/f_eff)**6)
            if f_lowpass:
                hplp1d *= np.exp(-1.*(f_eff/f_lowpass)**6)

        if ell_highpass or ell_lowpass:
            if ell_highpass:
                hplp1d *= np.exp(-1.*(ell_highpass/ellg_to_use)**6)
            if ell_lowpass:
                hplp1d *= np.exp(-1.*(ellg_to_use/ell_lowpass)**6)

        if not keepdc:
            # kill the dc component
            hplp1d[ellg_to_use == 0] = 0 # ell=0 is the dc component

        for _num in np.arange(0, ngridy):
            trans_func_grid[_num,:] = hplp1d
        if beam_grid != None:
            trans_func_grid *= beam_grid
        else:
            # multiply the beam into the trans_func_grid at every ell in ellg
            trans_func_grid *= np.reshape(bfilt[flat_ellg],np.shape(ellg))
    
    # get 2D ell-space representation of unwanted 1D spectra
    skysigvar2d = np.reshape(skysigvar[flat_ellg],np.shape(ellg))# skysigvar is 1D, make it 2D
    # skysigvar2d is in power, so multipy it by the trans function squared
    skysigvar2d *= abs(trans_func_grid)**2

    # calculate noise + cotaminant vaiance at all grid points
    skysigvar2d += noisepsd**2

    # ivarr is the noise spectra (cmb+point source+pixel noise)
    ivarr = skysigvar2d/sigspect**2 # Im not really sure what sigspect is, 
    szbftemp = abs(trans_func_grid*sigspect)**2 # sz beam temp?
    
    whsig = [szbftemp/np.max(szbftemp) > 1e-4][0]# where the values are over 1/10,000
    nsig = whsig.tolist().count(True)
    # only use ell values up to where the transfer function has values over 1/10,000
    ellgmax = max(abs(ellg[whsig]))
    ellmax2use = ellgmax
    if ellgmax < ellmax:
        ellmax2use = ellgmax
    ell2use = np.arange(0, ellmax2use+1) # these are the good ells
    
    # 1d case for debugging  - interp(xneed,xhave,yhave)
    if debug1d:
        noisepsd1d = np.interp(ell2use,ellg[0,0:ngridx/2.],noisepsd[0,0:ngridx/2])
        whgellg = np.where(ell2use > max(ellg[0, 0:ngridx/2.]))[0]
        ngellg = len(whgellg)
        if ngellg > 0:
            noisepsd1d = noisepsd1d[:whgellg[0]-1]
        ivarr1d = (skysigvar[:int(ellmax2use)]*bfilt[:int(ellmax2use)]**2 + noisepsd1d**2)/sigspect**2 
        ivarr1d[0] = 1e20
        #  note: ivarr1d only has the beam in the transfer function
        
    ivarr[0,0] = 1e20# kill the dc component (the first bin)


#----------------------------------------------------------------
    # up until this point all we have made is the noise and unwanted signal convolved with the 
    #   transfer function, now we will start making the optimal filter
#----------------------------------------------------------------
    # make optimal filter(s) to use on synthesized map (using
    # Haehnelt-Tegmark prescription).
    try:
        nprofs, profile_len = np.shape(profiles) # nprofs = the number of profiles            
    except ValueError:# above will crash if profiles is a single array
        nprofs = 1
    if whitenoise1d:
        optfilts = np.zeros([ellmax2use+1.,nprofs],dtype = np.complex)
    else:
        optfilts = np.zeros([ngridy,ngridx,nprofs],dtype = np.complex)

    results = {}

    for _nprof in range(0,nprofs):
        if nprofs==1:
            prof_now = profiles
        else: 
            prof_now = profiles[_nprof]
       
        results[_nprof] = {'profile':prof_now}

        # If we want to get the 1D case working:
        # we must take extra care to ensure the profile, beam, and transfer function 
        #    are all the same length, they will be truncated to whatever the shortest
        #    one of those variables is
        #max(len(prof_now),len(bfilt))
        
        if debug1d:
            tauprime = prof_now*bfilt[:len(prof_now)] # tauprime is the profile*beam
            if norm_unfilt: # this sets the normalization....not really sure why we do this
                if do_ext_norm:
                    area_eff_1d = external_profile_norm[_n_prof]
                else:
                    area_eff_1d = 1./np.sum((2.*ell+1)/4./np.pi*prof_now)
            else:
                area_eff_1d = 1./np.sum((2.*ell+1)/4./np.pi*tauprime)
            tauprime = tauprime*area_eff_1d # 1d version of tau*filtering*beam
        

        area_eff_1am = 1./np.sum((2.*ell+1)/4./np.pi*prof_now*gauss_beam(
            1*core.G3Units.arcmin/core.G3Units.rad, lmax = ellmax))
        # make 2d versions np.reshape(skysigvar[flat_ellg],np.shape(ellg))
        prof2d = np.reshape(prof_now[flat_ellg],np.shape(ellg)) # takes the 1d profile and makes it 2d
        tauprime_2d = prof2d*trans_func_grid # this is multipying the beam*transfer_func into the profile, both are in ell space
        du = (ellg[0,1] - ellg[0,0])/2./np.pi # du = dx -- although it really does not matter for this
        dv = (ellg[1,0] - ellg[0,0])/2./np.pi # dv = dy ---- anytime we use a du we use dv also

        if norm_unfilt:# this sets the normalization....not really sure why we do this
            if do_ext_norm:
                area_eff = external_profile_norm
            else:
                area_eff = 1./np.sum((2.*ell+1.)/4./np.pi*prof_now)
        else:
            area_eff = 1./np.sum(abs(tauprime_2d)*du*dv)

        tauprime_2d = tauprime_2d*area_eff

        optfilt = tauprime_2d/ivarr
        optfilt_norm = np.sum((np.conj(tauprime_2d)*optfilt).real)*du*dv
        optfilt = optfilt/optfilt_norm
        
        # 1d case for checking things and debugging
        if debug1d:
            optfilt1d = tauprime/ivarr1d
            optfilt1d_norm = np.sum(2*(ell+1)/4./np.pi*optfilt1d*tauprime)
            optfilt1d = optfilt1d/optfilt1d_norm
            sigma_1d = 1./sqrt(optfilt1d_norm)
            arcmin_sigma_1d = sigma_1d*area_eff_1d/area_1am                             
        # calculate variance of central signal value of features detected in
        # synthesized map using optimal filter (this is the variance on the
        # central value of the features after they have been smoothed by the
        # smallest beam among the bands).
        sigma = 1./np.sqrt(optfilt_norm)
        # scale to that same value if the minimum beam were 1 arcminute
        arcmin_sigma = sigma*area_eff/area_eff_1am
        # stuff values in results object
        if whitenoise1d: # do something different with this
            optfilts[:,_nprof] = optfilt1d
            sigmas[_nprof] = sigma_1d
            arcmin_sigmas[_nprof] = arcmin_sigma_1d
            area_effs[_nprof] = area_eff_1d
        
        results[_nprof]['optfilt'] = optfilt
        results[_nprof]['sigma'] = sigma
        results[_nprof]['arcmin_sigma'] = arcmin_sigma
        results[_nprof]['area_eff'] = area_eff
        results[_nprof]['trans_func_grid'] = trans_func_grid
        results[_nprof]['ell_grid'] = ellg
        results[_nprof]['beam_profile'] = bfilt
        results[_nprof]['skysigvar'] = skysigvar2d

    return results


def call_matched_filt_point(map_fr,
                            cmb = True,
                            ellmax = 30000,
                            skysigvar = None,
                            white_noise_level = None,
                            noisepsd = None, 
                            beam = None,
                            fwhm_beam_arcmin = 1.2,
                            ell_highpass = 300,
                            ell_lowpass = 20000,
                            apod_mask = None,
                            arcmin_sigma = None,
                            norm_unfilt = None,
                            external_profile_norm = None,
                            return_filtered_map = False,
                            verbose = False,
                            ):
    '''
    This function is just a wrapper for calc_optfilt_oneband_2d above. 
        
    You feed this function a map and it extracts parameters 
       to make the matched source filter.

    This code assumes your profile is for a point source. 

    Parameters:
    -----------
    map_fr: a Map Frame containing a temperature FlatSkyMap.
        If white_noise_level = True, then the white noise level will be the
        average(map_psd[6000 < ell < 10000]) of this map. 

    cmb [True]: (bool or file)
        If True use the Planck cmb spectra (base_plikHM_TTTEEE_lowl_lowE_lensing)
        in the unwanted map signal.
        If cmb = file, it is assumed that file is a camb
        file and will use this file as the cmb spectra in the unwanted
        map signal

    white_noise_level [None]: (bool or float)
        If a float, the float given  will be the white noise level given
        to the filter maker. If True, then the white_noise_level will be taken from
        the Map supplied and will be the average(map_psd[4000 < ell <6000]).
        In uK-arcmin

    skysigvar [None]: 

    plot_map [False]: (bool)

    debug [False]: (bool)

    returned_filtered_map [False]: (bool)
        Returns the m_filt dict as well as a map frame where the maps
        have been filtered (but only the T map has been properly filtered). 

    verbose [False]: (bool)
        Print extra information to the screen.

    Returns:
    --------
    The optimal filter for a point source and the noise in the given map. If 
    return_filtered_map = True then the filtered map object is also returned.
    '''
    wgt = None
    if isinstance(map_fr, str):
        obs = core.G3File(map_fr)
        for map_fr in obs:
            if map_fr.type == core.G3FrameType.Map:
                break
    if isinstance(map_fr, core.G3Frame):
        for k in map_fr.keys():
            if 'W' in k and 'pol' in k:
                wgt = map_fr[k]
        if map_fr['T'].is_weighted:
            if wgt is None:
                raise KeyError("'Wunpol' or 'Wpol' not found in frame.")
            Tmap = mapmaker.mapmakerutils.remove_weight_t(map_fr['T'], wgt)
        else:
            Tmap = map_fr['T']
    else:
        raise TypeError('Input must be a Map frame')
        
    if isinstance(apod_mask, str):
        if apod_mask == 'default':
            apod = apodmask.make_border_apodization(
                wgt, apod_type='cos', radius_arcmin=90.)
            ptsrc = apodmask.make_apodized_ptsrc_mask(
                Tmap, os.path.join(spt3g_software,'scratch/ddutcher',
                                   '1500d_3band_10sigma_ptsrc.txt'))
            apod_mask = apod
        else:
            apod_mask = files.load_pickle(apod_mask)
            
    # get the resolution from the T map
    reso_arcmin = Tmap.res/core.G3Units.arcmin
    reso_rad = Tmap.res/core.G3Units.rad

    if skysigvar is not None:
        ellmax = len(skysigvar)-1
    ell = np.arange(0,ellmax)

    if skysigvar is None:
        skysigvar = np.zeros(ellmax)
        if cmb != False:
            if cmb == True:
                cl_file = os.path.join(spt3g_software,
                    'simulations/camb/base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat')
                cls = files.read_camb_cls(cl_file, as_cls = True, extrapolate = True, 
                                    max_ell = ellmax-1)
            elif type(cmb) == str:
                cls = files.read_camb_cls(cmb, as_cls = True, extrapolate = True, 
                                    max_ell = ellmax-1)
            else:
                # assume cmb is a spectra that has already been read in and is the 
                #   standard form ---------------- make this more general
                cls = cmb
    
        skysigvar[cls['L']] = cls['TT']

    if white_noise_level:
        if white_noise_level==True:
            # get the noise level from the input map
            psd = map_analysis.calculateCls(map_fr, apod_mask = apod_mask*ptsrc, t_only=True,
                                            ell_max = 6000)
            psd['TT'] /= (core.G3Units.arcmin * core.G3Units.uK)**2
            white_noise_level = np.sqrt(
                np.mean(psd['TT'][np.where((psd['ell']>4000) &
                                           (psd['ell']<6000))[0]]))
            if verbose:
                print('White noise level is %.2f uK-arcmin'%white_noise_level)     

        # if white_noise_level is not a bool, assume it is a float of the noise level
        noisepsd = np.ones(np.shape(Tmap)) * white_noise_level

    elif white_noise_level==False and noisepsd==False:
        raise ValueError(
            'You must either set white_noise_level or supply a noisepsd.')
    if beam:
        # not dealing with beam files yet
        bfilt = None
#         try:
#             # if we were given a beam already read in this will work
#             beam_1d_ell = beam['B_ell']
#         except TypeError:
#             # we will now read the beam in 
#             beam_1d_ell,beam_1d = files.read(beam)['B_ell']
#         #get the beam in the same ells as everything else
#         bfilt = np.interp(ellx[0][:len(ellx)/2],beam_1d_ell,beam_1d) #interp(xneed,xhave,yhave)
    else:
        bfilt = None

    m_filt =calc_optfilt_oneband_2d([0,0],# profiles is not used with pointsource=True
                                    skysigvar, noisepsd, reso_arcmin,
                                    fwhm_beam_arcmin = fwhm_beam_arcmin, 
                                    ell_highpass = ell_highpass,
                                    ell_lowpass = ell_lowpass,
                                    bfilt = bfilt,
                                    # keeping below for all calls from the function
                                    pointsource = True,
                                    sigspect = 1,
                                    keepdc = False, 
                                    arcmin_sigma = None,
                                    norm_unfilt = None,
                                    external_profile_norm = None)

    if return_filtered_map:
        filt_map = map_analysis.filterMap(map_fr, m_filt[0]['optfilt'], apod_mask=apod)
        filt_map['filt_sigma']=m_filt[0]['sigma']
        filt_map['source_filtered_map']=True
#         if plot_map:
#             filt_map.pol_maps.T.drawImage(figure=1,mask=mask)
#             filt_map.pol_maps.T.drawImage(figure=2,mask=mask*(1./m_filt[0].sigma),)
#             pl.suptitle('units in S/N, not K_CMB')
#             plotting.kspaceImshow(m_filt[0].optfilt,reso_arcmin=reso_arcmin,figure=4,log=False,title='ell-space optfilt')
#             pl.figure(5);pl.clf()
#             plotting.plotRadialProfile(np.fft.fftshift(m_filt[0].optfilt), center_bin=[ngridx/2,ngridy/2], 
#                                        resolution=(ellx[0,1]-ellx[0,0]), 
#                                        histogram=True, range=[0,20000],
#                                        plot=True, bin_resolution_factor=0.1, 
#                                        figure=5,new=False)
#             pl.title('1D optfilt')
#             pl.xlabel('ell')
    if verbose:
        print('sigma of map = ',m_filt[0]['sigma'])

    if return_filtered_map:
        return m_filt, filt_map
    return m_filt
