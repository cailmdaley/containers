import numpy as np
from spt3g import core, coordinateutils, mapmaker
from spt3g.util.fitting import fit_gaussian2d, gaussfit_hist
from spt3g.util.maths import gaussian2d
from spt3g.sources.source_utils import find_sources_quick
from spt3g.util import files
import os
from scipy import ndimage

def check_astrometry_at20(map_in, pixel_mask=None, weight=None, 
                          gauss_fit_switch=False, 
                          field='tmp_field', nsigma=10, 
                          plot=True, verbose=False,
                          check_beam_size = True,
                          show_source_positions=True, 
                          rcrop=None, sf=30.,
                          close_cut_arcsec=120., 
                          use_external_list=False,
                          ra_external=None, 
                          dec_external=None, 
                          sn_external=None, 
                          power_weight=0.):

    '''
    This program looks for very bright point sources in a map and
    compares their RA/DEC to nearby sources in the AT20 catalog.
    
    Parameters
    ----------
    map_in: Either an unweighted T map (FlatSkyMap) or a Map frame containing
        T and weights OR a list of several such objects
    pixel_mask [None]: Mask applied to map before source-finding. If map_in
        contains a list of thumbnails, and pixel_mask is defined, it is
        assumed that the same pixel mask should be used for all maps.
    weight [None]: If pixel_mask is not provided, these weights will be
        used to create one. If neither is provided, the code will look
        in the map object(s) for weight arrays to use.
    nsigma [10.]: Minimum signal-to-noise required to detect a source
    check_beam_size [True]: Use brightest point sources in map to make
        measurement of instrument beam.
    rcrop [None]: Sets the size of the map cutout around the sources used
        for beam fitting, in arcmin. If not set, defaults to 11x11 pixels.
    use_external_list [False]: If True, does not search map for sources,
        instead using locations provided via ra_external, dec_external.
    close_cut_arcsec [120.]: How far away a detected source can be from a
        catalog source and still be matched with it.
    plot [True]: Generate plot of coordinate differences from map sources
        and catalog sources.
    sf [30.]: The min/max range of the coord difference plot, in arcsec.

    Output:
    -------
    Dictionary with keys 'ra','dec','sn','dra','ddec','dxdec','dr',
        corresponding to matched source coordinates, SN in the map and the
        differences between those coords and those in the catalog. 

    Notes:
    ------
    Created Mar2010, RK.
    Jan2019 DPD: Ported to python, spt3g_software from 
        spt_analysis/pointing/check_astrometry_at20.pro
    Aug2019 TC: Hacked to accept lists of source-centered
        thumbnail maps.
    '''

    if plot:
        import matplotlib.pyplot as plt

    dtor = core.G3Units.deg/core.G3Units.rad
    try:
        spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
    except KeyError:
        spt3g_software = os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))
    
    results = dict()
    
    if not use_external_list:

        ra = []
        dec = []
        xpeak = []
        ypeak = []
        sn = []
        sn_all = [] # for debugging

        if isinstance(map_in, list):
            maplist = map_in
            if check_beam_size and gauss_fit_switch:
                fwhms = []
                avg_fwhm = []
                ecc = []
        else:
            maplist = [map_in]

        for map in maplist:

            if isinstance(map, coordinateutils.FlatSkyMap):
                if (map.is_weighted or map.pol_type != coordinateutils.MapPolType.T):
                    raise ValueError("Input map must be unweighted T")
                else:
                    Tmap = map
            elif isinstance(map, core.G3Frame):
                for key in map:
                    if isinstance(map[key],coordinateutils.G3SkyMapWeights):
                        w = map[key]
                if weight is None: 
                    weight_from_map = np.asarray(w.TT)
                    weight = weight_from_map
                if map['T'].is_weighted:
                    Tmap = mapmaker.mapmakerutils.remove_weight_t(map['T'], w)
                else:
                    Tmap = map['T']
            else:
                raise ValueError("Input %s not supported; "%type(map)+
                                 "please use a Map frame or FlatSkyMap")
    
            thisxpeak = []
            thisypeak = []

            reso_arcmin = Tmap.res/core.G3Units.arcmin
            nx=Tmap.shape[1]
            ny=Tmap.shape[0]
            #         npixels=[nx,ny]
        
            # get pixel mask
            if pixel_mask is None:
                if weight is None:
                    print('You must provide a pixel mask or weight mask.')
                    print('Quitting!')
                    return

                wh_non_zero = np.where(weight != 0.)
                med_w = np.median(weight[wh_non_zero])
                whmask = np.where(weight > 0.5*med_w)

                pixel_mask=np.zeros(Tmap.shape)
                pixel_mask[whmask]=1
            
            if weight is None:
                weight = np.zeros(np.shape(pixel_mask))
                wh_non_zero = np.where(pixel_mask != 0.)
                weight[wh_non_zero] = 1
            # find bright point sources
            if verbose:
                print('...finding point sources...')
                print('...NSIGMA = %s'%nsigma)

            whn0 = np.where(pixel_mask > 0.)
            if gauss_fit_switch:
                ntmap = np.asarray(Tmap)
                if np.max(np.abs(ntmap)) > 0.:
                    ntmap -= np.mean(ntmap[0:10,0:10])
                    params = fit_gaussian2d(ntmap*pixel_mask,fit_offset=True)
                    fitfun = gaussian2d(*params)
                    ntmap_resid = ntmap - fitfun(*np.indices(ntmap.shape))
                    sigma_best = np.sqrt(np.abs(params[3]*params[4]))
                    mapsm = ndimage.gaussian_filter(ntmap,sigma_best)
                    mapsm_resid = ndimage.gaussian_filter(ntmap_resid,sigma_best)
                    temprms = np.std(mapsm_resid[whn0])
                    atemp = gaussfit_hist(mapsm_resid[whn0],1000,-8.*temprms,8.*temprms,
                                          do_plot = False)
                    maprms = atemp['params'][2]
                    try:
                        thissn = np.abs(mapsm[np.int(np.round(params[2])),np.int(np.round(params[1]))])/maprms
                        sn_all.append(thissn)
                        fwhms_temp = params[3:5]*reso_arcmin*(2*np.sqrt(2*np.log(2)))
                        if thissn > nsigma and np.min(fwhms_temp) > 0.5 and np.max(fwhms_temp) < 3.:
                            thisxpeak.append(params[1])
                            thisypeak.append(params[2])
                            xpeak.append(params[1])
                            ypeak.append(params[2])
                            sn.append(thissn)
                    except:
                        sn_all.append(0.)
                else:
                    sn_all.append(0.)
                if check_beam_size:
                    thesefwhms = params[3:5]*reso_arcmin*(2*np.sqrt(2*np.log(2)))
                    fwhms.append(thesefwhms)
                    avg_fwhm.append(np.sqrt(np.abs(thesefwhms[0]*thesefwhms[1])))
                    ecc.append(max(thesefwhms)/min(thesefwhms))
     

            else:
                temprms = np.std(np.asarray(Tmap)[whn0])
                atemp = gaussfit_hist(np.asarray(Tmap)[whn0],1000,-8.*temprms,8.*temprms,
                                      do_plot = False)
                maprms = atemp['params'][2]
                s_out = find_sources_quick(Tmap, pixel_mask=pixel_mask,nsigma=nsigma,
                                           reso_arcmin = reso_arcmin, maprms=maprms)

                if s_out is None and len(map_in) == 1:
                    print('CHECK_ASTROMETRY: no bright point sources were found.')
                    print('Quitting!')
                    return
        
                if s_out is not None:
                    for source,d in s_out.items():
                        thisxpeak.append(d['xpeak'])
                        thisypeak.append(d['ypeak'])
                        xpeak.append(d['xpeak'])
                        ypeak.append(d['ypeak'])
                        sn.append(d['peaksig'])

            for i in np.arange(len(thisxpeak)):
                coords = Tmap.xy_to_angle(
                    np.float(thisxpeak[i]), np.float(thisypeak[i]))/core.G3Units.deg
                tmp_ra = coords[0] if coords[0]>0 else coords[0]+360.
                ra.append(tmp_ra)
                dec.append(coords[1])
                
        if plot:
            if show_source_positions:
                plt.figure(num=1,figsize=(10,6))
                plt.imshow(Tmap*pixel_mask, origin='lower', cmap=plt.cm.gray,
                           vmin = -5.*maprms,
                           vmax = 5.*maprms)
                plt.plot(xpeak, ypeak,'ro',markersize=12, mfc='None')


        xpeak = np.array(xpeak)
        ypeak = np.array(ypeak)
        sn = np.array(sn)
        ra = np.array(ra)
        dec = np.array(dec)
            
        
        n_sources = len(xpeak)
        print(str(n_sources)+' bright point sources were found.')




        # while we're at it, let's get the fwhm of the brightest source(s)
        if check_beam_size:
            nind = nsources if (n_sources < 3) else 3
            ind = np.arange(nind)
            if gauss_fit_switch:
                ssn = np.argsort(sn)
                fwhms = np.asarray(fwhms)[ssn[-n_sources:],:]
                avg_fwhm = np.asarray(avg_fwhm)[ssn[-n_sources:]]
                ecc = np.asarray(ecc)[ssn[-n_sources:]]
            else:
                if rcrop is None:
                    ncrop = 11
                else:
                    ncrop = np.int(round(rcrop/reso_arcmin))
                fwhms = np.zeros((nind,2))
                avg_fwhm = np.zeros(nind)
                ecc = np.zeros(nind)
                for kk in np.arange(nind):
                    xc=xpeak[ind[kk]]
                    yc=ypeak[ind[kk]]
                    minx=int(max([0,xc-ncrop]))
                    maxx=int(min([nx-1,xc+ncrop]))
                    miny=int(max([0,yc-ncrop]))
                    maxy=int(min([ny-1,yc+ncrop]))
                    submap = np.array(Tmap)[miny:maxy,minx:maxx]
                    subweight = weight[miny:maxy,minx:maxx]
                    fit = fit_gaussian2d(submap)
                    fwhms[kk,:] = fit[3:5]*reso_arcmin*(2*np.sqrt(2*np.log(2)))
                    avg_fwhm[kk] = np.sqrt(np.abs(fwhms[kk,0]*fwhms[kk,1]))
                    ecc[kk] = max(fwhms[kk,:])/min(fwhms[kk,:])

            results['fwhms'] = fwhms
            results['avg_fwhm'] = avg_fwhm
            results['ecc'] = ecc

    else:
        ra = ra_external
        dec = dec_external
        n_sources = len(ra)
        if sn_external is None:
            sn_external = np.zeros(n_sources) + 10.
        sn = sn_external

    # compare to the AT20 catalog
    if verbose:
        print('...loading AT20 catalog...')
    at20g = files.load_pickle(os.path.join(spt3g_software,
                                           'sources','at20g.pkl'))
    close_cut_deg=close_cut_arcsec/3600.

    dra=np.zeros(n_sources)
    dxdec=np.zeros(n_sources)
    ddec=np.zeros(n_sources)
    dr=np.zeros(n_sources)-1
    for i in np.arange(n_sources):
        this_dra=ra[i]-at20g['ra']
        this_ddec=dec[i]-at20g['dec']
        this_dxdec=this_dra*np.cos(dec[i]*dtor)
        this_dr=np.sqrt(this_dxdec**2. + this_ddec**2.)
        wh_close=np.where(this_dr <= close_cut_deg)[0]
        if len(wh_close) == 0:
            continue
        if len(wh_close)>= 2:
            wh_close=np.where(this_dr == min(this_dr))[0]

        wh_close=wh_close[0]
        dra[i]=ra[i]-(at20g['ra'])[wh_close]
        ddec[i]=dec[i]-(at20g['dec'])[wh_close]
        dxdec[i]=dra[i]*np.cos(dec[i]*dtor)
        if verbose:
            print(i,dxdec[i]*3600.,ddec[i]*3600.)
    
    dr=np.sqrt(dxdec**2. + ddec**2.)

    whok=np.where(dr > 0.)[0]
    if np.size(whok[0]) == 0:
        print('CHECK_ASTROMETRY: There were no matches to the AT20 catalog.')
        print('Quitting!')
        return

    dra=dra[whok]
    dxdec=dxdec[whok]
    ddec=ddec[whok]
    dr=dr[whok]
    ra_out=ra[whok]
    dec_out=dec[whok]
    sn_out = sn[whok]

    if plot:
        if show_source_positions:
            if not use_external_list:
                plt.figure(1)
                plt.plot(xpeak[whok],ypeak[whok],'bo',markersize=12, mfc='None')

    # plot stuff
    # it would be good to plot realilstic error bars (quadrature sum of
    # SPT and AT20)
    if plot:
        fig, ax = plt.subplots()
#         window,/free,xs=600,ys=500
        ax.plot(dxdec*3600.,ddec*3600.,'k+')
        ax.set_xlim(sf*np.array([-1,1]))
        ax.set_ylim(sf*np.array([-1,1]))
        ax.set_aspect('equal')
                 
        ax.set_xlabel('(SPT_RA - AT20_RA)*cos(dec) [arcsec]') 
        ax.set_ylabel('(SPT_DEC - AT20_DEC) [arcsec]')
        ax.set_title('SPT Astrometry for '+field)

    # created a weighted average for these offsets.
    # the error should be dominated by SPT's
    # error, which ~1/(SPT_S/N), so we create weights that are
    # proportional to SPT_S/N^2.
    # the SPT positional uncertainty should be sigma ~ FWHM/(S/N).  RK has
    # tested that sigma = 60/(S/N) is quite accurate (better than 20%) for
    # SPT 150 GHz.  The uncertainty at 90 GHz (220 GHz) will be worse
    # (better) for fixed S/N.
    err_arcsec = 60./sn_out
    w = err_arcsec**(power_weight)
    avg_dxdec = np.sum(dxdec*w)/np.sum(w)
    avg_dra = np.sum(dra*w)/np.sum(w)
    avg_ddec = np.sum(ddec*w)/np.sum(w)

    err_dxdec = (np.sum(err_arcsec**(-2.)))**(-0.5)/3600.
    err_ddec = (np.sum(err_arcsec**(-2.)))**(-0.5)/3600.
    dec0 = Tmap.delta_center/core.G3Units.deg
    err_dra = err_dxdec/np.cos(dec0*dtor)

    if verbose:
        print(' ')
        print('SUBTRACT these numbers (D_RA, D_DEC) '+
              'from the nominal map centers:')
        print('<D_RA> = ',str(avg_dra*3600.),' +/- ',
              str(err_dra*3600.),' arcsec')
        print('<D_XDEC> = ',str(avg_dxdec*3600.),' +/- ',
              str(err_dxdec*3600.),' arcsec')
        print('<D_DEC> = ',str(avg_ddec*3600.),' +/- ',
              str(err_ddec*3600.),' arcsec')
        print(' ')

    results['ra'] = ra_out
    results['dec'] = dec_out
    results['sn'] = sn_out
    results['dra'] = dra
    results['ddec'] = ddec
    results['dxdec'] = dxdec
    results['dr'] = dr
    try:
        results['sn_all'] = sn_all
    except:
        pass

    return results
