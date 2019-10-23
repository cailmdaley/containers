'''
source_utils.py

Functions for finding, analyzing (point) sources in maps.
'''
import os
import numpy as np
from spt3g import core, coordinateutils, mapmaker
from spt3g.util.fitting import gaussfit_hist
from spt3g.util import files

def parse_point_source_file(filename, return_flux=False):
    '''
    Reads the spt3g/sptpol/sptsz point source file format.
    
    Arguments:
    -----------
    
    filename: str
        Path to point source file.
    
    return_flux: bool
        Return the point source flux as an array

    Returns:
    --------
    
    returns 4 or 5 arrays:
      -indices, ras, decs, source_radii [,flux]

    '''
    if not os.path.isfile(filename):
        filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), filename)
    if not os.path.isfile(filename):
        raise FileNotFoundError(filename)
    with open(filename) as f:
        indices = []
        ras  = []
        decs = []
        source_radii = []
        flux = []
        for line in f:
            ls = line.strip().split()
            if len(line.strip()) == 0 or line.strip()[0] == '#':
                continue
            if ls[0] == 'VALIDITYBEG' or ls[0] == 'VALIDITYEND':
                continue
            if len(ls) <= 3:
                continue
            indices.append( int( float(ls[0]) ))
            ras.append( core.G3Units.deg * float( ls[1] ))
            decs.append(core.G3Units.deg * float(ls[2]))
            source_radii.append(core.G3Units.deg * float(ls[3]))
            if return_flux:
                flux.append(float(ls[4]))
    if return_flux:
        return np.asarray(indices), np.asarray(ras), np.asarray(decs), np.asarray(source_radii), np.asarray(flux)
    else:
        return np.asarray(indices), np.asarray(ras), np.asarray(decs), np.asarray(source_radii)

def make_point_source_map(sky_map, point_source_file, mask_oob_pixels = True):
    '''
    Parses the point source file and fills a point source mask into the sky map.

    mask_oob_pixels: whether or not we mask timestream data outside the map
    '''

    indices, ras, decs, source_radii  = parse_point_source_file(point_source_file)
    mapmaker.make_point_source_mask_cpp(ras, decs, source_radii, mask_oob_pixels, sky_map)

    
    
def get_brightest_sources(filename, subfield=None, n=10, ra_lim=[-35,35], dec_lim=[-41,-70]):
    '''
    Reads in a point source file, filename, and returns the brightest n sources within the 
    user-defined region. If a subfield is specified, the dec_lims are set to that subfield.
    
    Arguments:
    -----------
    filename: str
        path to point source file.
        
    n: int
      The brightest n sources to return.
     
    subfield: str
        If you want to search a subfield, enter the subfield name here i.e. 'ra0hdec-44.75'
    
    ra_lim: list --> [float,float]
        The lower and upper RA extent of the field [deg]
    
    dec_lim: list --> [float,float]
        The lower and upper Dec extent of the field [deg]
    
    Returns:
    --------
        brightest_sources: dict
            A Dictionary of the n brightest sources in the region, keyed by source ID. 
            The dictionary contains RA, Dec and Flux for the source.
    
    '''
    ptsrc_info = parse_point_source_file(filename,return_flux=True)
    
    if subfield is not None:
        if subfield == 'ra0hdec-44.75':
            dec_lim[1] = -48.5
        elif subfield == 'ra0hdec-52.25':
            dec_lim[0] = -48.5
            dec_lim[1] = -56.
        elif subfield == 'ra0hdec-59.75':
            dec_lim[0] = -56.
            dec_lim[1] = -63.5
        else:
            dec_lim[0] = -63.5
    ID  = ptsrc_info[0]
    RA  = ptsrc_info[1]/core.G3Units.deg
    Dec = ptsrc_info[2]/core.G3Units.deg
    flux = ptsrc_info[4]
    
    ## Rotate ra so that RA range goes from - to + about 0 deg.
    RA[RA>180.] -= 360.
    ra_lim  = np.array(ra_lim)
    dec_lim = np.array(dec_lim)
    np.subtract(ra_lim[ra_lim>180.], 360.)
    
    ## since we're at the south pole... check this
    if dec_lim[0] > dec_lim[1]:
        deccut = ID[np.logical_and(Dec<dec_lim[0], Dec>dec_lim[1])]
    else:
        deccut = ID[np.logical_and(Dec>dec_lim[0], Dec<dec_lim[1])]
    racut = ID[np.logical_and(RA>ra_lim[0],RA<ra_lim[1])]
    sources_in_map = np.intersect1d(racut,deccut)
    
    brightest_sources = {}
    els = sources_in_map-1
    bright_els = np.argsort(flux[els])[-n:]
    for bs in bright_els:
        brightest_sources[ID[els][bs]] = {}
        brightest_sources[ID[els][bs]]['RA']   = RA[els][bs]
        brightest_sources[ID[els][bs]]['Dec']  = Dec[els][bs]
        brightest_sources[ID[els][bs]]['flux'] = flux[els][bs]
    
    return brightest_sources
    
    
def find_sources_in_map(map, pixel_mask = None,
                        default_mask = True,
                        nsigma=10, threshold_mjy=50,
                        frequency = 150*core.G3Units.GHz,
                        beamsize = 1.2*core.G3Units.arcmin,
                        omega_b = None,
                        plot_sources = False,
                        write_ptsrc_file = False,
                        filename='ptsrc_list.txt',
                        mask_radius = 5*core.G3Units.arcmin,
                        include_flux = False):
    '''
    Finds sources in a map detected at a specified significance
    above a given flux threshold. Optionally writes a point source
    file suitable for handing to the mapmaker.
    
    Parameters:
    -----------
    map: Either an unweighted temperature FlatSkyMap or a Map Frame
        containing a temperature FlatSkyMap.
        
    pixel_mask [None]: Mask applied to map before source finding.
    
    default_mask [True]: If True, a pixel mask will be calculated
        for you. Input map must be a Map frame with weights.
        
    nsigma [10]: Only find sources detected at this significance.
    
    threshold_mjy [50]: Only find sources above this flux threshold.
    
    frequency [150GHz]: Observation frequency of the input map.
        Used for flux calculation.
        
    beamsize [1.2']: FWHM of beam. Used for flux calculation.
    
    omega_b: Beam solid angle (steradians) if known.
        Used for flux calculation.
    
    plot_sources [False]: Set to True to plot the map with sources
        circled. Useful for checking quality of result.
        
    write_ptsrc_file [False]: Set to True to create a file of source
        positions, e.g. for handing to the mapmaker.
        
    filename: Where to save the point source list.
    
    mask_radius [5']: Exclusion radius around found sources.
    
    include_flux [False]: Set to True to include Flux column in 
        point source list. Mapmaker does not expect this column.
        
    Returns:
    --------
    Dictionary with source positions, detection significances, and fluxes.
        Sources labelled with sequential integers.
    '''
    wgt = None
    if isinstance(map, core.G3Frame):
        for k in map.keys():
            if 'W' in k and 'pol' in k:
                wgt = map[k]
        if map['T'].is_weighted:
            if wgt is None:
                raise KeyError("'Wunpol' or 'Wpol' not found in frame.")
            Tmap = map['T']/wgt.TT
        else:
            Tmap = map['T']
    elif isinstance(map, coordinateutils.FlatSkyMap):
        if map.is_weighted:
            raise TypeError('Input map must be unweighted')
        Tmap = map
    else:
        raise TypeError('Unsupported input type. '+
                        'Must be G3Frame or FlatSkyMap.')
        
    if pixel_mask is None:
        if default_mask == True:
            if wgt is None:
                raise TypeError('Input map must be a Map Frame with weights '+
                                   'to use default_mask')
            from spt3g.mapspectra import apodmask
            pixel_mask = apodmask.make_border_apodization(wgt, apod_type='cos',
                                                        radius_arcmin=15.)
        else:
            pixel_mask = np.ones(np.shape(Tmap))
            
    all_sources = find_sources_quick(Tmap/core.G3Units.uK, nsigma = nsigma,
                                     pixel_mask = pixel_mask,
                                     reso_arcmin=Tmap.res/core.G3Units.arcmin)
    out_sources = dict()
    i=0
    for src, info in all_sources.items():
        flux = uk_to_mjy(info['peakval'], frequency = frequency,
                         beamsize = beamsize, omega_b = omega_b)
        if flux < threshold_mjy:
            continue
        ra, dec = Tmap.pixel_to_angle(int(round(info['xpeak'])),
                                      int(round(info['ypeak'])))
        out_sources[i] = {'uk':info['peakval'], 'flux_mjy':flux,
                          'sn':info['peaksig'],
                          'xpeak':info['xpeak'], 'ypeak':info['ypeak'],
                          'ra':ra/core.G3Units.deg,
                          'dec':dec/core.G3Units.deg}
        i += 1
        
    if plot_sources:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,6))
        sdtemp = np.std(np.array(Tmap)[np.where(pixel_mask)])
        atemp = gaussfit_hist(np.array(Tmap)[np.where(pixel_mask)],
                              1000, -5.*sdtemp, 5.*sdtemp, do_plot = False)
        plt.imshow(Tmap*pixel_mask, origin='lower', cmap=plt.cm.gray,
                   vmin = -5.*atemp['params'][2],
                   vmax = 5.*atemp['params'][2])
        xpeaks = [info['xpeak'] for src, info in out_sources.items()]
        ypeaks = [info['ypeak'] for src, info in out_sources.items()]
        plt.plot(xpeaks, ypeaks,'ro', markersize=12, mfc='None')
            
    if write_ptsrc_file:
        f = open(filename,'w')
        header = "# Index\tRA\tDEC\tRadius"
        units="# \t(deg)\t(deg)\t(deg)"
        if include_flux:
            header += "\tFlux"
            units += "\t(mJy)"
        f.write(header+"\n")
        f.write(units+"\n")
        for src,info in out_sources.items():
            line = "%s\t%.5f\t%.5f\t%.4f"%(
                src, info['ra'],info['dec'],mask_radius/core.G3Units.deg)
            if include_flux:
                line += "\t%.4f"%info['flux_mjy']
            f.write(line+"\n")
        f.close()
    
    return out_sources
            
def find_sources_quick(map, maprms=None, mapmean=0.,
                        nsigma=5, minrad_arcmin=0.5,
                        sigma_thresh_for_minrad=0,
                        reso_arcmin=0.25, pixel_mask=None):
    '''
    Quick 'n' dirty source finding
    
    Parameters:
    -----------
    map: 2d-array representing an unweighted Temperature map.
    
    maprms [None]: The 1-sigma noise level in the map. If not provided, will
        be calculated.
        
    mapmean [0.]: The zero-point of the map.
    
    nsigma [5.]: Required signal-to-noise to detect a source.
    
    minrad_arcmin [0.5]: (float or array-like) The required separation between
        detected sources. If given as a list, provides different radii for
        different-sigma sources.
        
    sigma_thresh_for_minrad [0]: (float or array-like) The source detection
        strengths corresponding to different exclusion radii. Only used if
        more than one element in minrad_arcmin.
        
    reso_arcmin [0.25]: Resolution of map, in arcminutes.
    
    pixel_mask [None]: Optional mask applied to map before source finding.

    Returns:
    --------
    output_struct: (dict) Contains a dictionary for each source, labelled by
        sequential integers, with keys 'xpeak', 'ypeak','peakval','peaksig'.
    
    Jan 2019: DPD ported from spt_analysis/sources/find_sources_quick.pro
    '''
    
    map = np.asarray(map)
    if not isinstance(minrad_arcmin, list):
        minrad_arcmin = [minrad_arcmin]
    if not isinstance(sigma_thresh_for_minrad, list):
        sigma_thresh_for_minrad = [sigma_thresh_for_minrad]
        
    if len(sigma_thresh_for_minrad) != len(minrad_arcmin):
        raise ValueError( "If you are specifying multiple avoidance radii,"+
                         "please supply a threshold level for each one.")

    if pixel_mask is None:
        pixel_mask = np.ones(map.shape)

    # get rms in map if not supplied
    if maprms is None:
        whn0 = np.where(abs(map*pixel_mask) > 1.e-8)
        if len(whn0[0]) == 0:
            maprms = np.nanstd(map)
        else:
            temprms = np.std(map[whn0])
            atemp = gaussfit_hist(map[whn0],1000,-8.*temprms,8.*temprms,
                                          do_plot = False)
            maprms = atemp['params'][2]

    # simple source finder
    peaks = find_groups(map*pixel_mask, maprms, offset = mapmean,
                           nsigma=nsigma, minnum=1)
    npeaks = peaks['n_detected']
    if npeaks==0:
        print('No sources found')
        return

    # gather detected peaks into output structure, ignoring repeat
    # detections of same object
    xpeaks = peaks['xcen']
    ypeaks = peaks['ycen']
    peakvals = peaks['maxvals']
    peaksigs = (peakvals-mapmean)/maprms
    peak_assoc = np.zeros(npeaks)

    output_struct = dict()
    for i in np.arange(npeaks):
        output_struct[i] = {'xpeak':xpeaks[0],'ypeak':ypeaks[0],
                            'peakval':peakvals[0],'peaksig':peaksigs[0]}
        
    minrad_pix = minrad_arcmin[0]/reso_arcmin

    ksource = 1

    # different accounting if exclusion radius is specified as a function
    # of significance
    if len(minrad_arcmin) > 1:
        minrad_pix_all = np.zeros(npeaks)
        sthresh = sort(sigma_thresh_for_minrad)
        for j in np.arange(len(minrad_arcmin)):
            i = sthresh[j]
            whgthresh = np.where(peaksigs >= sigma_thresh_for_minrad[i])[0]
            if len(whgthresh) > 0:
                minrad_pix_all[whgthresh] = minrad_pix[i]
        minrad_pix_os = np.zeros(npeaks)
        minrad_pix_os[0] = minrad_pix_all[0]
        for j in np.arange(npeaks):
            prev_x = np.array([output_struct[n]['xpeak'] for n in \
                               np.arange(0,ksource)])
            prev_y = np.array([output_struct[n]['ypeak'] for n in \
                               np.arange(0,ksource)])
            distpix = np.sqrt((prev_x - xpeaks[j])**2+(prev_y - ypeaks[j])**2)
            whclose = np.where(distpix <= minrad_pix_os[0:ksource])[0]
            if len(whclose) == 0:
                output_struct[ksource] = {xpeak:xpeaks[j],ypeak:ypeaks[j],
                                          peakval:peakvals[j],
                                          peaksig:peaksigs[j]}
                peak_assoc[j] = ksource
                minrad_pix_os[ksource] = minrad_pix_all[j]
                ksource += 1
            else:
                mindist = min(distpix)
                peak_assoc[j] = distpix.argmin()
    else:
        for j in range(npeaks):
            prev_x = np.array([output_struct[n]['xpeak'] for n in \
                               np.arange(0,ksource)])
            prev_y = np.array([output_struct[n]['ypeak'] for n in \
                               np.arange(0,ksource)])
            distpix = np.sqrt((prev_x - xpeaks[j])**2+(prev_y - ypeaks[j])**2)
            mindist = min(distpix)
            if mindist > minrad_pix:
                output_struct[ksource] = {
                    'xpeak':xpeaks[j], 'ypeak':ypeaks[j],
                    'peakval':peakvals[j], 
                    'peaksig':(peakvals[j]-mapmean)/maprms}
                peak_assoc[j] = ksource
                ksource += 1
            else:
                peak_assoc[j] = distpix.argmin()
    for src in list(output_struct.keys()):
        if src>=ksource:
            del output_struct[src]
            
    return output_struct

def find_groups(Tmap, signoise = None, offset = 0,
                 minnum = 2.0, nsigma = 5.):
    '''
    Given a 2d array (a map), will find groups of elements (pixels) that
    are spatially associated (sources).

    Parameters:
    -----------
    Tmap: 2d array, e.g. representing an unweighted T map.
    
    signoise [None]: 1-sigma noise level in map. If not provided, is calculated
        from the map.
        
    offset [0]: zero point of map.
    
    minnum [2]: minimum number of pixels needed to form a group.
    
    nsigma [5.]: required detection threshold for a group.
        
    Output:
    -------
    Dictionary with following keys:
        xcen - array of x-location of group centers.
        ycen - array of Y-location of group centers.
        maxvals - array of heights (in map units) of found objects.
        sigvals - array of heights (in significance units) of found objects.
        n_detected - Number of detected sources.
       
    NOTES:
        Jan 2019: DPD ported this from
            spt_analysis/sources/find_groups.pro
            sptpol_software/scratch/tnatoli/transient/find_groups.py
    '''
    if not isinstance(Tmap, np.ndarray):
        Tmap = np.array(Tmap)
    assert len(Tmap.shape) == 2
        
    nxgrid = Tmap.shape[1]
    nygrid = Tmap.shape[0]
    
    # get rms in map if not supplied
    if signoise is None:
        whn0 = np.where(abs(Tmap) > 1.e-8)
        if len(whn0[0]) == 0:
            signoise = np.nanstd(Tmap)
        else:
            temprms = np.std(Tmap[whn0])
            atemp = gaussfit_hist(Tmap[whn0],1000,-8.*temprms,8.*temprms,
                                          do_plot = False)
            signoise = atemp['params'][2]
    # make a dictionary that found groups will be returned in
    #  also add dummy values to be returned if no peaks/groups found
    groups = dict()
    groups['maxvals'] = 0.
    groups['sigvals'] = 0.
    groups['xcen'] = 0.
    groups['ycen'] = 0.
    groups['n_detected'] = 0
    
    # finding points above detection_wanted_sigma*noise_level
    plist = np.where(np.asarray(Tmap)-offset >= nsigma*signoise)
    npix = len(plist[0])

    # if npix is less than minnum, then you have not found enough sources
    if npix < minnum:
        print('Not enough pixels above threshold')
        return groups
    
    xlist = plist[1] # X-coord of pix > thres
    ylist = plist[0] # Y-coord of pix > thres
    plist = list(zip(plist[0], plist[1])) # list of (y,x) coordinate pairs
    
    if npix == 1:
        gbounds = [0,0]
        groups['xcen'] = xlist
        groups['ycen'] = ylist
        groups['n_detected'] = 1
        groups['sigvals'] = (Tmap[ylist,xlist]-offset)/signoise
        groups['maxvals'] = Tmap[ylist,xlist]
        return groups
    
    # assiciate into groups
    gplist = []
    objsize = []
    npdone = 0 # number of pixels done

    # take first pixel above threshold and find other pixels above
    # threshold that touch first pixel (diagonal is OK).
    while (npdone < npix and len(plist) >= minnum):
        firstpix = plist[0]
        firstx = xlist[0]
        firsty = ylist[0]
        if len(xlist) > 1:
            plist = plist[1:]
            xlist = xlist[1:]
            ylist = ylist[1:]
        npdone += 1
        thisgplist = [firstpix]
        thisnp = 1
        # if another thresholded pixel is touching the first pixel,
        # include it in "who's associated"
        whassoc = np.where((abs(xlist-firstx) <= 1) &
                           (abs(ylist-firsty) <= 1))[0]
        countassoc = len(whassoc)
        # Also record "who's not associated"
        whnassoc = np.array([i for i in np.arange(len(xlist))\
                             if i not in whassoc])
        
        # now take each pixel above threshold that touched first pixel and
        # find pixels above threshold that touch them.  Repeat until we
        # find no new associated pixels.
        while countassoc != 0:
            templist = [plist[i] for i in whassoc]
            thisgplist += templist
            npdone += countassoc
            thisnp += countassoc
            if len(whnassoc)==0:
                # If every pixel is associated in one group, stop
                countassoc = 0
            else:
                plist = [plist[i] for i in whnassoc]
                xlist = np.array([coord[1] for coord in plist])
                ylist = np.array([coord[0] for coord in plist])
                whassoc = []
                catemp1 = 0
                for i in np.arange(countassoc):
                    thispix = templist[i]
                    thisx = thispix[1]
                    thisy = thispix[0]
                    whatemp = np.where((abs(xlist-thisx) <= 1) &
                                       (abs(ylist-thisy) <= 1))[0]
                    catemp2 = len(whatemp)
                    if catemp2 > 0:
                        whassoc += list(whatemp)
                        catemp1 += catemp2
                if catemp1 > 0:
                    whassoc = np.unique(whassoc)
                    countassoc = len(whassoc)
                    whnassoc = np.array(list(
                        set(whassoc)^set(np.arange(len(plist)))))
                else:
                    countassoc = 0
                    
        # If this number of pixels is big enough to be a group, add it  
        if thisnp >= minnum:
            gplist += thisgplist
            objsize.append(thisnp)
            
    # clean up lists
    if len(gplist) == 0:
        print('Whole map identified as one object.')
        print('Maybe your threshold is too low? Exiting...')
        return groups
    
    plist = gplist
    nobj = len(objsize)
    gbounds = np.zeros((nobj,2))
    gbounds[0] = [0, objsize[0]-1]
    for i in np.arange(1, nobj):
        gbounds[i,0] = gbounds[i-1,1]+1
        gbounds[i,1] = gbounds[i,0] + objsize[i]-1

    # find object centers and peak values
    xcen = np.zeros(nobj)
    ycen = np.zeros(nobj)
    maxvals = np.zeros(nobj)
    sigvals = np.zeros(nobj)
    
    for i in np.arange(0, nobj):
        pixtemp = plist[int(gbounds[i,0]):int(gbounds[i,1])+1]
        xtemp = [coord[1] for coord in pixtemp]
        ytemp = [coord[0] for coord in pixtemp]
        minx = min(xtemp)
        miny = min(ytemp)
        nxtemp = max(xtemp) - minx + 1
        nytemp = max(ytemp) - miny + 1
        objmtemp = np.zeros((nytemp,nxtemp))
        for kpix in np.arange(len(pixtemp)):
            objmtemp[ytemp[kpix]-miny, xtemp[kpix]-minx] = Tmap[pixtemp[kpix]]
        if nytemp > 1:
            xobj = np.sum(objmtemp, axis = 0)
        else:
            xobj = objmtemp
        yobj = np.sum(objmtemp, axis = 1)
        xcen[i] = np.sum(xobj*(np.arange(0,nxtemp)+1.))/np.sum(xobj) + minx -1
        ycen[i] = np.sum(yobj*(np.arange(0,nytemp)+1.))/np.sum(yobj) + miny -1
        maxvals[i] = objmtemp.max()
        sigvals[i] = objmtemp.max()/signoise
        
    rsmv = np.argsort(maxvals)[::-1]
    groups['maxvals'] = maxvals[rsmv]
    groups['sigvals'] = sigvals[rsmv]
    groups['xcen'] = xcen[rsmv]
    groups['ycen'] = ycen[rsmv]
    groups['n_detected'] = nobj
    
    return groups


def uk_to_mjy(input, invert=False, frequency = 150*core.G3Units.GHz,
              beamsize = 1.2*core.G3Units.arcmin, omega_b = None):
    '''
    Converts an intensity in uK_CMB to a flux in mJy and vice-versa.
    Taken from the IDL mjy_to_uk.pro
    
    mJy_per_uK obtained by taking derivative of Planck's Law with respect to
    temperature, evaluated at T_CMB and band center frequency, and multiplying
    by beam solid angle.

    Parameters:
    -----------
    input: Either a flux in mJy or a temperature in uK.
    
    invert [False]: (bool) If False, convert uK to mJy.
        If True, convert mJy to uK.
        
    frequency [150GHz]: Observing band center frequency, in G3Units.
    
    beamsize [1.2arcmin]: FWHM of the beam, in G3Units.
    
    omega_b [None]: Specify beam solid angle (steradians) if known.
    
    Returns:
    --------
    input quantity converted to mJy or uK
    '''
    frequency /= core.G3Units.Hz
    beamsize /= core.G3Units.rad

    tcmb = 2.725
    c = 2.9979e8
    h = 6.62607e-34
    kb = 1.38065e-23
    
    # x contains all frequency dependence
    x = (h/kb)*(frequency/tcmb)
    # i0 is prefactor that contains everything but freq dependence and beam
    i0 = 2.* kb * (kb*tcmb/(h*c))**2

    if omega_b is None:
        # gaussian beam solid angle where beamsize = FWHM
        omega_b = (np.pi/(4*np.log(2))) * beamsize**2

    # Numbers out front are to convert from SI to milliJanskys and uK to K
    mjy_per_uK = (1e29/1e6)*(i0*omega_b*np.exp(x)* x**4)/((np.exp(x)-1.)**2)

    if invert: 
        return input/mjy_per_uK
    else:
        return mjy_per_uK*input

def mjy_to_uk(input, **kwargs):
    return uk_to_mjy(input, invert=True, **kwargs)


def ptsrc_uk_to_mjy(year = 2018, band = '150', alpha = -0.5,
                    beamsize = 1.2*core.G3Units.arcmin, omega_b=None,
                    verbose = False):
    '''
    Calculates the conversion factor between uK and mJy for a specified source
    spectrum power law and the measured SPT bandpasses and beam.

    Parameters:
    ----------
    year [2018]: (int)
        The version of the focal plane used to take data

    band ['150']: (str)
        The detector observing band to use for the FTS spectra.
        The center frequency of the bandpass will be calculated based on
        the FTS spectra for the respective band and year.

    alpha [-0.5]: (float)
        The power for the power law of the population spectra.
        -0.5 = radio (most sources in the field)
        2 = Rayleigh-Jeans
        3.5 = Dusty

    beamsize [1.2arcmin]: (float)
        Beam FWHM in G3Units.

    omega_b [None]: (float)
        The solid angle of the detector in steradians.
        If None, a gaussian approximation of the beam will be integrated
        over for the solid angle of the detector. If using a matched source
        (optimal) filter, the filter creation code will spit out 'area_eff'
        which is the source profile (times the beam and filter tranfer function)
        times the optimal filter integrated over is k_x and k_y.

    Returns:
    --------
    mJy_per_uK
    '''
    # constants
    kb = 1.38065e-23 # m^2 kg /(s^2 K)
    c = 2.9979e8  # m/s
    h = 6.62607e-34 # m^2 kg/s
    tcmb = 2.725 # K

    # check inputs
    beamsize /= core.G3Units.rad
    if not isinstance(year, int):
        year = int(year)
    if isinstance(band, int):
        band = str(band)
    bands = ['90','150','220']
    if band not in bands:
        raise ValueError('Band must be one of '+str(bands))

    # choose the right bandpass file
    if year < 2012:
        raise NotImplementedError(year)
    if year < 2017:
        if band == '150':
            if year == 2012:
                spec_file = '/spt/user/ddutcher/fts/sptpol_150GHz_2012season.txt'
            if year >= 2013:
                spec_file = '/spt/user/ddutcher/fts/sptpol_150GHz_2013season.txt'
        elif band == '90':
            if year == 2012:
                spec_file = '/spt/user/ddutcher/fts/sptpol_90GHz_2012season.txt'
            if year >= 2013:
                spec_file = '/spt/user/ddutcher/fts/sptpol_90GHz_2013season.txt'
        else:
            raise ValueError('%s GHz did not exist in year %s'%(band, year))
        # read in the appropriate bandpass file
        dat = np.loadtxt(spec_file)
        freqs, amps = dat[:,0], dat[:,1]
    elif year == 2017:
        fts = files.load_pickle(
            '/spt/user/ddutcher/fts/2017_corrected_coadded_fts_spectra.pkl')
        freqs = fts[band]['freq']
        amps = fts[band]['spectrum']
    elif year == 2018:
        fts = files.load_pickle(
            '/spt/user/ddutcher/fts/2018_corrected_coadded_fts_spectra.pkl')
        freqs = fts[band]['freq']
        amps = fts[band]['spectrum']
    else:
        raise NotImplementedError(year)

    if freqs[0] < 1e9:
        # freqs was probably supplied in GHz already, so make it in Hz now
        freqs *= 1e9
    dnu = freqs[1:]-freqs[:-1]
    if len(set(dnu)) == 1:
        dnu = dnu[0]
        # Calculate the band center
    else:
        # this is used if freqs are not evenly spaced
        amps = amps[:-1]
        freqs = freqs[:-1]
    nu0 = np.sum(freqs*amps*dnu)/np.sum(amps*dnu)
    if verbose:
        print('band center = ',nu0/1e9, 'GHz')

    # there are no A*Omega (entendue) factors below because
    # that factor is already in the FTS spectrum data.
    def x(nu):
        # term in exponent in Planck's law
        return (h*nu/(kb*tcmb))
    def source_spec(nu,nu0,alpha):
        return (nu/nu0)**alpha
    def dBdT(nu):
        # derivative of Planck's law w.r.t temperature, evaluated at Tcmb
        return (1e29 * (2*kb*(kb*tcmb/(h*c))**2 * np.exp(x(nu)) * x(nu)**4) /
                ((np.exp(x(nu))-1.)**2))

    if omega_b is None:
        # gaussian beam solid angle where beamsize = FWHM
        omega_b = (np.pi/(4*np.log(2))) * beamsize**2
       # For a tophat beam: omega_b = np.pi*((FWHM/2)**2)

    mjy_per_uk = ( 1e-6 * np.sum(dBdT(freqs)*amps*dnu) /
                  np.sum(source_spec(freqs,nu0,alpha)*amps*dnu) )
    if verbose:
        print('omega_b = ',omega_b)
        print('mJy/uK/sr = ', mjy_per_uk)
        print('For point source, mJy/uK= ',mjy_per_uk*omega_b)

    return mjy_per_uk * omega_b

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
