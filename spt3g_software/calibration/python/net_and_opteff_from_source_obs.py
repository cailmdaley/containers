from spt3g import core
import numpy as np
from scipy import ndimage

#
# hard-coded for RCW38 for now
# 

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

    from spt3g.util.genericutils import shift
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
def net_and_opteff_wrapper(obsid):

    """
    wrapper for all the bits of the map-based NET calculation. returns
    frame to be written to file by a script in calibration/scripts
    """

    obsid_str = str(np.int(obsid))
    outdict = get_maps_and_call_net_calc(obsid)

# initialize our output calibration frame
    frame = core.G3Frame(core.G3FrameType.Calibration)
    data_names = outdict.keys()
    for s in data_names:
        frame[s] = core.G3MapDouble()

# stuff frame and write to file
    for k in outdict['Band'].keys():
        for s in data_names:
            frame[s][k] = outdict[s][k]

    return frame

                
#################################
def grab_noisedata(obsid, names, nframes=20):

    obsid_str = str(np.int(obsid))
    nbolo = len(names)

# get noise data & convert to watts
    try:
        ftemp = core.G3File('/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/0000.g3')
    except:
        try:
            ftemp = core.G3File('/spt/data/bolodata/downsampled/RCW38/'+obsid_str+'/0000.g3')
        except:
            core.log_warn('Cannot find obsid '+obsid_str+' in /spt/data/bolodata/downsampled/RCW38-pixelraster or /spt/data/bolodata/downsampled/RCW38-pixelraster.\n')
            
    noisedict1 = {}
    from spt3g.dfmux.unittransforms import ConvertTimestreamUnits
    for i in np.arange(nbolo):
        noisedict1[names[i]] = []
        converter = ConvertTimestreamUnits(Input='RawTimestreams_I')
    i = 0
    while i < nframes:
        frame = ftemp.next()
        if frame.type is core.G3FrameType.Wiring:
            converter(frame)
        if frame.type is core.G3FrameType.Scan:
            converter(frame)
            for name in names:
                try:
                    noisedict1[name].append(frame['CalTimestreams'][name])
                except:
                    pass
        i += 1

    noisedict = {}
    for name in names:
        ntemp = noisedict1[name]
        npts_tot = 0
        for i in np.arange(len(ntemp)):
            npts_tot += len(ntemp[i])
        noisedict[name] = np.zeros(npts_tot)
        npts_running = 0
        for i in np.arange(len(ntemp)):
            thisnpts = len(ntemp[i])
            noisedict[name][npts_running:npts_running+thisnpts] = ntemp[i]
            npts_running += thisnpts
        noisedict[name] = noisedict[name][npts_running/4:3*npts_running/4]

    return noisedict

#################################
def get_maps_and_call_net_calc(obsid):

    obsid_str = str(np.int(obsid))

# set up output dicts
    psd_dict = {}
    banddict = {}
    netdict = {}
    opteffdict = {}
    wndict = {}
    peakdict = {}
    rmsdict = {}
    mapdict = {}

# get bolo properties file
    try:
        ftemp2 = core.G3File('/spt/data/bolodata/downsampled/RCW38-pixelraster/'+obsid_str+'/offline_calibration.g3')
    except:
        try:
            ftemp2 = core.G3File('/spt/data/bolodata/downsampled/RCW38/'+obsid_str+'/offline_calibration.g3')
        except:
            core.log_warn('Cannot find obsid '+obsid_str+' in /spt/data/bolodata/downsampled/RCW38-pixelraster or /spt/data/bolodata/downsampled/RCW38-pixelraster.\n')
    for frame in ftemp2:
        if frame.type is core.G3FrameType.Calibration:
            bp = frame['BolometerProperties']
        
# get map file and extract maps
    obsid_str = str(np.int(obsid))
    try:
        ftemp = core.G3File('/spt/user/production/calibration/RCW38-pixelraster/maps/'+obsid_str+'.g3')
    except:
        try:
            ftemp = core.G3File('/spt/user/production/calibration/RCW38/maps/'+obsid_str+'.g3')
        except:
            core.log_warn('Cannot find obsid '+obsid_str+' in /spt/data/bolodata/downsampled/RCW38-pixelraster or /spt/data/bolodata/downsampled/RCW38-pixelraster.\n')

    for frame in ftemp:
        if frame.type is core.G3FrameType.Wiring:
            names = frame['WiringMap'].keys()
        if frame.type is core.G3FrameType.Map:
            name = frame['Id']
            if name != 'bsmap':
                if 'Wunpol' in frame:
                    twt = np.asarray(frame['Wunpol'].TT)
                    twt_orig = twt.copy()
                    twt[np.where(twt == 0.)] = 1.
                if name in names:
                    map_weighted = np.asarray(frame['T'])
                    reso_arcmin = (frame['T']).res/core.G3Units.arcmin
                    if np.max(np.abs(map_weighted)) > 0.:
                        map_unw = map_weighted.copy()
                        map_unw /= twt
                        mapdict[name] = map_unw
                        banddict[name] = np.int(bp[name].band/10.)

# get time-ordered data for noise
    noisedict = grab_noisedata(obsid, names)

# W/K for a perfect single-pol bolo with assumed bandwidths
    kb = 1.38e-23
    wpk_perfect = {}
    wpk_perfect[90] = kb*28.2e9 # new from Brad
    wpk_perfect[150] = kb*42.6e9
    wpk_perfect[220] = kb*51.9e9

# rough T_RJ to T_CMB conversion
    rj2cmb = {}
    rj2cmb[90] = 1.23
    rj2cmb[150] = 1.73
    rj2cmb[220] = 3.07

    for mname in mapdict.keys():
        map_unw = mapdict[mname]
        try:
            net, w_per_k, whitelevel = \
                quick_net_from_rcw38(map_unw,noisedict[mname], 
                                             invert_map=True,band=banddict[mname], 
                                             rate=76.3,frange=[21.,24.])
            psddicttemp = quick_pspec(noisedict[mname],rate=76.3)
            psd_dict[mname] = psddicttemp['psd']
            freqs = psddicttemp['freqs']
        except:
            net = 1e6
            w_per_k = 0.
            whitelevel = 1e6
        netdict[mname] = net*rj2cmb[banddict[mname]]
        opteffdict[mname] = w_per_k/wpk_perfect[banddict[mname]]
        wndict[mname] = whitelevel
        amap_sm = ndimage.gaussian_filter(-map_unw,2./reso_arcmin)
        ny,nx = np.shape(map_unw)
        ycenter, xcenter = np.unravel_index(np.argmax(amap_sm),[ny,nx])
        peakdict[mname] = map_unw[ycenter,xcenter]
        try:
            rmsdict[mname] = np.std(map_unw[ycenter+20:ycenter+60,xcenter-60:xcenter-20])
        except:
            rmsdict[mname] = 1e6

# output dict of dicts
    outdict = {}
    outdict['Band'] = banddict
    outdict['NETCMB'] = netdict
    outdict['OptEff'] = opteffdict
    outdict['WhiteLevel'] = wndict
    outdict['MapPeak'] = peakdict
    outdict['MapRMS'] = rmsdict

    return outdict
