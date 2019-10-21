from spt3g import core, dfmux, todfilter
from spt3g.timestreamflagging import add_flag
from spt3g.calibration.template_groups import get_template_groups

import numpy as np
from scipy import optimize
#todfilter.util.CutTimestreamsWithoutProperties dafuq

@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def FlagToBand(frame, band, flag_key = 'Flags', flag_reason = 'WrongBand',
               bolo_props=None):
    '''
    Flags any detectors with a band specified in the bolometer properties different from
    band
    '''

    assert(not bolo_props is None)
    bad_bolos = []
    for k in bolo_props.keys():
        if bolo_props[k].band != band:
            bad_bolos.append(k)
    add_flag(frame, flag_key, flag_reason, bad_bolos)


@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def FlagIncompletePixelPairs(frame, flag_key = 'Flags',
                             ts_key = 'RawTimestreams_I',
                             flag_reason = 'IncompletePixelPair',
                             bolo_props = None):
    '''
    For every polarization pair of detectors, checks whether there are two
    detectors present in the BolometerProperties, neither of which is
    permitted to be in the list of flagged bolos. Otherwise, flags the pair.
    This enforces having both polarizations for a band in a pixel live.
    '''
    if frame.type != core.G3FrameType.Scan:
        return

    if flag_key in frame.keys():
        flags = frame[flag_key]
    else:
        flags = []
    pixel_groups = get_template_groups(bolo_props, per_band = True,
                                       per_pixel = True, per_wafer = True,
                                       include_keys=True)

    bad_bolos = []
    for group, bolos_in_pixel in pixel_groups.items():
        good_bolos_in_pixel = [bolo for bolo in bolos_in_pixel \
                               if bolo not in flags and \
                               bolo in frame[ts_key]]
        bolos_in_pixel_in_frame = [bolo for bolo in bolos_in_pixel \
                                   if bolo in frame[ts_key]]
        if len(good_bolos_in_pixel) != 2:
            bad_bolos += bolos_in_pixel_in_frame

    add_flag(frame, flag_key, flag_reason, bad_bolos)


@core.indexmod
def FlagHighQ(frame, i_data, q_data, variance_threshold=5, flag_key='Flags',
              flag_id='HighQ'):
    '''
    Flag detectors in i_data for which q_data has a total variance
    more than variance_threshold times higher. Because this is
    total variance, prefiltering the data is a good idea.
    '''
    if i_data not in frame:
        return

    bad_bolos=[]
    q = frame[q_data]
    for k,i_ts in frame[i_data].iteritems():
        if np.var(q[k]) > np.var(i_ts)*variance_threshold:
            bad_bolos.append(k)
    add_flag(frame, flag_key, flag_id, bad_bolos)

@core.pipesegment
def FlagHighQWithPolyFilter(pipe, i_data='RawTimestreams_I',
                            q_data='RawTimestreams_Q', variance_threshold=5,
                            poly_order = 1, flag_key='Flags', flag_id='HighQ'):
    '''
    Run FlagHighQ after polyfiltering timestreams at the order given.
    '''
    pipe.Add(todfilter.MaskedPolyHpf, in_ts_map_key=i_data,
             out_ts_map_key='__Detrended' + i_data, poly_order=poly_order)
    pipe.Add(todfilter.MaskedPolyHpf, in_ts_map_key=q_data,
             out_ts_map_key='__Detrended' + q_data, poly_order=poly_order)
    pipe.Add(FlagHighQ, i_data='__Detrended' + i_data,
             q_data='__Detrended' + q_data,
             variance_threshold=variance_threshold,
             flag_key=flag_key, flag_id=flag_id)
    pipe.Add(core.Delete, keys=['__Detrended' + i_data, '__Detrended' + q_data])
  

@core.indexmod
def FlagMissing(frame, key_ts_key, flagged_ts_key, 
                flag_key = 'Flags', flag_id = 'MissingInOther'):
    '''
    Flags any keys missing from frame[key_ts_key] that are present in frame[flagged_ts_key] 
    '''

    if frame.type != core.G3FrameType.Scan:
        return    
    bad_bids = []
    oks = frame[key_ts_key].keys()
    for k in frame[flagged_ts_key].keys():
        if not k in oks:
            bad_bids.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_bids)

@core.indexmod
def FlagNotPresentInAny(frame, key_ts_keys, flagged_ts_key,
                        flag_key = 'flags', flag_id = 'MissingInOthers'):
    '''
    Flags any keys not present in any of [frame[k] for k in key_ts_keys] that
    are present in frame[flagged_ts_key].  Any of key_ts_keys may or may not
    be present in the frame.
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    # Removing any keys found in any of the key_ts_key map objects
    # from the set of keys in the flagged_ts_key map object
    bad_ks = set(frame[flagged_ts_key].keys())
    for key_ts_key in key_ts_keys:
        if key_ts_key in frame:
            bad_ks -= set(frame[key_ts_key].keys())
    # remaining keys are to be flagged
    bad_ks = sorted(bad_ks)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_ks)

@core.indexmod
def FlagNaNs(frame, ts_key, flag_key = 'Flags', flag_id = 'HasNaNs'):
    '''
    Flags any timestreams with NaNs in timestream.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    bad_ks = []
    tsm = frame[ts_key]
    for k in tsm.keys():
        sm = np.sum(tsm[k])
        if ( not np.isfinite(sm) or sm == 0):
            bad_ks.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_ks)

@core.scan_func_cache_data(wiring_map = 'WiringMap')
def FlagNegativeDANChannels(frame, ts_key, flag_key = 'Flags', flag_id = 'NegativeDAN', wiring_map=None):
    '''
    Flags any timestreams from DfMux channels with DAN on and negative values.
    These can't happen unless (a) the channels is secretly unbiased, (b) DAN
    phasing is horribly wrong, or (c) there is RF noise larger than the carrier.
    In any of those cases, the data are garbage.
    '''
    if frame.type != core.G3FrameType.Scan:
        return
    bad_ks = []
    hk = frame['DfMuxHousekeeping']
    tsm = frame[ts_key]
    for k,ts in tsm.items():
        h = dfmux.HousekeepingForBolo(hk, wiring_map, k)
        if not h.dan_streaming_enable:
            continue
        assert(ts.units == core.G3TimestreamUnits.Counts)
        if ( (ts < 0).any() ):
            bad_ks.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_ks)
    
@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def FlagTimestreamsWithoutProperties(frame, ts_key, flag_key = 'Flags', flag_id = 'NoProperties',
                                    bolo_props = None, flag_nan_offsets = False):
    '''
    Flags any timestreams without a key in bolo_props
    '''
    assert(not bolo_props is None)
    tsm = frame[ts_key]
    bad_ks = []
    nans = []
    for k in tsm.keys():
        if not k in bolo_props:
            bad_ks.append(k)
        else:
            if flag_nan_offsets and (np.isnan(bolo_props[k].x_offset) or np.isnan(bolo_props[k].y_offset)):
                nans.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_ks)
    if flag_nan_offsets:
        add_flag(frame, flag_key, 'NanOffsets', timestreams_to_flag = nans)
    
@core.scan_func_cache_data(RCW38_cal = 'RCW38FluxCalibration',
                          MAT5A_cal = 'MAT5AFluxCalibration')
def FlagMissingFluxCalibration(frame, ts_key, flag_key = 'Flags', flag_id = 'MissingFluxCalibration', 
                               RCW38_cal = None, MAT5A_cal = None):
    if RCW38_cal is not None:
        fluxcal = RCW38_cal
    elif MAT5A_cal is not None:
        fluxcal = MAT5A_cal
    else:
        raise AssertionError('No flux calibration')
    tsm = frame[ts_key]
    bad_ks = []
    for k in tsm.keys():
        if not k in fluxcal:
            bad_ks.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_ks)    

@core.scan_func_cache_data(bolo_props = 'BolometerProperties', wiring_map = 'WiringMap',
                           hk = 'DfMuxHousekeeping')
def FlagBadHousekeeping(frame, ts_key, flag_key = 'Flags', flag_id ='BadHk',
                        bolo_props = None, wiring_map = None, hk = None ):

    bad_bolos = []
    for k in frame[ts_key].keys():
        if not k in wiring_map:
            bad_bolos.append(k)
            continue
        a = dfmux.HousekeepingForBolo(hk,wiring_map,k, all_hk = True)
        mod_hk = a[2]
        chan_hk = a[3]
        if (chan_hk.dan_railed or mod_hk.demod_railed or 
            mod_hk.carrier_railed or mod_hk.nuller_railed or
            chan_hk.carrier_amplitude == 0 or
            chan_hk.carrier_frequency == 0 or chan_hk.demod_frequency == 0):            
            bad_bolos.append(k)
    add_flag(frame, flag_key, flag_id, timestreams_to_flag = bad_bolos)

@core.indexmod
class FlagBadRfrac(object):
    '''
    Flags bolometers with R/Rnorm < latch_co, R/Rnorm > overbias_co,
    or abs(R/Rnorm - Rcal/Rnorm) > delta_rfrac_co
    '''

    def __init__(self, ts_key, latch_co=0.5, overbias_co=0.98, delta_rfrac_co=0.03, flag_key='Flags',flag_on_delta_rfrac='False'):

        self.ts_key = ts_key
        self.latch_co = latch_co
        self.overbias_co = overbias_co
        self.delta_rfrac_co = delta_rfrac_co
        self.flag_key = flag_key
        self.cal_rfrac = None
        self.wiring_map = None
        self.hk = None
        self.rts_key = '_ResTimestreams'
        self.flag_on_delta_rfrac = flag_on_delta_rfrac
        self.convert = dfmux.ConvertTimestreamUnits(Input=self.ts_key, Output=self.rts_key,
                                                    Units=core.G3TimestreamUnits.Resistance)

    def __call__(self, frame):
        if 'CalibratorResponseRfrac' in frame:
            self.cal_rfrac = frame['CalibratorResponseRfrac']
        if 'WiringMap' in frame:
            self.wiring_map = frame['WiringMap']
        if 'DfMuxHousekeeping' in frame:
            self.hk = frame['DfMuxHousekeeping']

        self.convert(frame)

        if frame.type == core.G3FrameType.Scan:
            tsm = frame[self.ts_key]
            rtsm = frame[self.rts_key]
            latched_bolos = []
            overbiased_bolos = []
            changing_responsivity_bolos = []

            for k in tsm.keys():
                mean_ts = todfilter.get_mean(tsm[k]) # np.mean is slowwwwwww
                if mean_ts == 0:
                    latched_bolos.append(k)
                    continue
                chan_hk = dfmux.HousekeepingForBolo(self.hk, self.wiring_map, k)
                rnorm = chan_hk.rnormal
                if not rnorm:
                    continue
                res = todfilter.get_mean(frame[self.rts_key][k])
                if res/rnorm < self.latch_co:
                    latched_bolos.append(k)
                if res/rnorm > self.overbias_co:
                    overbiased_bolos.append(k)
                self.delta_rfrac_co = float(self.delta_rfrac_co)
                if self.cal_rfrac is not None:
                    if np.abs(res/rnorm - self.cal_rfrac[k]) > self.delta_rfrac_co:
                        changing_responsivity_bolos.append(k)
            add_flag(frame, self.flag_key, 'Latched', timestreams_to_flag=latched_bolos)
            add_flag(frame, self.flag_key, 'Overbiased', timestreams_to_flag=overbiased_bolos)
            if self.flag_on_delta_rfrac and self.cal_rfrac is not None:
                add_flag(frame, self.flag_key, 'ChangingRfrac',
                         timestreams_to_flag=changing_responsivity_bolos)

            core.Delete(frame, keys=[self.rts_key], type=core.G3FrameType.Scan)


@core.pipesegment
def FlagSaturatedBolos(pipe, ts_key, thresh=1.15, flag_key='Flags'):
    """
    Flags bolometers whose peak Mars scan is saturated. Determined by
    fitting a gaussian to 1) the first half of the Mars bump, 2) the second
    half of the Mars bumped, and 3) the full Mars bump. If the FWHM
    of 3 is at least thresh*mean(fwhm_1, fwhm2) bigger, flag the scan.
    """
    def FlagSaturatedBolosInner(frame):
        if frame.type != core.G3FrameType.Scan:
            return

        def gauss(x, *p):
            y = p[0]*np.exp(-np.power((x - p[1]), 2.)/(2. * p[2]**2.)) + p[3]
            return y

        def sig_to_arcmin(sig):
            am = 60*sig*2.*np.sqrt(2.*np.log(2.))*180/np.pi
            return am

        tsm = frame[ts_key]
        sat_bolos = []

        # the values below need not be 100% accurate-- ballpark to get close
        # to right starting gaussian fit params
        el_mars = 22.75 *np.pi / 180. 
        scanspeed_onsky = 0.3*np.cos(el_mars)
        arcmin_per_sample = 60. * scanspeed_onsky / 152.6
        nx = 300 # samples used for the gaussian fit
        xtemp = (np.arange(nx)-nx/2)*arcmin_per_sample

        for k in tsm.keys():
            d = -np.asarray(tsm[k]) # minus sign to make mars bump positive

            # first, subtract linear drift
            y0 = np.mean(d[0:50])
            y1 = np.mean(d[-50:])
            slope = np.arange(len(d)) / np.float(len(d))*(y1-y0) + y0
            d -= slope

            # check for high s/n peak not too close to scan edges
            imax = np.argmax(np.abs(d))
            sd = np.std(d[0:50])
            if np.abs(d[imax])/sd > 100. and imax > nx/2 and imax < len(d)-nx/2:
                # Mars detected!
                d2 = d[imax-nx//2:imax+nx//2] # data over which to fit gaussians
                p0 = [d[imax], 0.0, 0.6, 0.0] # starting parameter guesses
                try:
                    # Left half
                    errs0 = np.tile(sd, nx) # uniform weights base
                    errs0[nx//2+1:] *= 1000. # kill weight on right half
                    popt1, pcov1 = optimize.curve_fit(gauss, xtemp, d2, 
                                                      p0=p0, sigma=errs0)
                    y_fit_1 = gauss(xtemp, *popt1)

                    # Right half
                    errs0 = np.tile(sd, nx) # uniform weights base
                    errs0[:nx//2] *= 1000. # kill weight on left half
                    popt2, pcov2 = optimize.curve_fit(gauss, xtemp, d2, 
                                                      p0=p0, sigma=errs0)
                    y_fit_2 = gauss(xtemp, *popt2)

                    # All the data
                    errs0 = np.tile(sd, nx) # uniform weights
                    popt3, pcov3 = optimize.curve_fit(gauss, xtemp, d2, 
                                                      p0=p0, sigma=errs0)
                    y_fit_3 = gauss(xtemp, *popt3)
                    lr_mean = np.mean([popt1[2], popt2[2]])
                    if np.any(np.isinf([pcov1,pcov2,pcov3])):
                        continue
                    if popt3[2] / lr_mean >= thresh:
                        # call it saturated
                        sat_bolos.append(k)
                except(RuntimeError):
                    pass
        add_flag(frame, flag_key, 'Saturated', timestreams_to_flag=sat_bolos)
    pipe.Add(FlagSaturatedBolosInner)
