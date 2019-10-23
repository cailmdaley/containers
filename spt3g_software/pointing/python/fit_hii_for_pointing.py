# DEPRECATED
#
# Pointing corrections from HII region very_fast_point observations
# now calculated along with the rest of the VFP analysis in
# opacity.py.

import numpy as np
import pickle as pk
import pylab as py
import sptpol_software.util.time as time
import sptpol_software.util.files as files
from time import time as comptime
from sptpol_software.observation.receiver import Receiver
import sptpol_software.analysis.offline_pointing as op
import sptpol_software.analysis.pointing.pointing_tools as pt
import sptpol_software.util.mpfit as mpfit

def fit_hii_for_pointing(summary_filename, savename=True, datestart=True, datestop=True, quick=False, cut_c4=True, \
                             dodebug=True, do_plot=False, return_reco=False, test_mode=False, dosave=True, \
                             dowrite=True, writename=True, skip_cal=True, skip_elnod=True, skip_ampVsNoise=True, \
                             skip_response=False, use_this_source='rcw38', use_this_a5=None, summary_pkl=False, \
                             skip_position_fit=False, resid_limit=20., a4_limit=0.65, \
                             response_cut_low_long=1.25e-13, response_cut_high_long=1.7e-13, \
                             response_cut_low_short=1.25e-13, response_cut_high_short=1.7e-13, \
                             aztilt_config_file='sptpol_offline_tilts', use_this_a6=None, fit_a5=False, \
                             focus_position_input=None, filename_flag='', \
                             alter_HII_positions=True, a4_offset=0.0, a6_offset=0.0):
    """
    The purpose of this program is to fit a telescope imperfection model to obseravations of 
    HII regions (RCW38 and Mat5a).

    INPUTS:
        None

        a4_offset, a6_offset: The measured global offsets from point sources. (These are SUBTRACTED
                              from the best-fit a4 and a6 below.
    
    OUTPUTS:
        A txt file or pickle file can be saved, if requested by the user, which contains the 
        parameters for the telescope imperfection model fit to each HII observation.
    """
    #Bump the response cuts up if you're looking at mat5a.
    if use_this_source =='mat5a':
        response_cut_low_long += 0.25e-13
        response_cut_high_long += 0.25e-13

    if use_this_source == 'cena':
        response_cut_low_long = 0.3e-13
        response_cut_high_long = 0.6e-13

    #Check if certain keyword flags are turned on.
    if savename==True: savename = 'fit_hii_'+use_this_source+'.pkl'
    if datestart==True: datestart = '01-Jan-2011:00:00:00'
    if datestop==True: datestop = '01-Jan-2097:00:00:00'

    mjdstart = time.SptDatetime(datestart).mjd
    mjdstop = time.SptDatetime(datestop).mjd

    #If the HII observation happened after the last az tilt measurement, just skip it.
    tilt_date, meas_a2, meas_a3, tilt_ang = np.array(files.readConfig(aztilt_config_file, has_header=False))
    last_aztilt_mjd = np.max(tilt_date)
    mjdstop = np.min([mjdstop,last_aztilt_mjd])

    #Read in the summary pickle file with all of the source scan observations' information.
    with open(summary_filename, 'rb') as filein:
        s = np.array(pk.load(filein))

    #Do a bit of housekeeping on the structure.
    wh_mjd = np.nonzero(np.array([item['starttime_mjd'] for item in s]) != 0.0)[0]
    wh_nscans = np.nonzero(np.array([item['nscans']for item in s]) > 15.)[0]
    wh_goodsource1 = np.nonzero(np.array([item['source'] for item in s]) =='rcw38')[0]
    wh_goodsource2 = np.nonzero(np.array([item['source'] for item in s]) =='mat5a')[0]
    wh_goodsource3 = np.nonzero(np.array([item['source'] for item in s]) =='cena')[0]
    wh_goodsource_all1 = np.union1d(wh_goodsource1, wh_goodsource2)
    wh_goodsource = np.union1d(wh_goodsource_all1, wh_goodsource3)

    ##Right now nscans for all rcw38 observations are 0, so let's skip this cut. -- Fixed now.
    wh_ok1 = np.intersect1d(wh_mjd, wh_nscans)
    wh_ok2 = np.intersect1d(wh_ok1, wh_goodsource)
    #wh_ok2 = np.intersect1d(wh_mjd, wh_goodsource)
    s = s[wh_ok2]

    #Now sort the good s list by starttime mjd.
    slength = len(s)
    sorted_s_indices = sorted(range(slength), key=lambda k: ([item['starttime_mjd'] for item in s][k]))
    s = s[sorted_s_indices]

    n = len([item['source'] for item in s])

    #Create a dictionary to hold our results
    t = []
    for i in range(n):
        t.append({'starttime':'', 'starttime_mjd':0., 'mean_mjd':0.,
             'xoff':np.zeros(1599), 'yoff':np.zeros(1599), 'gb':np.zeros(1599, dtype=np.int16), 
             'xoff_err':np.zeros(1599), 'yoff_err':np.zeros(1599), 
             'focus_position':np.zeros(6),
             'a0':0., 'a1':0., 'a4':0., 'a5':0., 'a6':0., 
             'a0_err':0., 'a1_err':0., 'a4_err':0., 'a5_err':0., 'a6_err':0.,
             'a2':0., 'a3':0., 'az_tilt_mag':0., 'az_tilt_deg':0., 
             'chisq':0., 'dof':np.NaN,
             'temp_avg':0., 'wind_speed_avg':0., 'wind_dir_avg':0., 'pressure_avg':0.,
             'stddev_resx':0., 'stddeb_resy':0.,
             'med_resr':0., 'mean_resr':0.,
             'max_resx':0., 'max_resy':0., 'max_resr':0., 'mean_az':0.,
             'record_daz':0., 'record_del':0., 'record_det':0.,
             'use_daz':0., 'use_del':0., 'use_det':0.,
             'ngb':0., 'nscans':0.,
             'meanx':0., 'meany':0.,
             'daz':0., 'del':0., 'det':0.,
             'th_det':0., 'th_del':0., 
             'med_refract':0., 'stddev_refract':0.,
             'source':''})

    #Now loop over each source observation.
    print('...Looping over the HII observations...')

    for i in range(0,n):
    #for i in range(3,4):
        if s[i]['source'] != use_this_source: 
            #print('')
            #print('')
            #print('Source: ', s[i]['source'], 'does not equal input source: ', use_this_source)
            #print('Skipping index ', i, '...')
            continue
        
        #Get this observation's focus position.
        focus_position = s[i]['focus_position'][0]

        #if focus_position_input != None:
        #    if np.abs(focus_position - focus_position_input) > 0.001:
        #        print('\n\nWrong focus position.  Skipping index ', i, '...\n')
        #        continue
        if (quick == True) and (s[i]['source']=='rcw38') and (i%5 != 0): continue
        timeb = comptime()
        print('')
        print('')
        print('*******************     ', str(i),'/',str(n), '     ********************')

        this_s = s[i]
        print(this_s['starttime'])
        this_startdate = time.SptDatetime(this_s['mean_mjd'])
        print(this_startdate)
        this_mjd = this_s['mean_mjd']

        # Only consider HII observations that are between DATESTART and DATESTOP.
        if this_mjd < mjdstart: 
            print('1: This observation not between valid dates... skipping to next observation...')
            continue

        if this_mjd > mjdstop:
            print('1: This observation after the last tilt measurement... skipping to next observation...')
            continue

        # Get the pixel offsets valid for this observation
        rec = Receiver(date=this_startdate)

        #bolo_list = rec.getBoloIds()
        print('Making good bolo list...')
        bolo_list = this_s['id']
        common_xoff = []
        common_yoff = []
        for j in range(len(bolo_list)):
            common_xoff.append(rec[bolo_list[j][:-2]].pointing_offset[1]/60.) #Offsets read in as arcmin, not deg.
            common_yoff.append(rec[bolo_list[j][:-2]].pointing_offset[0]/60.) #Offests read in as arcmin, not deg.
        common_xoff = np.array(common_xoff)
        common_yoff = np.array(common_yoff)


        #Let's make some cuts on the bolos so the fits we make are only on the really
        #good ones.  We can be pretty harsh on the cuts.  We just need to make sure
        #we don't accidentally include any bad bolos.
        rres_cut = 0.008 #degrees
        cal2wnoise_cut = [2., 50.]
        elnod2wnoise_cut = [50., 70000.]
        amp2wnoise_cut = [100., 2000.]

        elnod_med_cut = 0.1


        # ==============================================
        # Find good bolos.
        # ==============================================

        this_cal = this_s['cal_response']
        this_elnod = this_s['elnod_response']
        this_elnod_med = np.median(this_s['elnod_response'][this_s['elnod_response'] != 0.])
        this_wnoise = this_s['wnoise']
        this_amp = this_s['amp']
        xdeg = this_s['xdeg']
        ydeg = this_s['ydeg']
        ID = this_s['id']

        #if test_mode == True:
        #    return xdeg, ydeg, ID, rec

        if skip_position_fit == False:
            print('Performing position fit cut...')
            wh_pos = np.where((xdeg == xdeg) & (ydeg == ydeg))[0]
            npos = len(wh_pos)

        #Calibrator cut
        if skip_cal == False:
            print('Performing calibrator cut...')
            if np.std(this_cal[this_cal==this_cal]) != 0.:
                wh_cal = np.array(np.nonzero((this_cal/this_wnoise > cal2wnoise_cut[0]) & \
                                             (this_cal/this_wnoise < cal2wnoise_cut[1]))[0])
                ncal = len(wh_cal)
            else:
                wh_cal = np.arange(1599, dtype=np.float)
                ncal = 1599
        else: 
            wh_cal = np.arange(1599, dtype=np.float)
            ncal = 1599

        #El nod cut
        if skip_elnod == False:
            print('Performing elnod cut...')
            #if np.std(this_elnod[this_elnod==this_elnod]) != 0.:
            #    wh_elnod = np.array(np.nonzero((this_elnod/this_wnoise > elnod2wnoise_cut[0]) & \
            #                                   (this_elnod/this_wnoise < elnod2wnoise_cut[1]))[0])
            #    nelnod = len(wh_elnod)
            wh_elnod = np.array(np.nonzero(np.abs(this_elnod - this_elnod_med) < elnod_med_cut*np.abs(this_elnod_med))[0])
            nelnod = len(wh_elnod)
            #else:
            #    wh_elnod = np.arange(1599, dtype=np.float)
            #    nelnod = 1599
        else: 
            wh_elnod = np.arange(1599, dtype=np.float)
            nelnod = 1599

        #source_scan amp response cut
        if skip_response == False:
            print('Performing source amplitude response cut...')
            if s[i]['nscans'] > 50:
                wh_good_response = np.arange(1599, dtype=np.float)
                n_goodResponse = 1599
                print('This is a FAST scan...')
                response_cut_low = response_cut_low_long
                response_cut_high = response_cut_high_long
            else:
                print('This is a VERY FAST scan...')
                response_cut_low = response_cut_low_short
                response_cut_high = response_cut_high_short

            if np.std(this_amp[this_amp==this_amp]) != 0.:
                wh_good_response = np.array(np.nonzero((np.abs(this_amp) > response_cut_low) & \
                                                       (np.abs(this_amp) < response_cut_high))[0])
                n_goodResponse = len(wh_good_response)
            else:
                wh_good_response = np.arange(1599, dtype=np.float)
                n_goodResponse = 1599
        else: 
            wh_good_response = np.arange(1599, dtype=np.float)
            n_goodResponse = 1599

        #source_scan amp vs noise level
        if skip_ampVsNoise == False:
            print('Performing amplitude vs white noise cut...')
            if np.std(this_amp[this_amp==this_amp]) != 0.:
                wh_amp = np.array(np.nonzero((np.abs(this_amp/this_wnoise) > amp2wnoise_cut[0]) & \
                                             (np.abs(this_amp/this_wnoise) < amp2wnoise_cut[1]))[0])
                namp = len(wh_amp)
            else:
                wh_amp = np.arange(1599, dtype=np.float)
                namp = 1599
        else: 
            wh_amp = np.arange(1599, dtype=np.float)
            namp = 1599

        #Find where the various cuts intersect.
        wh_gb1 = np.intersect1d(wh_amp, np.intersect1d(wh_cal, wh_elnod))
        wh_gb2 = np.intersect1d(wh_gb1, wh_good_response)
        wh_gb = np.intersect1d(wh_gb2, wh_pos)
        ngb = len(wh_gb)
        if dodebug==True:
            print('NCAL: ', ncal)
            print('NELNOD: ', nelnod)
            print('NAMP: ', namp)
            print('N_goodResponse: ', n_goodResponse)
            print('NPOS: ', npos)
            print('N_LIVE_INTERSECT: ', ngb)

        ngb_cut = 5
        if ngb < ngb_cut:
            print('Only ', ngb, ' good bolos.')
            print('We require ', ngb_cut, ' good bolos, so skipping this observation.')
            print('This occured after ELNOD/CAL/AMP cuts, but before SPATIAL cuts.')
            continue

        #Now do some spatial cuts (on XOFF and YOFF).
        #xc = -0.16
        #yc = 0.22
        #dr = np.sqrt((xdeg-xc)**2. + (ydeg-yc)**2.)
        #wh_dr = np.nonzero(dr < 0.65)[0]

        ##If NSCANS < 50, consider this a "very_fast" "stripe scan", and
        ##only pay attention to those detectors near the stripes.
        ##if this_s['nscans'] < 50.:
        ##    stripe1 = 0.55
        ##    stripe2 = 0.00
        ##    d_stripe = 0.1 + 0.15
        ##    wh_stripe = np.nonzero((np.abs(ydeg-stripe1) < d_stripe) | \
        ##                           (np.abs(ydeg-stripe2) < d_stripe))[0]
        ##else: wh_stripe = np.arange(1599, dtype=np.float)
        
        #wh_stripe = np.arange(1599, dtype=np.float)


        #Let's only include bolos with 0.5 deg of nominal boresight.
        #This is a static cut, and doesn't care about the response to sources on 
        # any given observation.
        dr = np.sqrt((common_xoff)**2. + (common_yoff)**2.)
        wh_dr = np.nonzero(dr < 0.5)[0]


        #If NSCANS < 50, consider this a "very_fast" "stripe scan", and
        #only pay attention to those detectors near the stripes.
        if this_s['nscans'] < 50.:
            stripe1 = 0.28 #Relative to boresight, NOT source
            stripe2 = -0.28 #Relative to boresight, NOT source
            d_stripe = 0.02#25
            wh_stripe = np.nonzero((np.abs(common_yoff-stripe1) < d_stripe) | \
                                   (np.abs(common_yoff-stripe2) < d_stripe))[0]
        else: wh_stripe = np.arange(1599, dtype=np.float)

        wh_spatial = np.intersect1d(wh_dr, wh_stripe)

        print(len(wh_spatial), ' bolos passing spatial cuts.')
        wh_gb = np.intersect1d(wh_spatial, wh_gb)
        ngb = len(wh_gb)
        
        if ngb < ngb_cut:
            print('Only ', ngb, ' good bolos.')
            print('We require ', ngb_cut, ' good bolos, so skipping this observation.')
            print('This occured after ELNOD/CAL/AMP/RADIUS cuts, but before fit to PHYSICAL MODEL.')
            continue
            
        this_gb = np.array(wh_gb, dtype=np.int)
        
        #Make sure to cut bolos from module C4.
        if cut_c4==True:
            wh_not_c4 = np.nonzero([bolo_id[0:2] != 'C4' for bolo_id in ID])[0]
            ngb = len(wh_not_c4)
            if ngb > 0. : this_gb = np.intersect1d(this_gb, wh_not_c4)
        ngb = len(this_gb)

        if ngb < ngb_cut:
            print('Only ', ngb, ' good bolos.')
            print('We require ', ngb_cut, ' good bolos, so skipping this observation.')
            print('This occured after ELNOD/CAL/AMP/RADIUS/MODULE cuts.')
            continue

        #Finally, let's not look at heaters and dark detectors for this analysis.
        wh_not_HDJ = np.nonzero([(bolo_id[-1] == 'X') | (bolo_id[-1] == 'Y') for bolo_id in ID])[0]
        ngb = len(wh_not_HDJ)
        if ngb > 0. : this_gb = np.intersect1d(this_gb, wh_not_HDJ)
        ngb = len(this_gb)

        if ngb < ngb_cut:
            print('Only ', ngb, ' good bolos.')
            print('We require ', ngb_cut, ' good bolos, so skipping this observation.')
            print('This occured after ELNOD/CAL/AMP/RADIUS/MODULE/HDJ cuts.')
            continue
        
        # ================================================
        # END finding good bolos.
        # ================================================

        # =================================================
        # Fit to a telescope imperfection model
        # =================================================
        mysort = np.sort(this_gb)
        all_xdeg = this_s['xdeg'][mysort]
        all_ydeg = this_s['ydeg'][mysort]
        all_id = mysort
        all_azact = this_s['az_actual'][mysort]
        all_elact = this_s['el_actual'][mysort]

        all_source = []
        for j in range(ngb):
            all_source.append(this_s['source'])
        #all_source = np.array(all_source)

        wh_rcw38 = np.nonzero(all_source == 'rcw38')[0]
        n_rcw38 = len(wh_rcw38)
        wh_mat5a = np.nonzero(all_source == 'mat5a')[0]
        n_mat5a = len(wh_mat5a)

        uniq_id = all_id
        n_uniq_id = len(all_id)

        nparam_base = 12
        nparam = nparam_base

        pa = []
        for j in range(nparam):
            pa.append({'value':1., 'fixed':0, 'mpside':2, 'relstep':0.01, \
                       'limited':[0,0], 'limits':[-1000., 1000.]})

        #Get the guesses for the fit parameters.
        guess = np.zeros(nparam)

        #These are boom flexure terms used for late SPTsz.  Also used for the first map run of Year-1 SPTPol.
        #a0_guess = -0.0098537632 
        #a1_guess = -0.021804936

        #2012 Point Source residual fits (boom_flex_fit.py). (Just new guard ring. no side shields).
        if (this_mjd <= 56262.):
            a0_guess = -0.00549382
            a1_guess = -0.02351464

        #2013 Point source residual fits. (After new guard ring and side shield).
        # 56262 == 01-Dec-2012 (after year-1 observations stopped).
        elif (this_mjd > 56262.) and (this_mjd < 57471.):
        
            #Determined from 2014 CMB field and radio point source residuals 
            #when determining new online parameters for EHT 20150109.
            a0_guess = 0.00549382
            a1_guess = -0.01758231

        #Calculated from 2016 point source offsets from Lindsey.
        elif (this_mjd > 57471.):
            a0_guess = 0.08068834
            a1_guess = 0.02130123



        #az encoder offset
        az0_guess = -0.304527
        pa[5]['fixed'] = 1 #Fix the az encoder offset value.

        #Do we want to fix the boom flex terms?
        if use_this_a6 == None:
            #If a6 isn't fixed (i.e., use_this_a6==None), then we want to fix the boom flex terms.
            pa[0]['fixed'] = 1 # Fix a0?
            pa[1]['fixed'] = 1 # Fix a1?
        else:
            pa[0]['fixed'] = 0 #Should be 0.
            pa[1]['fixed'] = 0 #Should be 0.
            pa[0]['limited'] = [1,1]
            pa[0]['limits'] = [-0.5, 0.5]
            pa[1]['limited'] = [1,1]
            pa[1]['limits'] = [-0.5, 0.5]

        a4_guess = -0.45
        pa[2]['fixed'] = 0 # Fix a4?
        pa[2]['limited'] = [1,1]
        pa[2]['limits'] = [-a4_limit, a4_limit]

        seed = comptime()
        ramp = 5./3600.

        if alter_HII_positions:
            #2012 deepfield cold focus
            if this_mjd <= 56220.3:
                offset_ra_sign = 1.
                offset_dec_sign = -1.

                #TARGET a4 = -0.00480843116706
                #TARGET a6 = -0.24813778

                #rcw38_x_guess = offset_ra_sign*(0.0) yields
                #med_a4 = -0.014181695
                #rcw38_y_guess = offset_dec_sign*(0.0) yields
                #med_a6 = -0.252596765


                #rcw38_x_guess = offset_ra_sign*(0.012) yields
                #med_a4 = -0.025227525
                #rcw38_y_guess = offset_dec_sign*(0.012) yields
                #med_a6 = -0.264596765


                #Original offsets were (-3.66", -6.22"), but I want (2.3", -0.1") for PS1.
                rcw38_x_guess = offset_ra_sign*(-0.0101829528)
                rcw38_y_guess = offset_dec_sign*(-0.001 + 0.0017)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)


            #2012 deepfield warm focus
            elif (this_mjd > 56220.3) and (this_mjd <= 56293):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = -0.0118614709
                #TARGET a6 = -0.23674354

                #rcw38_x_guess = offset_ra_sign*(0.0) yields
                #med_a4 = -0.02286047
                #rcw38_y_guess = offset_dec_sign*(0.0)
                #med_a6 = -0.24293826

                #rcw38_x_guess = offset_ra_sign*(0.005) yields
                #med_a4 = -0.0182708
                #rcw38_y_guess = offset_dec_sign*(0.005) yields
                #med_a6 = -0.23793826

                #Original PS1 offsets were (0.76", -1.3"), which is good enough!
                #Use sptpol_off_hii_params_20150302_201202.config
                rcw38_x_guess = offset_ra_sign*(0.0119823419767)
                rcw38_y_guess = offset_dec_sign*(-0.001)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)



            #2013 deepfield observations.
            #01-Jan-2013 - 30-Apr-2013
            elif (this_mjd > 56293.) and (this_mjd <= 56412.25):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = 0.0166711881015
                #TARGET a6 = -0.25329623

                rcw38_x_guess = offset_ra_sign*(0.010089560771416 - 0.00083)
                rcw38_y_guess = offset_dec_sign*(0.0)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #500d observations, starting 30-Apr-2013 through 29-Oct-2013:04:48:00.
            elif (this_mjd > 56412.25) and (this_mjd <= 56594.2):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = -0.028855554
                #TARGET a6 = -0.23026365

                #Corresponds to hii_config_20150204_201300
                rcw38_x_guess = offset_ra_sign*(0.010654954983)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.001480447)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #500d observations, 29-Oct-2013:04:48:00 through 26-Nov-2013:00:00:00.
            elif (this_mjd > 56594.2) and (this_mjd <= 56622.0):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = -0.02706584
                #TARGET a6 = -0.2250026

                #Corresponds to hii_config_20150204_201300
                rcw38_x_guess = offset_ra_sign*(0.013065831408)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.001480447)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #Summer field observations, 26-Nov-2013:00:00:00 through 07-Mar-2014:12:00:00.
            elif (this_mjd > 56622.0) and (this_mjd <= 56723.5):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = -0.0161141527531
                #TARGET a6 = -0.22705381

                #rcw38_x_guess = offset_ra_sign*(0.00823243) yields
                #med_a4 = -0.01359006

                #Corresponds to hii_config_20150204_201301
                rcw38_x_guess = offset_ra_sign*(0.00549943235)
                rcw38_y_guess = offset_dec_sign*(-0.002173 - 0.00085476)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #Summer field observations, 07-Mar-2014:12:00:00 through 25-Mar-2014:03:16:59.
            elif (this_mjd > 56723.5) and (this_mjd <= 56741.136794):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = 0.0076873926078
                #TARGET a6 = -0.2328389

                #rcw38_x_guess = offset_ra_sign*(0.00823243) yields
                #med_a4 = 0.00496712

                #Corresponds to hii_config_20150204_201301
                rcw38_x_guess = offset_ra_sign*(0.0111778442006)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.0040818)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #500d observations in 2014 after 25-Mar-2014:03:16:59.
            elif (this_mjd > 56741.136794) and (this_mjd <= 57023.):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #I'm trying to match the median a4 and a6 to what was fitted in Case 2
                #for EHT pointing.
                #a4_case2 = -0.0135806240929
                #a6_case2 = -0.21985015

                #rcw38_x_guess = offset_ra_sign*(0.00823243) yields
                #med_a4 = -0.013083635

                #Corresponds to hii_config_20150204_201400
                rcw38_x_guess = offset_ra_sign*(0.00769430793)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.0010465)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #Summer field observations 01-Jan-2015 - 27-Mar-2015.
            elif (this_mjd > 57023.) and (this_mjd <= 57108.308333333):
                
                #TARGET a4 = -0.00948234966
                #TARGET a6 = -0.22941241

                #rcw38_x_guess = offset_ra_sign*(0.0)
                #rcw38_y_guess = offset_dec_sign*(0.0)
                #together yield:
                

                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #rcw38_x_guess = offset_ra_sign*(0.00769430793 + 0.001112)
                #rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.0010465 + 0.0057523)

                rcw38_x_guess = offset_ra_sign*(0.00769430793)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.0010465 - 0.004808)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #Winter observations after 27-Mar-2015:07:24:00 and before 26-Oct-2015:14:21:46.
            #elif (this_mjd > 57108.308333333) and (this_mjd <= 57321.5984491):
            #    offset_ra_sign = -1.
            #    offset_dec_sign = 1.

            #    #TARGET a4 = -0.02407588285
            #    #TARGET a6 = -0.23175581

            #    #(I want residuals for ps1_1 of [4.6", 1.3"])

            #    rcw38_x_guess = offset_ra_sign*(0.010654954983)
            #    rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.001480447)
            #    mat5a_x_guess = offset_ra_sign*(0.0)
            #    mat5a_y_guess = offset_dec_sign*(0.0)

            #Winter observations after 27-Mar-2015:07:24:00 and before 26-Oct-2015:14:21:46.
            #Re-doing 2015 pointing after fixing position cut in very_fast_scan RCW38 fits.
            elif (this_mjd > 57108.308333333) and (this_mjd <= 57321.5984491):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                #TARGET a4 = -0.0223199298891
                #TARGET a6 = -0.23301919

                #(I want residuals for ps1_1 of [1.7", 0.25"])
                rcw38_x_guess = offset_ra_sign*(0.010654954983)
                rcw38_y_guess = offset_dec_sign*(-0.002173 + 0.001480447)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            #2016 observations:
            elif (this_mjd > 57321.5984491):
                offset_ra_sign = -1.
                offset_dec_sign = 1.

                rcw38_x_guess = offset_ra_sign*(0.0107)
                rcw38_y_guess = offset_dec_sign*(-0.001)
                mat5a_x_guess = offset_ra_sign*(0.0)
                mat5a_y_guess = offset_dec_sign*(0.0)

            

        else:
            rcw38_x_guess = 0.0
            rcw38_y_guess = 0.0
            mat5a_x_guess = 0.0
            mat5a_y_guess = 0.0

        pa[6]['fixed'] = 1 #Fit for HII region offsets? (0 means yes).
        pa[7]['fixed'] = 1 #Fit for HII region offsets? (0 means yes).
        pa[8]['fixed'] = 1 #Fit for HII region offsets? (0 means yes).
        pa[9]['fixed'] = 1 #Fit for HII region offsets? (0 means yes).

        #We fix a5 (the cross-dec collimation term) for each period of "constant focus"
        #or "constant optical bench position".  This needs to be udpdated each season!
        #A crude method for finding your new a5 is to run this program a few times with different
        #a5's using the 'quick' keyword and settling on the one with the following properties:
        #1) mean(t[where_rcw38]['a4']) ~= mean(t[where_mat5a]['a4'])
        #2) abs(t['a4']) ~< 0.01

        if use_this_a5 != None:
            a5_guess = use_this_a5

        #Choose a5 based on a4-vs-a5 source interception at given f0 focus position.
        if use_this_a5 == None:
            #100d 2012 observations.
            #01-Jan-2012 - 01-Dec-2012
            if (this_mjd <= 56262.) and (this_mjd >= 55927.):
                if np.abs(focus_position - -5.0) < 0.01:
                    a5_guess = 0.145835017388 # corresponds to a4 = -0.00480843116706
                elif np.abs(focus_position - -6.1) < 0.01:
                    a5_guess = 0.139079294488 # corresponds to a4 = -0.0118614709
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #100d 2013 observations.
            #26-Mar-2013 through 30-Apr-2013:06:00:00
            elif (this_mjd > 56293.) and (this_mjd < 56412.25):
                if np.abs(focus_position - -5.0) < 0.01:
                    a5_guess = 0.159899339897 # corresponds to a4 = 0.0166711881015
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #500d 2013 observations.
            #30-Apr-2013::06:00:00 through 26-Nov-2013
            elif (this_mjd > 56412.25) and (this_mjd < 56622):
                if np.abs(focus_position - -4.3) < 0.01:
                    a5_guess = 0.130988464637 # corresponds to a4 = -0.0288555536595
                elif np.abs(focus_position - -5.6) < 0.01:
                    a5_guess = 0.130327576827 # corresponds to a4 = -0.0270658402067
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #Summer field observations 2013-2014.
            #26-Nov-2013 - 27-Mar-2014:03:16:59
            elif (this_mjd > 56622.) and (this_mjd <= 56743):
                if np.abs(focus_position - -6.4) < 0.01:
                    a5_guess = 0.143134469079 # corresponds to a4 = -0.0161141527531
                elif np.abs(focus_position - -6.8) < 0.01:
                    a5_guess = 0.15832207864 # corresponds to a4 = 0.0076873926078
                elif np.abs(focus_position - -4.3) < 0.01:
                    a5_guess = 0.143133871539 # corresponds to a4 = -0.0135806240929
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #500d 2014 observations.
            #27-Mar-2014:00:00:00 - 01-Jan-2015
            elif (this_mjd > 56743.) and (this_mjd <= 57023.):
                if np.abs(focus_position - -4.3) < 0.01:
                    a5_guess = 0.143133871539 # corresponds to a4 = -0.0135806240929
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #Summer field observations 2015.
            #01-Jan-2015 - 27-Mar-2015.
            elif (this_mjd > 57023.) and (this_mjd <= 57108.308333333):
                if np.abs(focus_position - -6.4) < 0.01:
                    a5_guess = 0.148218267455 # corresponds to a4 = -0.00948234966081
                elif np.abs(focus_position - -6.8) < 0.01:
                    a5_guess = 0.148218267455 # corresponds to a4 = -0.00948234966081
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            ##500d 2015 observations (all cold focus)
            #elif (this_mjd > 57108.308333333) and (this_mjd <= 57321.5984491):
            #    if np.abs(focus_position - -4.3) < 0.01:
            #        a5_guess = 0.133716140436 # corresponds to a4 = -0.02407588285
            #    else:
            #        print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
            #        continue

            #500d 2015 observations (all cold focus)
            elif (this_mjd > 57108.308333333) and (this_mjd <= 57321.5984491):
                if np.abs(focus_position - -4.3) < 0.01:
                    a5_guess = 0.132872503136 # corresponds to a4 = -0.0223199298891
                elif np.abs(focus_position - -4.7) < 0.01:
                    a5_guess = 0.132872503136 # corresponds to a4 = -0.0223199298891
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

            #500d 2016 cold focus observations
            elif (this_mjd > 57473.25) and (this_mjd <= 57640.):
                if np.abs(focus_position - -4.3) < 0.01:
                    a5_guess = 0.090608 # matches to a4 = -0.0631255
                elif np.abs(focus_position - -4.7) < 0.01:
                    a5_guess = 0.090608 # matches to a4 = -0.0631255
                else:
                    print("Pointing paramter a5 has not be defined for focus position %f" % focus_position)
                    continue

        if use_this_a5== None:
            if fit_a5:
                pa[3]['fixed'] = 0 # a5 is free.
                pa[3]['limited'] = [1,1]
                pa[3]['limits'] = [-0.75, 0.75]
            else:
                pa[3]['fixed'] = 1 # Fix a5 with the value from the correct focus position.
        else:
            pa[3]['fixed'] = 1 # Fix a5 with the use_this_a5 value.



        #Set fit parameters for a6
        if use_this_a6 == None:
            pa[4]['fixed'] = 0 # Fix a6?
            pa[4]['limited'] = [1,1]
            pa[4]['limits'] = [-0.5, 0.0]
            pa[4]['limits'] = [-1.0, 1.0]
            a6_guess = -0.32
        else:
            pa[4]['fixed'] = 1
            a6_guess = use_this_a6

        
        #Now grab the az tilt measurements for this day (a2 and a3).
        fixed_aztilt = False
        if fixed_aztilt == False:
            print('... getting the interpolated aztilt...')
            aztilt = op.get_az_tilt_params(np.array([this_s['mean_mjd']]), config_file=aztilt_config_file)
            a2_guess = aztilt[0][0]
            a3_guess = aztilt[1][0]

            for item in pa[10:12]:
                item['fixed'] = 1 #Fix the az tilt parametes in the fit.

        #Get the linear sensor data.
        #Why was this set to 0?  It was neglecting the daz, del, and det terms below.
        do_linsens = 1
        daz = np.zeros(len(all_source)) + np.array([this_s['med_daz']])
        DEL = np.zeros(len(all_source)) + np.array([this_s['med_del']])
        det = np.zeros(len(all_source)) + np.array([this_s['med_det']])

        #daz -= np.median(daz)
        #DEL -= np.median(DEL)
        #det -= np.median(det)

        #Convert to degrees.
        daz /= 3600.
        DEL /= 3600.
        det /= 3600.

        #We'd like to record the linear sensor data, whether we use it or not.
        if len(daz) > 1:
            record_daz = np.median(daz)
            record_del = np.median(DEL)
            record_det = np.median(det)
        else:
            record_daz = daz[0]
            record_del = DEL[0]
            record_det = det[0]

        #If we're not using LINSENS, then set the values to zero.
        if do_linsens == 0:
            daz = np.zeros(len(all_source))
            DEL = np.zeros(len(all_source))
            det = np.zeros(len(all_source))

        #Record the LINSENS values that are actually being used.
        if len(daz) > 1:
            use_daz = np.median(daz)
            use_del = np.median(DEL)
            use_det = np.median(det)
        else:
            use_daz = daz[0]
            use_del = DEL[0]
            use_det = det[0]

        #Get the structure thermometry data and the pointing corrections
        # associated with it.
        do_thermo = 1

        #get the thermometry data
        if do_thermo:
            this_scu_temp = this_s['med_scu_temp']
        
            #get the linear sensor data (which goes into a joint thermo+linear sensor 
            #pointing model).
            daz2 = this_s['med_daz']
            del2 = this_s['med_del']
            det2 = this_s['med_det']
            lin = np.array([daz2, del2, det2])
            th = op.thermo2pointing(np.concatenate((this_scu_temp, lin)), [this_s['mean_mjd']])

            th_det = np.array(th['det'])[0]
            th_del = np.array(th['del'])[0]
        else:
            th_det = 0.
            th_del = 0.


        print('*** THERMO ***')
        print('TH_DET: ', th_det)
        print('TH_DEL: ', th_del)

        #Convert to degrees.
        th_del /= 3600.
        th_det /= 3600.

        print('TH_DET: ', th_det)
        print('TH_DEL: ', th_del)

        #Get the temperature and pressure for the day.
        airtemp = this_s['temp_avg']
        pressure = this_s['pressure_avg']

        sign = 1

        guess = [a0_guess, a1_guess, a4_guess, a5_guess, a6_guess, az0_guess,
                 rcw38_x_guess, rcw38_y_guess, mat5a_x_guess, mat5a_y_guess,
                 a2_guess, a3_guess]

        tofit = [all_xdeg, all_ydeg]
        myx = [all_azact, all_elact]

        for k in range(len(pa)):
            pa[k]['value'] = guess[k]

        #Define weight arrays based on the frequency of the bolos.
        tofite = np.zeros((2,len(all_xdeg)))
        wh_90 = np.nonzero([bolo_id[0:3] == 'ARG' for bolo_id in ID[this_gb]])[0]
        n90 = len(wh_90)
        wh_150 = np.nonzero([bolo_id[0:3] != 'ARG' for bolo_id in ID[this_gb]])[0]
        n150 = len(wh_150)
        
        if n90 > 0: 
            tofite[0][wh_90] = 1.5/3600.
            tofite[1][wh_90] = 1.5/3600.
        if n150 > 0: 
            tofite[0][wh_150] = 1.5/3600.
            tofite[1][wh_150] = 1.5/3600.

        #do_plot = 1
        xtol = 1e-7
        maxiter = 1000

        #Set up some variables to have in memory for fitting to the tele_imperfection model.
        com_airtemp = airtemp
        com_pressure = pressure
        com_sign = sign
        com_daz = daz
        com_del = DEL + th_del
        com_det = det + th_det
        com_source = all_source

        airtempK = airtemp + 273.15
        refract = quick_refr(all_elact, airtempK, pressure)
        
        final_xoff = common_xoff
        final_yoff = common_yoff

        #if test_mode == True:
        #    print('This is test mode.  Returning...')
        #    print('Parameter guesses: ', guess)
        #    return [xdeg, ydeg, common_xoff, common_yoff, this_gb, azc, elc, this_s, rec, ID]

        id_fit = all_id
        uniq_id_fit = uniq_id

        fa = {'x':myx, 'y':tofit, 'err':tofite, 'do_plot':do_plot}
        timea = comptime()
        print('...Fitting telescope imperfection model...')

        #Actually do the fit !
        result = mpfit.mpfit(mpfitfun(x=np.array(fa['x']), y=np.array(fa['y']), err=np.array(fa['err']), \
                                      com_airtemp=com_airtemp, com_pressure=com_pressure, \
                                      com_sign=com_sign, com_daz=com_daz, com_del=com_del, \
                                      com_det=com_det, com_source=com_source, id_fit=id_fit, \
                                      uniq_id_fit=uniq_id_fit, final_xoff=final_xoff, \
                                      final_yoff=final_yoff, \
                             do_plot=do_plot, return_reco=return_reco), \
                             parinfo=pa,xtol=xtol, maxiter=maxiter, quiet=True)

        #print out the results.
        print('------------------------------------------------------------------------------')
        print('Fitting to pair ', i, 'took ', comptime() - timea, ' seconds.')
        print('There were ', len(all_id), ' bolo-obs.')
        print('------------------------------------------------------------------------------')

        p1 = result.params
        ti = tele_imperfection(np.array(fa['x']),p1[0],p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8], \
                                        p1[9],p1[10],p1[11], \
                                        com_airtemp=com_airtemp, com_pressure=com_pressure, \
                                        com_sign=com_sign, \
                                        com_daz=com_daz, com_del=com_del, com_det=com_det, \
                                        com_source=com_source, id_fit=id_fit, uniq_id_fit=uniq_id_fit, \
                                        final_xoff=final_xoff, final_yoff=final_yoff, \
                                        xoff_meas=tofit[0], yoff_meas=tofit[1], \
                                        do_plot=do_plot, return_reco=True)

        resx = ti['x'] - all_xdeg
        resy = ti['y'] - all_ydeg
        resr = np.sqrt(resx**2. + resy**2.)

        #py.plot(resx, resy, 'k+')
        #py.show()

        sfr = 3600.

        mean_x = np.mean(resx)*sfr
        mean_y = np.mean(resy)*sfr

        print('STDDEV RESX: ', np.std(resx, ddof=1)*3600.)
        print('STDDEV RESY: ', np.std(resy, ddof=1)*3600.)
        print('Pointing Resid in arcsec: ', np.sqrt((np.std(resx)*3600.)**2 + (np.std(resy)*3600.)**2))

        #Skip this observation if for some reason either residuals are NaN.
        if np.std(resx, ddof=1)!=np.std(resx, ddof=1) or np.std(resy, ddof=1)!=np.std(resy, ddof=1): continue

        if test_mode == True:
            return {'az_c':ti['x'], 'el_c':ti['y'], 'xdeg':xdeg, 'ydeg':ydeg, \
                    'pix_xoff':common_xoff, 'pix_yoff':common_yoff, 'gb':this_gb}

        stddev_resx = np.std(resx)
        stddev_resy = np.std(resy)
        med_resr = np.median(resr)

        mean_resr = np.mean(resr)

        max_resx = np.max(np.abs(resx))
        max_resy = np.max(np.abs(resy))
        max_resr = np.max(resr)

        fit_error = np.zeros(nparam)
        for l in range(nparam):
            fit_error[l] = np.sqrt(result.covar[l][l])

        #Store the telescope imperfection params and pixel offsets into a 
        #big structure, along with the fit_errors.
        t[i]['starttime'] = this_s['starttime']
        t[i]['starttime_mjd'] = this_s['starttime_mjd']
        t[i]['mean_mjd'] = this_s['mean_mjd']
        
        t[i]['a0'] = result.params[0]
        t[i]['a1'] = result.params[1]
        t[i]['a4'] = result.params[2] - a4_offset

        print('This a0: ', result.params[0])
        print('This a1: ', result.params[1])
        print('This a4: ', result.params[2])
        print('This a5: ', result.params[3])
        print('This a6: ', result.params[4])
        print('This az0: ', result.params[5])

        t[i]['a5'] = result.params[3]
        t[i]['a6'] = result.params[4] - a6_offset
        t[i]['az0'] = result.params[5]

        t[i]['a0_error'] = fit_error[0]
        t[i]['a1_error'] = fit_error[1]
        t[i]['a4_error'] = fit_error[2]
        t[i]['a5_error'] = fit_error[3]
        t[i]['a6_error'] = fit_error[4]
        t[i]['az0_error'] = fit_error[5]

        t[i]['a2'] = result.params[10]
        t[i]['a3'] = result.params[11]
        t[i]['az_tilt_mag'] = np.sqrt(result.params[10]**2. + result.params[11]**2.)
        t[i]['az_tilt_deg'] = np.arctan(result.params[11]/result.params[10])*180./np.pi
        
        t[i]['temp_avg'] = this_s['temp_avg']
        t[i]['wind_speed_avg'] = this_s['wind_speed_avg']
        t[i]['wind_dir_avg'] = this_s['wind_dir_avg']
        t[i]['pressure_avg'] = this_s['pressure_avg']
        t[i]['focus_position'] = this_s['focus_position']

        t[i]['stddev_resx'] = stddev_resx
        t[i]['stddev_resy'] = stddev_resy
        t[i]['med_resr'] = med_resr
        t[i]['mean_resr'] = mean_resr
        t[i]['max_resx'] = max_resx
        t[i]['max_resy'] = max_resy
        t[i]['max_resr'] = max_resr

        gbtemp = np.zeros(1599, dtype=np.int16)
        gbtemp[uniq_id] = 1
        t[i]['gb'] = gbtemp

        t[i]['chisq'] = result.fnorm
        t[i]['dof'] = 2.*len(uniq_id) - (len(result.params) - len(np.where(result.perror != 0.)[0]))

        t[i]['mean_az'] = this_s['mean_az']

        t[i]['record_daz'] = record_daz
        t[i]['record_del'] = record_del
        t[i]['record_det'] = record_det

        t[i]['use_daz'] = use_daz
        t[i]['use_del'] = use_del
        t[i]['use_det'] = use_det

        t[i]['ngb'] = ngb
        t[i]['nscans'] = this_s['nscans']

        t[i]['meanx'] = mean_x
        t[i]['meany'] = mean_y

        t[i]['th_det'] = th_det
        t[i]['th_del'] = th_del
        
        t[i]['med_refract'] = np.median(refract)
        t[i]['stddev_refract'] = np.std(refract, ddof=1)

        t[i]['source'] = this_s['source']

        if summary_pkl == True:
            if quick == True:
                if (dosave == True) and (i%50==0.):
                    try:
                        with open(savename, 'ab') as savfile:
                            pk.dump(t, savfile, protocol=pk.HIGHEST_PROTOCOL)
                    except NameError:
                        with open(savename, 'wb') as savfile:
                            pk.dump(t, savfile, protocol=pk.HIGHEST_PROTOCOL)

                    else:
                        if (dosave == True) and (i%10==0.):
                            try:
                                with open(savename, 'ab') as savfile:
                                    pk.dump(t, savfile, protocol=pk.HIGHEST_PROTOCOL)
                            except NameError:
                                with open(savename, 'wb') as savfile:
                                    pk.dump(t, savfile, protocol=pk.HIGHEST_PROTOCOL)

    #If desired, write the fit results to a txt file.
    if dowrite==True:
        #First clean up t a little.
        tt = np.array(t).copy()
        print(type(tt))
        wh_chisq = np.array(np.where([item['chisq']/item['dof'] < 2000. for item in tt])[0])
        tt = tt[wh_chisq]

        if type(tt) != np.ndarray:
            if tt['starttime_mjd'] != 0.:
                if tt['mean_mjd'] != 0.:
                    if np.abs(tt['a4']) < a4_limit:
                        if (tt['stddev_resx'] < resid_limit) and (tt['stddev_resy'] < resid_limit):
                            pass
                        else: return t
                    else: return t
                else: return t
            else: return t
        else:
            wh_starttime_mjd = np.array(np.where([np.array(item['starttime_mjd']) != 0. for item in tt])[0])
            tt = tt[wh_starttime_mjd]
                
        if type(tt) != np.ndarray:
            if tt['mean_mjd'] != 0.:
                if np.abs(tt['a4']) < a4_limit:
                    if (tt['stddev_resx'] < resid_limit) and (tt['stddev_resy'] < resid_limit):
                        pass
                    else: return t
                else: return t
            else: return t
        else:
            wh_mean_mjd = np.array(np.where([item['mean_mjd'] != 0. for item in tt])[0])
            tt = tt[wh_mean_mjd]

        if type(tt) != np.ndarray:
            if np.abs(tt['a4']) < a4_limit:
                if (tt['stddev_resx'] < resid_limit) and (tt['stddev_resy'] < resid_limit):
                    pass
                else: return t
            else: return t
        else:
            wh_a4 = np.array(np.where([np.abs(item['a4']) < a4_limit for item in tt])[0])
            tt = tt[wh_a4]
            
        if type(tt) != np.ndarray:
            if (tt['stddev_resx'] < resid_limit) and (tt['stddev_resy'] < resid_limit):
                pass
            else: return t
        else:
            wh_goodstddev = np.array(np.nonzero([(item['stddev_resx'] < resid_limit) & \
                                                 (item['stddev_resy'] < resid_limit) for item in tt])[0])
            tt = tt[wh_goodstddev]

        if type(tt) == np.ndarray:
            print('Length of cleaned t: ', len(tt))

            utc = [item['mean_mjd'] for item in tt]
            a0 = [item['a0'] for item in tt]
            a1 = [item['a1'] for item in tt]
            a4 = [item['a4'] for item in tt]
            a5 = [item['a5'] for item in tt]
            a6 = [item['a6'] for item in tt]
            az0 = [item['az0'] for item in tt]
            f0 = [item['focus_position'][0] for item in tt]
            f1 = [item['focus_position'][1] for item in tt]
            f2 = [item['focus_position'][2] for item in tt]
            f3 = [item['focus_position'][3] for item in tt]
            f4 = [item['focus_position'][4] for item in tt]
            f5 = [item['focus_position'][5] for item in tt]
                
            tab = np.array(np.matrix(np.array([utc,a0,a1,a4,a5,a6,az0,f0,f1,f2,f3,f4,f5])).T)

            if writename == True: writename = ('pointing_param_'+filename_flag+'_'+use_this_source+'_a5_'+ 
                                                np.string0(use_this_a5)+'_a6_'+np.string0(use_this_a6)+'.txt')
            with open(writename, 'a') as out_file:
                for m in range(len(utc)):
                    out_file.write('%.4f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n' %
                                   (utc[m],a0[m],a1[m],a4[m],a5[m],a6[m],az0[m], f0[m],f1[m],f2[m],f3[m],f4[m],f5[m]))
                print("Wrote output to %s.\n" % writename)

        else:
            if len(tt) > 0:
                print('Length of cleaned t: ', len(tt))

                utc = tt['mean_mjd']
                a0 = tt['a0']
                a1 = tt['a1']
                a4 = tt['a4']
                a5 = tt['a5']
                a6 = tt['a6']
                az0 = tt['az0']
                f0 = tt['focus_position'][0]
                f1 = tt['focus_position'][1]
                f2 = tt['focus_position'][2]
                f3 = tt['focus_position'][3]
                f4 = tt['focus_position'][4]
                f5 = tt['focus_position'][5]
            
                tab = np.array(np.matrix(np.array([utc,a0,a1,a4,a5,a6,az0,f0,f1,f2,f3,f4,f5])))

                if writename == True: writename = ('pointing_param_'+filename_flag+'_'+use_this_source+'_a5_'+ 
                                                   np.string0(use_this_a5)+'_a6_'+np.string0(use_this_a6)+'.txt')
                with open(writename, 'a') as out_file:
                    for m in range(len(utc)):
                        out_file.write('%.4f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8\t%.8\tf%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n' % 
                                       (utc[m],a0[m],a1[m],a4[m],a5[m],a6[m],az0[m], f0[m],f1[m],f2[m],f3[m],f4[m],f5[m]))

                    print("Wrote output to %s.\n" % writename)
    #return t, p1
    return writename

                


def tele_imperfection(x,a0,a1,a4,a5,a6,az0,rcw38_x,rcw38_y,mat5a_x,mat5a_y,a2,a3, \
                      com_airtemp, com_pressure, com_sign, com_daz, com_del, com_det, \
                      com_source, id_fit, uniq_id_fit, final_xoff, final_yoff, \
                      xoff_meas, yoff_meas, \
                      do_plot=False, return_reco=False):

    uniq_id_fit = np.array(uniq_id_fit).tolist()
    n_uniq_id = len(uniq_id_fit)

    #Apply the thermometry + linear sensor correction to the EL TILT.
    #The sign is determined empirically. (Opposite that expected for SPTsz)
    #a4 += com_det

    #Unpack some of the variables.
    nbolo = len(x[0])
    az = x[0]
    el = x[1]

    airtemp = com_airtemp
    airtempK = airtemp + 273.15
    pressure = com_pressure
    refract = quick_refr(el, airtempK, pressure)

    xoff_model = np.zeros(nbolo)
    yoff_model = np.zeros(nbolo)

    #Create the source-dependent offset (b/c we were maybe able to
    #measure the position of a source, say mat5a, bettern than a previous experiment.
    source_offset_x = np.zeros(nbolo)
    source_offset_y = np.zeros(nbolo)
    whrcw38 = np.array(np.nonzero([item=='rcw38' for item in com_source])[0]).tolist()
    nrcw38 = len(whrcw38)
    whmat5a = np.array(np.nonzero([item=='mat5a' for item in com_source])[0]).tolist()
    nmat5a = len(whmat5a)

    if nrcw38 > 0:
        source_offset_x[whrcw38] += rcw38_x
        source_offset_y[whrcw38] += rcw38_y

    if nmat5a > 0:
        source_offset_x[whmat5a] += mat5a_x
        source_offset_y[whmat5a] += mat5a_y

    #Loop over each bolo, applying the pointing model & pixel offset relevant to each bolo.
    for k in range(n_uniq_id):
        thisaz_deg = np.array([az[k]])
        thisel_deg = np.array([el[k]])
        thisaz = thisaz_deg*np.pi/180.
        thisel = thisel_deg*np.pi/180.
        cos_az = np.cos(thisaz)
        sin_az = np.sin(thisaz)
        cos_el = np.cos(thisel)
        sin_el = np.sin(thisel)
        tan_el = np.tan(thisel)

        thisid = uniq_id_fit[k]

        sign = com_sign

        #This must match what is in the corresponding model in offline_pointing.py.
        az2 = thisaz_deg - (a2*cos_az + a3*sin_az)*tan_el \
                          + a4*tan_el - a5/cos_el - az0 \
                          + source_offset_x[k] \
                          - com_det[k]*tan_el

        el2 = thisel_deg + a0*sin_el + a1*cos_el \
                         + a2*sin_az  - a3*cos_az \
                         + source_offset_y[k] \
                         - com_del[k] - refract[k] \
                         - a6

        #Rotate boresight to (0,0)
        xel3, el3 = pt.rotate_boresight_to_zero_zero(thisaz_deg, thisel_deg, az2, el2)
        
        #Now add the pixel offsets
        xoff_model[k] = xel3 + final_xoff[thisid]
        yoff_model[k] = el3 + final_yoff[thisid]

    model = np.concatenate((xoff_model, yoff_model))
    meas = np.concatenate((xoff_meas, yoff_meas))

    #Plotting?
    if do_plot == True:
        py.plot(az, xoff_model - xoff_meas, 'b+')
        py.plot(az, yoff_model - yoff_meas, 'r+')
        py.show()
        
    if return_reco==True:
        reco = {'x':model[0:nbolo], 'y':model[nbolo:], 'daz':xel3, 'del':el3}
        return reco

    return model, meas



def mpfitfun(x,y,err,com_airtemp,com_pressure,com_sign, \
                     com_daz, com_del, com_det, com_source, \
                     id_fit, uniq_id_fit, final_xoff, final_yoff, \
             do_plot, return_reco):
    def f(p, fjac=None):
        a0 = p[0]
        a1 = p[1]
        a4 = p[2]
        a5 = p[3]
        a6 = p[4]
        az0 = p[5]
        rcw38_x = p[6]
        rcw38_y = p[7]
        mat5a_x = p[8]
        mat5a_y = p[9]

        a2 = p[10]
        a3 = p[11]

        xoff_meas = y[0]
        yoff_meas = y[1]

        concat_y = np.concatenate((xoff_meas, yoff_meas))
        concat_err = np.concatenate((err[0], err[1]))
        
        model, meas = tele_imperfection(x,a0,a1,a4,a5,a6,az0,rcw38_x,rcw38_y,mat5a_x,mat5a_y,a2,a3, \
                                        com_airtemp=com_airtemp, com_pressure=com_pressure, \
                                        com_sign=com_sign, \
                                        com_daz=com_daz, com_del=com_del, com_det=com_det, \
                                        com_source=com_source, id_fit=id_fit, uniq_id_fit=uniq_id_fit, \
                                        final_xoff=final_xoff, final_yoff=final_yoff, \
                                        xoff_meas=xoff_meas, yoff_meas=yoff_meas, \
                                        do_plot=do_plot, return_reco=return_reco)
    
        return [0, (np.array(concat_y) - np.array(model))/np.array(concat_err)]
    return f
        

def quick_refr(el, T, P):
    n0_minus_1 = (77.6e-6)*(P/T)
    delta_el = n0_minus_1/np.tan(el*np.pi/180.)*180./np.pi

    return delta_el
