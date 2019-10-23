import pdb
import numpy as np
import pylab as py
import cPickle as pk
import os
from glob import glob
from time import time as comptime
from sptpol_software.util import fits, time, tools
import sptpol_software.util.files as files
from sptpol_software.data.readout import SPTDataReader
import sptpol_software.analysis.offline_pointing as op

def consolidate_source_scan(start_date=None, stop_date=None, fitsdir='/data/sptdat/auxdata/', \
                            sources_to_use=['rcw38','mat5a'], filename=None, savename=None, doall=True, dosave=True):

    """
    Translated from RK's original consolidate_source_scan.pro.
    
    The purpose of this program is to consolidate the information contained in the
    source_scan fits files which is relevant to pointing.  The relevant info is stuffed
    into a dictionary and saved as a pkl file.

    ===Things that need to be saved in this pickle file===
    XDEG, YDEG, AMP, DC, ID
    STARTTIME,STARTTIME_MJD, STOPTIME, STOPTIME_MJD, MEAN_MJD
    SOURCE
    FILENAME
    TEMP_AVG, WIND_SPEED_AVG, WID_DIR_AVG, PRESSURE_AVG, TIPPER_TAU_MEAN
    AZ_ACTUAL, EL_ACTUAL
-    MEAN_AZ, MEAN_EL
    WNOISE, ELNOD_RESPONSE, ELNOD_SIGMA, CAL_RESPONSE, CAL_SIGMA
    LINEAR_SENSOR data (DET, DEL, DAZ, L1, L2, R1, R2)
    NSCANS, NFRAMES, NSAMPLES
    STRUCTURE THERMOMETRY
    ===

    Originally created by RK, 24 Sep 2008 in IDL.
    Translated to python by JWH, Oct 2012.
    """

    nbolo = 1599

    # Grab the filenames. Convert name strings to dates, and sort the filenames by date.
    filelist = np.array(tools.flatten([[glob(os.path.join(fitsdir, 'source/*'+srcname+'*.fits')) for srcname in sources_to_use], \
                                      [glob(os.path.join(fitsdir, 'source_vfast/*'+srcname+'*.fits')) for srcname in sources_to_use]]))
    filelist_times = time.filenameToTime(filelist)
    time_argsort = np.argsort(filelist_times)
    filelist = filelist[time_argsort]
    filelist_times = filelist_times[time_argsort]
    
    if start_date == None: 
        start_date = filelist_times[0]
    else:
        start_date = time.SptDatetime(start_date)
    if stop_date == None: 
        stop_date = filelist_times[-1]
    else:
        stop_date = time.SptDatetime(stop_date)
    
    #Figure out which filenames to read given the start and stop dates.
    times_to_use = np.array(map(lambda x: start_date <= x <= stop_date, filelist_times))
    print('Consolidating source scans between ', start_date, ' and ', stop_date, '.')

    if not times_to_use.any():
        print('There are no source scans in that date range...')
        return
    else:
        print('There are %d source observations in that date range.' % np.sum(times_to_use))
        filelist = filelist[times_to_use]
    nlist = len(filelist)
    
    #Record the system time
    clocka = comptime()

    #Restore the old dictionary if you gave filename to the function.
    if filename == None:
        filename = '/home/jhenning/sptpol_code/sptpol_software/analysis/pointing/source_scan_structure.pkl'
    if savename == None:
        savename = filename

    if doall == False:
        sold = pk.load(open(filename, 'r'))

    
    #Create an empty list of dictionaries to fill with the data that we need.
    s = []
    for i in range(nlist):
        
        s.append({'xdeg':np.zeros(nbolo), 'ydeg':np.zeros(nbolo),
         'amp':np.zeros(nbolo), 'dc':np.zeros(nbolo), 'id':np.zeros(nbolo, dtype=np.int16),
         'starttime':'', 'starttime_mjd':0.0, 
         'stoptime':'', 'stoptime_mjd':0.0, 'mean_mjd':0.0, 
         'source':'', 'filename':'',
         'temp_avg':0.0, 'wind_speed_avg':0.0, 'wind_dir_avg':0.0, 'pressure_avg':0.0,
         'tipper_tau_mean':0.0,
         'az_actual':np.zeros(nbolo), 'el_actual':np.zeros(nbolo),
         'mean_az':0.0, 'mean_el':0.0,
         'focus_position':np.zeros(6),
         'wnoise':np.zeros(nbolo), 'elnod_response':np.zeros(nbolo), 'elnod_sigma':np.zeros(nbolo),
         'cal_response':np.zeros(nbolo), 'cal_sigma':np.zeros(nbolo),
         'med_daz':0.0, 'med_del':0.0, 'med_det':0.0,
         'med_l1':0.0, 'med_l2':0.0, 'med_r1':0.0, 'med_r2':0.0,
         'nscans':0.0, 'nframes':0.0, 'nsamples':0.0,
         'med_scu_temp':np.zeros(60), 'mean_scu_temp':np.zeros(60)})

    #Now loop over each source scan and grab each of these data.
    nnew=0
    for i in range(nlist):
        timeb = comptime()
        print(str(i), '/', str(nlist-1))
        print('...checking ', filelist[i])

        #Check to see if the old dictionary list has this file's information.
        if doall == True: pass
        else:
            wh_already = np.where(filelist[i] == np.array([item['filename'] for item in sold]))[0]
            if len(wh_already) > 1:
                print('This file is in the source scan summary more than once!... ')
                print('Pausing so you can do something about it.')
                raw_input('Press ENTER to continue...')

            if len(wh_already)==1:
                print('This file has already been incorporated into the source scan summary... ')
                print('Packing old results and moving to next source scan.')
                #Stuff the new dictionary with the old information.
                s[i] = sold[wh_already[0]]
                continue


        #Okay, let's start grabbing information from new source scans.
        nnew += 1

        data = 0
        print('...reading FITS file...')
        print('...loading ', filelist[i])
        data = files.read(filelist[i])
        #nframes = data.observation.n_frames
        #nscans = data.observation.n_scans
        #nsamples = data.observation.n_samples

        #For now, load in the data from the archive for the observation time to grab
        #raw pointing information etc.  Eventually, these should be saved with the source fits.
        #extra = SPTDataReader(start_date=time.SptDatetime(data.header.start_date), \
        #                      stop_date = time.SptDatetime(data.header.stop_date))
        #extra.readData(verbose=False, timestream_units='watts', correct_global_pointing=False)

        #Double-check that the data has as observation structure.
        if hasattr(data, 'observation') == False:
            #print("   Missing an observation structure! Skipping observation %s." % data.header.start_date.arc)
            print("   Missing an observation structure! Skipping observation %s." % data.from_filename)
            nnew -= 1
            continue

        #Fill in auxiliary data for this scan.
        #Creat an SPTDataReader object to get ancillary information about the observation.
        da = SPTDataReader(start_date=data.header.start_date,
                           stop_date=data.header.stop_date, quiet=True)
        data.auxdata_keys = da.auxdata_keys
        data.rec = da.rec
        data.cuts = da.cuts
        data.telescope = da.telescope
        data.config = da.config
        try:
            #Let's not read in noise data.
            #files.fillAuxiliaryData(data, obs_types=['calibrator','elnod','noise'],
            files.fillAuxiliaryData(data, obs_types=['calibrator','elnod'],
                                    same_fridge_cycle=True, overwrite=False)
        except Exception as err:
            print("   Couldn't copy in auxiliary data: %s" % str(err))
            nnew -= 1
            continue
        
        #Finally, let's start pulling out some of the values we're looking for.
        xdeg = data.observation.bolo_x_offset
        ydeg = data.observation.bolo_y_offset
        bolo_id = data.observation.bolo_id_index_in_arc_file
        starttime = data.observation.start_time
        stoptime = data.observation.stop_time
        starttime_mjd = time.SptDatetime(starttime).mjd
        stoptime_mjd = time.SptDatetime(stoptime).mjd
        mean_mjd = np.mean([starttime_mjd, stoptime_mjd])

        #Is this a valid source scan?
        wh = np.where(xdeg == xdeg)
        nwh = len(wh)
        if nwh == 0:
            print("   Not a valid source scan. Skipping observation %s." % data.header.start_date.arc)
            nnew -= 1
            continue
        
        #Make sure calibrator and noise auxiliary info is present.
        #Be sure to skip nan's otherwise np.std() will return nan no matter how many 
        #good bolo data there are.
        try:
            #if ( np.std(data.observation.bolo_cal_response[np.isfinite(data.observation.bolo_cal_response)]) == 0. or
            #     np.std(data.observation.bolo_psd_white_noise_level[np.isfinite(data.observation.bolo_psd_white_noise_level)] == 0.0)):
            if ( np.std(data.observation.bolo_cal_response[np.isfinite(data.observation.bolo_cal_response)]) == 0.):
                #print("   Missing calibrator or noise data. Skipping this observation."
                print("   Missing calibrator data. Not worrying about noise... Skipping this observation.")
                nnew -= 1
                continue
        except AttributeError as err:
            # If we couldn't get calibrator and/or noise data, skip this observation.
            print(err)
            nnew -= 1
            continue
            

        #Now grab the linear sensor and thermometer data for this observation time.
        ls_time1 = time.SptDatetime(mean_mjd)
        ls_time2 = time.SptDatetime(mean_mjd + 0.5/60./24.)
        ls_date = np.array([ls_time1,ls_time2])

        print('... getting the linear sensor and structure thermometry data...')
        try:
            ls = op.get_lin_sens(ls_date, HII_fit=True)
        except (ValueError, TypeError) as err:
            print(err)
            nnew -= 1
            continue
        
        nls = len(ls['daz'])
        if nls > 1:
            med_daz = np.median(ls['daz'])
            med_del = np.median(ls['del'])
            med_det = np.median(ls['det'])
            med_l1 = np.median(ls['l1'])
            med_l2 = np.median(ls['l2'])
            med_r1 = np.median(ls['r1'])
            med_r2 = np.median(ls['r2'])
            med_scu_temp = np.zeros((60))
            mean_scu_temp = np.zeros((60))
            for k in range(len(ls['temp'])):
                med_scu_temp[k] = np.median(ls['temp'][k])
                mean_scu_temp[k] = np.mean(ls['temp'][k])
            
        else:
            med_daz = 0.
            med_del = 0.
            med_det = 0.
            med_l1 = 0.
            med_l2 = 0.
            med_r1 = 0.
            med_r2 = 0

        #Fill the output dictionary.
        s[i]['xdeg'] = xdeg
        s[i]['ydeg'] = ydeg
        s[i]['amp'] = data.observation.bolo_source_response
        s[i]['dc'] = data.observation.bolo_source_dc
        s[i]['id'] = data.observation.bolo_id_ordered
        s[i]['starttime'] = starttime
        s[i]['starttime_mjd'] = starttime_mjd
        s[i]['stoptime'] = stoptime
        s[i]['stoptime_mjd'] = stoptime_mjd
        s[i]['mean_mjd'] = mean_mjd
        s[i]['source'] = data.header.source
        s[i]['filename'] = filelist[i]
        s[i]['temp_avg'] = data.observation.temp_avg
        s[i]['wind_speed_avg'] = data.observation.wind_speed_avg
        s[i]['wind_dir_avg'] = data.observation.wind_dir_avg
        s[i]['pressure_avg'] = data.observation.pressure_avg 
        s[i]['tipper_tau_mean'] = data.observation.tipper_tau_mean
        s[i]['focus_position'] = data.observation.scu_benchoff

        s[i]['az_actual'] = data.observation.bolo_az_actual
        s[i]['el_actual'] = data.observation.bolo_el_actual

        s[i]['mean_az'] = data.observation.mean_az
        s[i]['mean_el'] = data.observation.mean_el
        #s[i]['wnoise'] = data.observation.bolo_psd_white_noise_level
        s[i]['elnod_response'] = data.observation.bolo_elnod_response
        s[i]['elnod_sigma'] = data.observation.bolo_elnod_sigma
        s[i]['cal_response'] = data.observation.bolo_cal_response
        s[i]['cal_sigma'] = data.observation.bolo_cal_sigma
        s[i]['med_daz'] = med_daz
        s[i]['med_del'] = med_del
        s[i]['med_det'] = med_det
        s[i]['med_l1'] = med_l1
        s[i]['med_l2'] = med_l2
        s[i]['med_r1'] = med_r1
        s[i]['med_r2'] = med_r2
        s[i]['nscans'] = data.observation.n_scans
        s[i]['nframes'] = data.observation.n_frames
        s[i]['nsamples'] = data.observation.n_samples
        s[i]['med_scu_temp'] = med_scu_temp
        s[i]['mean_scu_temp'] = mean_scu_temp

        print('Source observation data from %s: %d scans.' % (data.header.start_date.archive, s[i]['nscans']))

        #Let's save every few observations in case of a crash or a break or something.
        if dosave==True and (nnew%5)==0:
            print('...saving...')
            #wh_save = np.where([item['filename'] for item in s] != '')[0]
            #n_save = len(wh_save)
            #if n_save > 0:
            #    s=s[wh_save]
            with open(savename,'wb') as pickle_file:
                pk.dump(s, pickle_file, protocol=pk.HIGHEST_PROTOCOL)

        timec = comptime()
        print('That observation took ', str(timec-timeb), ' seconds.')

    if dosave==True:
        print('...saving...')
    #    wh_save = np.where([item['filename'] for item in s] != '')[0]
    #    n_save = len(wh_save)
    #    if n_save > 0:
    #        s=s[wh_save]
        with open(savename,'wb') as pickle_file:
            pk.dump(s, pickle_file, protocol=pk.HIGHEST_PROTOCOL)

    print('Added ', str(nnew), 'new observations.')
    print('Making the big structure took ', \
          str(comptime() - clocka), ' seconds.')

    return
            
            
