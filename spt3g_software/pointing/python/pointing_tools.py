import numpy as np
import os
from scipy.optimize import fmin
from time import localtime as lt
import datetime as dt
from glob import glob
import pickle as pk
import subprocess as sp
import shlex
import re

from spt3g import core, gcp

def find_newest_tilts(tilts_dir='/poleanalysis/sptdaq/azTilts'):
    """
    Look through the tilts directory, grabbing MJD of each measurement
    from the directory structure.  Return the most recent MJD.

    Arguments
    ---------
    tilts_dir : string, optional
        The absolute path to a directory containing az tilt data.

    Returns
    -------
    newest_tilts : tuple, 2 elements
        A 2-element tuple consisting of newest_time (spt3g.core.G3Time),
        the (start?) time of the newest tilt measurement, and obs_id
        (string), the observation ID of that measurement.

    TODO
    ----
    * Confirm that newest_time is the start time of the newest az tilt
        measurement
    """

    from spt3g.std_processing import obsid_to_g3time
    d = np.array(glob(os.path.join(tilts_dir, '*')))
    obsids = []
    dates = []
    for index, item in enumerate(d):
        if item[-5:] != 'plots':
            thisobsid = item.split('/')[-1]
            obsids.append(thisobsid)
            dates.append(obsid_to_g3time(thisobsid))

    if not len(obsids):
        return (core.G3Time('20170101_000000'), '0')

    mjds = np.asarray([date.mjd for date in dates])
    smjd = np.argsort(mjds)
    obs_id = obsids[smjd[-1]]
    newest_time = dates[smjd[-1]]
    newest_tilts = (newest_time, obs_id)

    return newest_tilts
    

def find_tilt_times(start_time='20170101_000000',
                    stop_time='20580101_000000',
                    runlog_location='/spt_data/gcplogs'):
    """
    Find start and stop times of all valid az tilt measurements between the given
    starting and ending times.
    """

    start_data = sp.check_output(
        'grep -h -A 100 "Performing az tilt measurement" %s'
        % os.path.join(runlog_location, '*.log'), shell=True)
    start_data = start_data.decode().split('--\n')

    tilt_times = []

    # az order
    azs = np.arange(360, -1, -22.5) - 270
    azs[azs < 0] += 360
    az_list = ['45:00:00', '90:00:00']
    az_list += ['%d:%02d:00' % (int(az), int(60 * (az - int(az)))) for az in azs]

    for chunk in start_data:

        # extract the start time
        start = dt.datetime.strptime(chunk.partition(' Performing')[0], '%y%m%d %H:%M:%S')
        start = core.G3Time(start.isoformat())

        # skip data outside of the requested data range
        if start < core.G3Time(start_time):
            continue
        if start > core.G3Time(stop_time):
            continue

        # select only data between two incrementObsId commands
        if chunk.count('incrementObsId') < 2:
            continue
        chunk = chunk.partition(' incrementObsId\n')[-1].partition(' incrementObsId\n')[0]
        lines = chunk.split('\n')

        # extract the stop time
        stop = dt.datetime.strptime(lines[-1], '%y%m%d %H:%M:%S')
        stop = core.G3Time(stop.isoformat())

        # find all the slews
        slews = [line.split(' ', 2)[-1] for line in lines if 'slew az=' in line]
        if len(slews) != len(az_list):
            continue

        # check that all the slews are in the expected order
        slews = [re.search('slew az=([0-9:\.]+),', slew).group(1) for slew in slews]
        if slews != az_list:
            continue

        # store
        print('Found az tilt at: ', start, stop)
        tilt_times.append((start, stop))

    return tilt_times

def get_az_tilt_start_files(start_time='20000101_000000', 
                            stop_time='21000101_000000', 
                            runlog_location='/spt_data/gcplogs',
                            az_tilt_fit_dir=''):
    """
    Take a list of start and stop dates of tilt measurements and ensure
    the dates match up (end date is ~20 minutes after start date).

    Arguments
    ---------
    start_time : string, optional
        The starting time of an azimuth tilt observation.
    stop_time : string, optional
        The ending time of an azimuth tilt observation.
    runlog_location : string, optional
        The absolute path to a directory containing GCP runlogs.
        The files in this directory are grepped to find all azimuth tilt
        measurements after the year 2012.
    az_tilt_fit_dir : string, optional

    Returns
    -------
    cleaned_start_name, cleaned_stop_name : tuple, two elements
        The names of the cleaned start and stop date files. Both
        elements of the returned tuple are strings.

    TODO
    ----
    * Figure out what to do with this:
        Writes four lists of start and stop times for az tilt measurements (two each).
        One is every measurement found in the logs, the other is cleaned to only include
        matching dates from 2012 onward. (This information would be better
        suited for inline comments rather than inside the docstring. Low
        level descriptions of operations within the function are more useful
        to the programmer than the end-user.)
    * Possibly remove this since it's now included in the docstring above:
        Returns the names of the cleaned start and stop date files.
    """

    #First pull a list of start dates for az tilt measurements from the logs
    startname = 'az_tilts_start_'+str(lt().tm_year)+ \
                str(lt().tm_mon).zfill(2)+str(lt().tm_mday).zfill(2)+'.txt'
    stopname = 'az_tilts_end_'+str(lt().tm_year)+ \
                str(lt().tm_mon).zfill(2)+str(lt().tm_mday).zfill(2)+'.txt'

    os.system('awk "/Performing az tilt measurement/" %s > %s' %(os.path.join(runlog_location, '*.log'), os.path.join(az_tilt_fit_dir, startname)))
    os.system('awk "/incrementObsId/ {cnt=2} (cnt && cnt--) && /slew az=270/" %s > %s' %(os.path.join(runlog_location, '*.log'), os.path.join(az_tilt_fit_dir, stopname))) 

    #Open the start and stop files and only keep entries from the years requested.
    dstart = open(startname,'r').read().split('\n')[:-1]
    dstop = open(stopname,'r').read().split('\n')[:-1]

    start_time = core.G3Time(str(start_time))
    stop_time = core.G3Time(str(stop_time))

    f = open(os.path.join(az_tilt_fit_dir,startname[:-4]+'_cleaned.txt'),'w')
    for i in range(len(dstart)):
        if int(dstart[i][0:2]) >= 16:
            f.write(dstart[i]+'\n')
    f.close()

    f = open(os.path.join(az_tilt_fit_dir,stopname[:-4]+'_cleaned.txt'),'w')
    for i in range(len(dstop)):
        if int(dstop[i][0:2]) >= 16:
            f.write(dstop[i]+'\n')
    f.close()

    #Now read in the cleaned start and stop times and make sure the dates line up.
    dstart = open(startname[:-4]+'_cleaned.txt','r').read().split('\n')[:-1]
    dstop = open(stopname[:-4]+'_cleaned.txt','r').read().split('\n')[:-1]

    #Loop over the longer list.
    iterate_length = np.max([len(dstart),len(dstop)])# - 1
    for i in range(iterate_length):
        good_index = False
        while not good_index:
            try:
                start_mjd = core.G3Time('20'+dstart[i][0:6]+'_'+dstart[i][7:9]+dstart[i][10:12]+dstart[i][13:15]).mjd
                stop_mjd = core.G3Time('20'+dstop[i][0:6]+'_'+dstop[i][7:9]+dstop[i][10:12]+dstop[i][13:15]).mjd
            except IndexError:
                break
            #Only accept date pairs longer than 11 minutes and shorter than 20.5 minutes.
            if ((stop_mjd - start_mjd) < 0.0142361111) and ((stop_mjd - start_mjd) > 0.007639):
                good_index = True
                continue
            else:
                print('Duration of observation date pair is bad!')

            #If the start and stop times don't match up, then there is at least one extra line in either
            #the start or stop times list.  We need to find it(them).

            #1) First scan the next 20 start dates to see if they match the current stop date.
            start_counter = 0
            limit = max(21,iterate_length)
            for j in range(i+1,i+limit):
                start_counter += 1
                try:
                    test_start_mjd = core.G3Time('20'+dstart[j][0:6]+'_'+dstart[j][7:9]+dstart[j][10:12]+dstart[j][13:15]).mjd
                except IndexError:
                    test_start_mjd = core.G3Time('20'+dstop[j][0:6]+'_'+dstop[j][7:9]+dstop[j][10:12]+dstop[j][13:15]).mjd
                if ((stop_mjd - test_start_mjd) < 0.0142361111) and ((stop_mjd - test_start_mjd) > 0.007639): 
                    print('Had to skip', start_counter, ' start dates to match tilt times for index', i)
                    break
                
            #If the start_counter stops at anything other than 20, then remove that many lines from
            #the start date list.
            if start_counter < 20:
                print('Removing ', start_counter, ' start dates.')
                counter = 0
                while counter < start_counter:
                    dstart.pop(i)
                    counter += 1

            #If the start_counter reaches 20, no start dates match this stop date.
            #Instead remove the stop date.
            if start_counter == 20:
                print('Index ',i, ' has a bad stop date.  Removing it.')
                dstop.pop(i)

    #If the logs are incomplete for some reason and we have extra start times, 
    #remove them.  Start and end times are sorted so extra start times should be 
    #after the last end time.
    if len(dstart) > len(dstop):
        dstart = dstart[0:len(dstop)]

    #Now save the cleaned dates to 'cleaned' date files.
    f = open(startname[:-4]+'_cleaned.txt','w')
    for i in range(len(dstart)):
        f.write(dstart[i]+'\n')
    f.close()

    f = open(stopname[:-4]+'_cleaned.txt','w')
    for i in range(len(dstop)):
        f.write(dstop[i]+'\n')
    f.close()

    return startname[:-4]+'_cleaned.txt', stopname[:-4]+'_cleaned.txt'


def get_az_dates(start_time_file, end_time_file):
    """
    Obtain a list of az bearing tilt measurement observation start and
    stop times to feed to az_tilt.
 
    Author: Jason W. Henning - Dec 14, 2012.

    Arguments
    ---------
    start_time_file : string
        The absolute path to a file that contains the log entry start
        lines for az tilt observations. See get_az_tilt_start_files().
    end_time_file : string
        The absolute path to a file that contains the log entry end
        lines for az tilt observations. See get_az_tilt_start_files().

    Returns
    -------
    start_times : array-like
        An array of SptDatetime start dates.
    stop_times : array-like
        An array of SptDatetime stop dates.

    TODO
    ----
    * I don't want to remove this function's author attribution because
        that feels dirty, but do we do this anywhere else in 3g
        software? We should do here whatever is standard elsewhere.
        But it only seems right for Jason to take his name off himself.
    """
    
    start_info = open(start_time_file, 'r').read().split('\n')[:-1]
    stop_info = open(end_time_file, 'r').read().split('\n')[:-1]

    start_times = []
    stop_times = []

    for i in range(len(start_info)):
        start_times.append(core.G3Time(start_info[i].split(' ')[0][-6:]+' '+start_info[i].split(' ')[1][0:8]))
        stop_times.append(core.G3Time(stop_info[i].split(' ')[0][-6:]+' '+stop_info[i].split(' ')[1][0:8]))

    start_times = np.array(start_times)
    stop_times = np.array(stop_times)
    
    return start_times, stop_times


def interp_over_dropouts(signal, whnodrop=None):
    """
    Interpolate over points where signal is identically zero.

    Arguments
    ---------
    signal : array-like
        An array-like element comprised of (floats?).
    whnodrop : array-like, optional
        Description

    Returns
    -------
    signal : int or array-like
        Description

    TODO
    ----
    * Improve descriptions of arguments
    * Write description of returned value
    """

    npts = np.size(signal)
    if whnodrop is None:
        whnodrop = np.where(np.ravel(signal != 0.))[0]
        nnodrop = np.size(whnodrop)
    else:
        nnodrop = np.size(whnodrop)
    if nnodrop < npts/2.:
        signal = -1
    else:
        indtemp = np.arange(npts)
        signal = np.interp(indtemp, indtemp[whnodrop], signal[whnodrop])

    return signal

def readPointingRegisters(start_time, stop_time, arcdir='/spt_data/arc/'):
    """
    Return a dictionary whose keys correspond to GCP registers that record data
    needed to calculate the pointing model.

    Arguments
    ---------
    start_time : type?
        Description
    stop_time : type?
        Description
    arcdir : string, optional
        The absolute path to the directory containing GCP arcfile data.

    Returns
    -------
    d : dictionary
        A dictionary whose keys contain data necessary for pointing
        corrections. The returned keys are:
        'horiz_mount', 'horiz_off', 'tilts', 'linsens_avg', 'utc',
        'endcoder_off', 'features', 'pressure', 'temp', 'scu_temp'

    TODO
    ----
    * Track down argument types
    * Write descriptions for arguments
    * Write general description of what this function does
    """
    from spt3g.std_processing import ARCTimerangeReader
    pipe = core.G3Pipeline()
    pipe.Add(ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe.Add(gcp.ARCExtract)
    
    d = {}
    d['horiz_mount_x'] = []
    d['horiz_mount_y'] = []
    d['horiz_off_x'] = []
    d['horiz_off_y'] = []
    d['tilts_x'] = []
    d['tilts_y'] = []
    d['linsens_avg_l1'] = []
    d['linsens_avg_l2'] = []
    d['linsens_avg_r1'] = []
    d['linsens_avg_r2'] = []
    d['encoder_off_x'] = []
    d['encoder_off_y'] = []
    d['features'] = []
    d['scu_temp'] = []
    d['pressure'] = []
    d['temp'] = []
    d['utc'] = []
    d['obs_id'] = None

    def fillPointingStructure(fr):
        if fr.type != core.G3FrameType.GcpSlow:
            return

        if d['obs_id'] == None:
            d['obs_id'] = fr['antenna0']['tracker']['obs_id'].value

        d['horiz_mount_x'].append(fr['TrackerPointing'].horiz_mount_x)
        d['horiz_mount_y'].append(fr['TrackerPointing'].horiz_mount_y)
        d['horiz_off_x'].append(fr['TrackerPointing'].horiz_off_x)
        d['horiz_off_y'].append(fr['TrackerPointing'].horiz_off_y)
        d['tilts_x'].append(fr['TrackerPointing'].tilts_x)
        d['tilts_y'].append(fr['TrackerPointing'].tilts_y)
        d['linsens_avg_l1'].append(fr['TrackerPointing'].linsens_avg_l1)
        d['linsens_avg_l2'].append(fr['TrackerPointing'].linsens_avg_l2)
        d['linsens_avg_r1'].append(fr['TrackerPointing'].linsens_avg_r1)
        d['linsens_avg_r2'].append(fr['TrackerPointing'].linsens_avg_r2)
        d['encoder_off_x'].append(fr['TrackerPointing'].encoder_off_x)
        d['encoder_off_y'].append(fr['TrackerPointing'].encoder_off_y)
        d['features'].append(fr['TrackerPointing'].features)
        d['scu_temp'].append(fr['TrackerPointing'].scu_temp)
        d['pressure'].append(fr['TrackerPointing'].telescope_pressure)
        d['temp'].append(fr['TrackerPointing'].telescope_temp)
        d['utc'].append(fr['TrackerPointing'].time)

    
    pipe.Add(fillPointingStructure)
    pipe.Run()

    d['horiz_mount'] = np.array([np.array(d['horiz_mount_x']).flatten()[::100], 
                              np.array(d['horiz_mount_y']).flatten()[::100]])/core.G3Units.deg
    d['horiz_off'] = np.array([np.array(d['horiz_off_x']).flatten()[::100], 
                            np.array(d['horiz_off_y']).flatten()[::100]])/core.G3Units.deg
    d['tilts'] = np.array([np.array(d['tilts_x']).flatten()[::100], 
                        np.array(d['tilts_y']).flatten()[::100]])/core.G3Units.arcsec
    d['linsens_avg'] = np.array([np.array(d['linsens_avg_l1']).flatten()[::100], 
                              np.array(d['linsens_avg_l2']).flatten()[::100], 
                              np.array(d['linsens_avg_r1']).flatten()[::100], 
                              np.array(d['linsens_avg_r2']).flatten()[::100]])
    d['encoder_off'] = np.array([np.array(d['encoder_off_x']).flatten(), 
                              np.array(d['encoder_off_y']).flatten()])/core.G3Units.deg
    d['features'] = np.array(d['features']).flatten()
    d['pressure'] = np.array(d['pressure']).flatten()/core.G3Units.mb
    d['temp'] = np.array(d['temp']).flatten()
    d['utc'] = np.array([i.mjd for i in np.array(d['utc']).flatten()[::100]])
    d['scu_temp'] = np.array(d['scu_temp']).T
    
    
    d.__delitem__('horiz_mount_x')
    d.__delitem__('horiz_mount_y')
    d.__delitem__('horiz_off_x')
    d.__delitem__('horiz_off_y')
    d.__delitem__('tilts_x')
    d.__delitem__('tilts_y')
    d.__delitem__('linsens_avg_l1')
    d.__delitem__('linsens_avg_l2')
    d.__delitem__('linsens_avg_r1')
    d.__delitem__('linsens_avg_r2')
    d.__delitem__('encoder_off_x')
    d.__delitem__('encoder_off_y')

    return d


def az_tilt(tilt_times, batch_mode=None, remove_first_nsample_90=1,
            only_use_last_sample=None, test_save=True, plotting=1,
            skip_to='01-Feb-2012:00:00:00', save_tilt_mag=False,
            az_tilt_fit_dir='/poleanalysis/sptdaq/azTilts',
            wait_for_data_retrieval=True,
            arcdir=None, skip_day_fraction=0.5):
    """
    Calculate the az axis tilt from bearing tiltmeter readings.

    This is a port of spt_analysis/pointing/az_tilt.m.

    Arguments
    ---------
    tilt_times : list of 2-tuples of G3Times
        List of start and end times for each az tilt measurement to analyze
    batch_mode : type?, optional
        Description
    remove_first_sample_90 : integer, optional
        Description
    only_use_last_sample : type?, optional
        Description
    test_save : boolean, optional
        Description
    plotting : integer, optional
        Description
    skip_to : string, optional
        Description
    save_tilt_mag : boolean, optional
        Description
    az_tilt_fit_dir : string, optional
        Description
    wait_for_retrieval : boolean, optional
        Description
    arcdir : type?, optional
        Description
    skip_day_fraction : float, optional
        Description

    Returns
    -------
    several_things : tuple, five elements
        Description

    TODO
    ----
    * Figure out what the original docstring is trying to convey, then
        add it to the docstring if it is relevant. OG docstring follows:
            This is a port of spt_analysis/pointing/az_tilt.m.
            It's used to calculate the az axis tilt from bearing tiltmeter readings.

            tstart,tend are starting/ending dates (UTC) in control system form
            - obtain using the "logfile" program as e.g.:
            logfile sch list 2006'*' | grep point
            e.g.:
            az_tilt('09-Aug-2006:19:07:08', '09-Aug-2006:19:21:07')
    * Figure out missing argument types
    * Write missing argument descriptions
    """
    import time as mytime
    #Get the local time
    rundate = mytime.ctime()

    #Check for function keywords
    print('Checking for function keywords')
    if batch_mode != None: batch_mode=1
    if remove_first_nsample_90 != None: remove_first_nsample_90=1
    if only_use_last_sample != None: only_use_last_sample=1
    if test_save != None: test_save = 1

    print(len(tilt_times), 'tilt measurements found in logs.  Fitting only measurements made after', skip_to)

    #Angular unit conversion.
    d2r = np.pi/180.

    output = []
    #Now let's loop over each date.
    counter = 0
    bad_counter = 0
    skip_to_counter = 0

    #Today's date in MJD
    today = core.G3Time()
    all_mjd = []
    all_tilt_lat = []
    all_tilt_ha = []
    all_tilt_mag = []
    all_tilt_angle = []
    for tstart, tend in tilt_times:
        print(tstart)
        if tstart < core.G3Time(skip_to):
            skip_to_counter += 1
            continue

        # Ignore any measurements taken today or yesterday: Assume it takes 48 hours for
        # data to be sent up North.
        if wait_for_data_retrieval:
            if core.G3Time(tstart).mjd > today.Now().mjd-skip_day_fraction: 
                print('Skipping measurement within 1/2 day of today:', core.G3Time(today))
                continue
        else:
            print('Considering all measurements...')
        
        # Load the registers we care about.
        try:
            if arcdir != None:
                d = readPointingRegisters(str(tstart),str(tend), arcdir=arcdir)
            else:
                d = readPointingRegisters(str(tstart),str(tend))
        except ValueError:
            #The list of start az dates is unordered.  If we find a start date where the data don't
            #exist in the archive, just skip the date and move on.
            print('Did not find tilt data starting', tstart)
            continue

        #Get rid of any bogus archive points.  We'll look for feature=1: that's good.
        print('Removing bad frames')
        d = framecut(d,d['features'] > 0)

        mjd = d['utc']
        mjd_avg = np.sum(mjd)/len(mjd)

        #Actual mount pointing directions not recorded.
        #Reconstruct by adding in encoder offsets.
        print('Reconstructing pointing direction.')
        obsaz = d['horiz_mount'][0] + d['encoder_off'][0]
        obsel = d['horiz_mount'][1] + d['encoder_off'][1]
        az = obsaz
        el = obsel

        #Put az values into 0-360 range.
        print('Correct az range')
        az[az < 0.] += 360.
        az[az > 360.] -= 360.

        #Tiltmeter readings in arcsec.
        #tilts is in mas but cal file scales this to arcsec.
        tilt_x = d['tilts'][0]
        tilt_y = d['tilts'][1]

        #Get the number of samples per azimuth.
        fi45 = np.where(np.abs(az-45) <= 3.0)[0]
        nsample = len(fi45)

        #If any of the tiltmeter values are exactly zero, consider this a failed measurement.
        print('Checking for tiltmeter values that are exactly zero.')
        n0 = len(np.where(tilt_y == 0.0)[0]) + len(np.where(tilt_y == 0.0)[0])
        if n0 != 0:
            bad_output = (np.zeros((5,1)) + np.nan)
            print('fail: 1 - At least one tiltmeter value is exactly zero.')
            print('Failed Date: ', tstart)
            bad_counter += 1
            print('Bad counter: ', bad_counter)
            continue

        #Remove first nsample samples at AZ=90, when the tiltmeters are still recovering from a big slew.
        if remove_first_nsample_90 == 1:
            print('Removing first nsample samples at Az=90')
            fi90 = np.where(np.abs(az-90) <= 2.0)[0]
            if len(fi90) == nsample*2:
                fi_bad = fi90[0:nsample]
                el[fi_bad] = -999
                fi_ok = np.where(el > -100)[0]

                el = el[fi_ok]
                az = az[fi_ok]
                tilt_x = tilt_x[fi_ok]
                tilt_y = tilt_y[fi_ok]
                mjd = mjd[fi_ok]
                obsaz = az
                obsel = el
            else:
                bad_output = (np.zeros((5,1))+np.nan)
                print('fail:1 - Tiltmeter values are exactly zero.')
                bad_counter += 1
                print('Bad counter: ', bad_counter)
                #if batch_mode == 1:
                    #exit()
                continue
          
        if remove_first_nsample_90==1 and only_use_last_sample==1:
            print('Only using last sample')
            ntmp = len(az)
            if ntmp % nsample == 0:
                #MATLAB is column-major while pythin is row-major.  To get reshape to work, use order='F'
                az = np.transpose(np.reshape(az, (nsample,ntmp/nsample), order='F'))
                el = np.transpose(np.reshape(el, (nsample, ntmp/nsample), order='F'))
                tilt_x = np.transpose(np.reshape(tilt_x, (nsample, ntmp/nsample), order='F'))
                tilt_y = np.transpose(np.reshape(tilt_y, (nsample, ntmp/nsample), order='F'))
                mjd = np.transpose(np.reshape(mjd, (nsample, ntmp/nsample), order='F'))

                az = az[nsample,:]
                el = el[nsample,:]
                tilt_x = tilt_x[nsample,:]
                tilt_y = tilt_y[nsample,:]
                mjd = mjd[nsample,:]

                obsaz = az
                obsel = el

            else:
                bad_output = (np.zeros((5,1))+np.nan)
                print('fail:1')
                bad_counter += 1
                print('Bad counter: ', bad_counter)
                continue

        #Display some info:
        print('Read in data points')
        print(' ')
        print(' ')
        print('utc (days)  obsaz    el     tilt_x   tilt_y')
        print(' ')
        for i in range(len(az)):
            print('%.5f %7.3f %7.3f %7.4f %7.4f' % (d['utc'][i], az[i], el[i], tilt_x[i], tilt_y[i]))
        print('Read in %d data points' % len(az))
        print(' ')

        #Fit tilt model
        #Start parameters for fit.
        #[x_amplitude, y_amplitude/x_amplitude, x_offset, y_offset, angle(deg), 2az_x_amplitude, 2az_angle(deg)]
        print('Fitting the data')

        par0 = [200.,1.,0.,0.,90.,5.,0.]
        data = [tilt_x, tilt_y, obsaz]

        SumSqTiltError = makeSS(data)
        fit_count=0
        warnflags=0.1
        while fit_count < 10 and warnflags > 0:
            par, fval, it, funcalls, warnflags = fmin(SumSqTiltError, par0, maxiter=1000,xtol=1e-8, full_output=1)
            par0 = par
            fit_count += 1

        #Get the rms error.
        print('Fit done.  Find rms error')
        rms = np.sqrt(fval/len(az))

        #What are the residuals?
        print('Calculating residuals')
        rx = tilt_x - (par[0]*np.sin((az-par[4])*d2r) + par[2]) - (par[5]*np.sin((2.*az - par[6])*d2r))
        ry = tilt_y - (par[0]*par[1]*np.sin((az - 90. - par[4])*d2r) + par[3]) - (par[5]*par[1]*np.sin((2.*(az - 90.)-par[6])*d2r))

        #If the residuals are exactly zero, consider this a failed measurement.
        print('Checking for zero residuals.')
        if np.sum(rx) == 0. or np.sum(ry) == 0.:
            bad_output = (np.zeros((5,1))+np.nan)
            print('fail:1 - Fit residuals are exactly zero.')
            bad_counter += 1
            print('Bad counter: ', bad_counter)
            continue

        #Convert tilt (mag, ang) to (lat, ha)
        print('Converting tilt units')
        tm = ((par[0] + par[0]*par[1])/2.)*(d2r/3600.) # Tilt magnitude in radians.
        ta = par[4]*d2r #tilt angle in radians.
        tilt_lat = np.arcsin(np.tan(tm*np.sin(ta))/np.tan(ta))
        tilt_ha = tm * np.sin(ta)/np.cos(tilt_lat)
        tilt_lat /= d2r
        tilt_ha /= d2r

        #Display some info
        print(' ')
        print('tilt_mag tilt_ang tilt_lat tilt_ha (all in deg)')
        print('%8.4f %8.4f %8.4f %8.4f' % (tm/d2r, ta/d2r, tilt_lat, tilt_ha))
        print(' ')
        print('2az_mag 2az_ang (all in deg)')
        print('%8.4f %8.4f' % (par[5]/3600., par[6]))
        print(' ')

        #Plot model, data, and residuals if turned on.
        if plotting==1:
            import pylab as py
            print('Plotting results')
            maz = np.array(range(360))
            mx = par[0]*np.sin((maz - par[4])*d2r) + par[2] + par[5]*np.sin((2.*maz - par[6])*d2r)
            my = par[0]*par[1]*np.sin((maz - 90. - par[4])*d2r) + par[3] + par[5]*par[1]*np.sin((2.*(maz - 90.) - par[6])*d2r)
            py.figure(figsize=(10,7), dpi=80)
            py.subplot(3,3,1)
            py.plot(az/10., tilt_x, 'k+')
            py.xlabel('az/10deg')
            py.ylabel("x tilt ('')")
            py.plot(maz/10., mx, 'k-')
            py.xlim([0,36])

            py.subplot(3,3,2)
            py.plot(az/10., rx, 'k+')
            py.xlabel('az/10deg')
            py.ylabel("x resid ('')")
            py.xlim([0,36])
            py.title('az tilt fit - %s to %s' % (tstart, tend) )
        
            py.subplot(3,3,4)
            py.plot(az/10., tilt_y, 'k+')
            py.xlabel('az/10deg')
            py.ylabel("y tilt ('')")
            py.plot(maz/10., my, 'k-')
            py.xlim([0,36])

            py.subplot(3,3,5)
            py.plot(az/10.,ry, 'k+')
            py.xlabel('az/10deg')
            py.ylabel("y resid ('')")
            py.xlim([0,36])

            #Shift y residuals -90 deg in az and plot with s residuals.
            saz = az-90.

            #Put saz values in 0-360 deg range.
            saz[saz < 0.] += 360.
            saz[saz > 360.] -= 360.

            #Plot x & y residuals on the same graph
            py.subplot(3,3,6)
            py.plot(az/10., rx, 'k+')
            py.xlabel('az/10deg')
            py.ylabel("resid ('')")
            py.plot(saz/10., ry, 'ko')
            py.xlim([0,36])

            #Print(Fit Summary)
            modpar = ['tilt mag', 'tilt ang', 'tilt lat', 'tilt ha', '2az mag', '2az ang','rms']
            disp_vals = [tm/d2r, ta/d2r, tilt_lat, tilt_ha, par[5]/3600., par[6], rms]
            for i in range(len(modpar)):
                py.figtext(0.7, 0.85-i*0.025, '%s = %f' % (modpar[i], disp_vals[i]))
        
            py.subplots_adjust(wspace=0.5, hspace=0.5)
            #py.show()
            if not os.path.exists(os.path.join(str(az_tilt_fit_dir),'plots')):
                os.makedirs(os.path.join(str(az_tilt_fit_dir),'plots'))
            py.savefig(os.path.join(str(az_tilt_fit_dir),'plots','az_tilt_fit'+str(tstart)+'_'+str(tend)+'.png'))

    
        #Display some info
        print(' ')
        print('brg_tilt=[%.4f %7.4f %7.4f]' % (mjd_avg, tilt_lat, tilt_ha))
        print(' ')

        #If test_save is on, spit out results in a txt file.
        if test_save == 1:
            if save_tilt_mag:
                out_file = open('az_tilt_fits_'+rundate+'.txt', 'a')
                out_file.write('%.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n' % (mjd_avg,tm/d2r,tilt_lat,tilt_ha,ta/d2r))
                out_file.close()
            else:
                out_file = open('az_tilt_fits_'+rundate+'.txt', 'a')
                out_file.write('%.4f\t%7.4f\t%7.4f\t%7.4f\n' % (mjd_avg,tilt_lat,tilt_ha,ta/d2r))
                out_file.close()
        
            print(' ' )
            print('brg_tilt=[%.4f %7.4f %7.4f]' % (mjd_avg, tilt_lat, tilt_ha))
            print('fail: 0')

    
        output.append([mjd_avg, tilt_lat, tilt_ha, rms, np.max([np.max(np.abs(rx)), np.max(np.abs(ry))])])
        print('fail: 0')
        counter += 1
        print('Good counter: ', counter)

        all_mjd.append(mjd_avg)
        all_tilt_lat.append(tilt_lat)
        all_tilt_ha.append(tilt_ha)
        all_tilt_mag.append(tm/d2r)
        all_tilt_angle.append(ta/d2r)

        #Pack results in a G3 file.
        frame = core.G3Frame(core.G3FrameType.Calibration)
        frame['mjd'] = mjd_avg
        frame['tiltLat'] = tilt_lat*core.G3Units.degrees
        frame['tiltHA'] = tilt_ha*core.G3Units.degrees
        frame['tiltMag'] = tm*core.G3Units.rad
        frame['tiltAngle'] = ta*core.G3Units.rad

        end_frame = core.G3Frame(core.G3FrameType.EndProcessing)
            
        tilt_obsID = int(tstart.time/core.G3Units.s + 1 - core.G3Time('20170101_000000').time/core.G3Units.s)

        if not os.path.exists(os.path.join(str(az_tilt_fit_dir),str(tilt_obsID))):
            os.makedirs(os.path.join(str(az_tilt_fit_dir),str(tilt_obsID)))
        a = core.G3Writer(filename=os.path.join(str(az_tilt_fit_dir),str(tilt_obsID),str(1).zfill(4)+'.g3'))
        a(frame)
        a(end_frame)

    print('Observations before last processed date: ', skip_to_counter)
    print('Bad observations: ', bad_counter)
    print('Good observations: ', counter)
    print('Observations unaccounted for: ', len(tilt_times) - skip_to_counter - bad_counter - counter)

    return all_mjd, all_tilt_lat, all_tilt_ha, all_tilt_mag, all_tilt_angle


def framecut(d, ind):
    """
    Remove all frames in a data structure d where the given logical
    condition 'ind' is false. If a frame is only partially correct,
    remove it.

    Arguments
    ---------
    d : type? (data structure?)
        Description
    ind : type? (logical condition?)
        Description

    Returns
    -------
    thing? : type?
        Description

    TODO
    ----
    * Track down argument types
    * Write description of arguments
    * Track down return value type
    * Write description of returned value
    """

    #Modify each key and remove the bad frames by indexing the boolean array.
    for i, key in enumerate(d.keys()):
        #To avoid dealing with keys that are different shapes, transpose the array to make the "frame"
        #dimension first.  All registers have number of frames as their last dimension, or their
        #first after the transpose.  This makes indexing the frames much easier at the end.
        d[key] = np.transpose(d[key])

        #For the first key, determine which frames are good and bad by making a boolean index array.
        if i==0:
            ind = np.transpose(ind)
            nframes = ind.shape[0]
            good_frame = np.zeros((nframes), dtype=bool)

            for j in range(nframes):
                #For conditions where we're filtering on timestream properties, each frame will have
                #have more than one value.
                if len(ind.shape) > 1:
                    test = ind[j][ind[j] == True]
                else:
                    #When filtering on flags/conditions, each frame will be have one value.
                    if ind[j] == True:
                        test = True
                    else:
                        test = False
                #If a frame is only partially good (True), remove it.
                if len(ind.shape) > 1:
                    if len(test) != len(ind[j]):
                        good_frame[j] = False
                    else:
                        #Only keep frames where all the points evaulate True.
                        good_frame[j] = True
                else:
                    if test == True:
                        good_frame[j] = True
                    else:
                        good_frame[j] = False

        #Take the good_frame index array and apply it to the key.
        if key != 'obs_id':
            d[key] = d[key][good_frame]
            #Transpose back to the original shape.
            d[key] = np.transpose(d[key])
    
    return d


def SumSqTiltError1(data, p):
    """
    Placeholder for function description.

    A conversion of the original SumSqTiltError.m Matlab code. It is
    kept here for ease of use of makeSS, which is hard to read.

    Arguments
    ---------
    data : type?
        Description
    p : type?
        Description

    Returns
    -------
    thing? : type?
        Description

    TODO
    ----
    * Write thorough description of function
    * Track down types for arguments
    * Write descriptions for arguments
    * Track down what is returned and its type
    * Write description for what is returned and its type
    """
    tilt_x = data[0]
    tilt_y = data[1]
    az = data[2]
    el = data[3]
    d2r = np.pi/180.

    #p is a vector of model parameters.
    #[x_amplitude, y_amplitude/x_amplitude, x_offset, y_offset, angle(deg),2az_x_amplitude, 2az_angle(deg)]

    #Function returns the sum of squares of tilt errors.

    ex = p[0]*np.sin((az - p[4])*d2r) + p[2] #Expected x tilt
    ex += (p[5]*np.sin((2.*az - p[6])*d2r))

    ey = p[0]*p[1]*np.sin((az - 90. - p[4])*d2r) + p[3]
    ey += (p[5]*p[1]*np.sin((2.*(az-90.)-p[6])*d2r))

    d2 = (tilt_x - ex)**2. #Vector of x tilt errors.
    SS = np.sum(d2)

    d2 - (tilt_y - ey)**2.
    SS += np.sum(d2)

    return SS


def makeSS(data):
    """
    I need to be able to pass both data and parameters, p, but I only
    want to optimize parameters. I need to end up with a function with
    only parameters, p, as input.

    In the original Matlab script, the data were passed as global
    variables to a function that calculated the sum of the squared
    errors given a set of parameters. Instead, I use an intermediate
    lambda function to which I pass the measured data. A function is
    returned that uses these data built into the definition. I make the
    actual SumSqTiltError function by calling makeSS and giving it the
    measured data: SumSqTiltError = makeSS(data). I can then call the
    output function, passing a single variable, the parameter vector p.

    p is a vector of model parameters:
        [x_amplitude, y_amplitude/x_amplitude, x_offset,
        y_offset, angle(deg), 2az_x_amplitude, 2az_angle(deg)]

    data is a list of measured data, [tilt_x, tilt_y, az]. The model uses
    the p parameters to transform az into tilts. This function outputs a
    function that gives the sum of the squares of the errors between the
    model and the actual measured tilt values. We pass the output
    function to a minimizer to find the best values for p.

    Arguments
    ---------
    data : list
        A list of measured data in the form [tilt_x, tilt_y, az].

    Returns
    -------
    totError : type?
        Description

    TODO
    ----
    * Write the (already thorough) docstring in the more traditional
        dicstring form
    * Track down type of returned value
    * Write description of the returned value
    """
    tilt_x = data[0]
    tilt_y = data[1]
    az = data[2]

    d2r = np.pi/180.

    totError = lambda p: np.sum((tilt_x - (p[0]*np.sin((az - p[4])*d2r) + p[2] + (p[5]*np.sin((2.*az - p[6])*d2r))))**2.) \
        + np.sum((tilt_y - (p[0]*p[1]*np.sin((az - 90. - p[4])*d2r) + p[3] + (p[5]*p[1]*np.sin((2.*(az-90.)-p[6])*d2r))))**2.)

    return totError


def testError(x):
    """
    Placeholder description

    Arguments
    ---------
    x : number
        Description

    Returns
    -------
    thing? : number
        Description

    TODO
    ----
    * Write function description
    * Write argument description
    * Write returned value description
    """
    return x + (x-2)**2.

def quick_refr(el, T, P):
    """
    Calculate atmospheric refraction based on elevation, temperature,
    and pressure.

    Arguments
    ---------
    el : float
        The elevation of the telescope in units of degrees.
    T : float
        The outdoor temperature in units of C.
    P : float
        The barometric pressure in units of (what?).

    Returns
    -------
    delta_el : float
        The change in elevation due to atmospheric refraction.

    TODO
    ----
    * Improve function description
    * Improve argument descriptions
    * Improve returned value description
    """
    n0_minus_1 = (77.6e-6)*(P/T)
    delta_el = n0_minus_1/np.tan(el*np.pi/180.)*180./np.pi

    return delta_el


def rotate_boresight_to_zero_zero(az_bs, el_bs, az_in, el_in, inverse=False, static=True):
    """
    Placeholder description

    Arguments
    ---------
    az_bs : type?
        Description
    el_bs : type?
        Description
    az_in : type?
        Description
    el_in : type?
        Description
    inverse : boolean, optional
        Description
    static : boolean, optional
        Description

    Returns
    -------
    az_out, el_out : tuple, two elements
        Decription

    TODO
    ----
    * Write function description
    * Track down argument types
    * Write argument descriptions
    * Write returned value description
    """
    n_samples = len(az_bs)
    
    if static==True:
        n_bolos = len(az_in)
        az_in2 = np.zeros((n_bolos, n_samples))
        el_in2 = np.zeros((n_bolos, n_samples))
        
        for k in range(n_bolos):
            az_in2[k] = az_in[k]
            el_in2[k] = el_in[k]

        az_in = az_in2
        el_in = el_in2
        az_in2 = 0
        el_in2 = 0

    n_bolos = len(az_in)

    az_out = np.zeros((n_bolos, n_samples))
    el_out = np.zeros((n_bolos, n_samples))

    #Convert from az/el to phi/theta.  Our coordinate system is a standard
    # physics-stype spherical coordinate system.  Phi measures right-handed rotation
    # about the z-axis, starting at the x-axis.  Theta measures the angle from the z-axis
    # to the x-y plane.

    phi_bs = -az_bs
    theta_bs = 90. - el_bs
    phi_in = -az_in
    theta_in = 90. - el_in

    #Convert to radians
    phi_bs *= np.pi/180.
    theta_bs *= np.pi/180.
    phi_in *= np.pi/180.
    theta_in *= np.pi/180.

    #Convert (phi_in, theta_in, 1) to a 3D cartesian vector.
    x = np.cos(phi_in)*np.sin(theta_in)
    y = np.sin(phi_in)*np.sin(theta_in)
    z = np.cos(theta_in)

    x = x.reshape(-1, order='F')
    x = x.reshape(-1, n_samples, order='C')
    y = y.reshape(-1, order='F')
    y = y.reshape(-1, n_samples, order='C')
    z = z.reshape(-1, order='F')
    z = z.reshape(-1, n_samples, order='C')


    xvec = np.array([x, y, z])
    xvec = xvec.reshape(-1, order='F')
    xvec = xvec.reshape(n_samples, n_bolos,3, order='C')
                    

    #Create the rotation matrix
    if inverse == True:
        #Rotate coordinate system about the z-axis by phi_bs
        ang = phi_bs
        cs = np.cos(ang)
        sn = np.sin(ang)

        blank = np.zeros(len(cs))
        
        r1 = np.array([cs, -sn, blank,
                       sn, cs, blank,
                       blank, blank, blank+1.])

        #Rotate coordinate system about the Y-axis by -(pi/2 - theta_bs)
        ang = -(np.pi/2. - theta_bs)
        cs = np.cos(ang)
        sn = np.sin(ang)

        r2 = np.array([cs, blank, sn,
                       blank, blank+1., blank,
                       -sn, blank, cs])

    else:
        #Rotate coordinate system about the Y-axis by (pi/2 - theta_bs)
        ang = (np.pi/2. - theta_bs)
        cs = np.cos(ang)
        sn = np.sin(ang)
        blank = np.zeros(len(cs))

        r1 = np.array([cs, blank, sn,
                       blank, blank+1., blank,
                       -sn, blank, cs])

        #Rotate coordinate system about the z-axis by -phi_bs
        ang = -phi_bs
        cs = np.cos(ang)
        sn = np.sin(ang)
        
        r2 = np.array([cs, -sn, blank,
                       sn, cs, blank,
                       blank, blank, blank+1.])


    r1 = r1.reshape(-1,order='F')
    r1 = r1.reshape(-1,3,3, order='C')
        
    r2 = r2.reshape(-1,order='F')
    r2 = r2.reshape(-1,3,3, order='C')

    #Loop through time samples
    for i in range(n_samples):
        #Create the composite rotation matrix for this time sample
        thisr = np.matrix(r1[i])*np.matrix(r2[i])

        #Loop over the different bolos' pointing
        for j in range(n_bolos):
            outvec = np.array((thisr*np.matrix(xvec[i,j]).reshape(3,-1)).tolist())
            
            #Transform back to spherical coordinates
            thisaz = -np.arctan(outvec[1]/outvec[0])
            thisel = np.arcsin(outvec[2])

            #Stuff it
            az_out[j][i] = thisaz*180./np.pi
            el_out[j][i] = thisel*180./np.pi

    #print('Rotated az/el correction: ', np.array([az_out[0]]), np.array([el_out[0]]))

    return az_out, el_out

def get_pointing_fit_params(fit_file):
    '''
    Grab the HII pointing fit results for a particular .txt output file.

    Arguments
    ---------
    fit_file : string
        The absolute path to a file containing (something?).

    Returns
    -------
    lots_of_things : tuple, 13 elements
        Description

    TODO
    ----
    * Improve argument description
    * Write return value description
    '''

    data = open(fit_file, 'r').read().split('\n')[:-1]

    mjd = []
    a0 = []
    a1 = []
    a4 = []
    a5 = []
    a6 = []
    az0 = []
    f0 = []
    f1 = []
    f2 = []
    f3 = []
    f4 = []
    f5 = []

    for i in range(len(data)):
        mjd.append(np.float(data[i].split('\t')[0]))
        a0.append(np.float(data[i].split('\t')[1]))
        a1.append(np.float(data[i].split('\t')[2]))
        a4.append(np.float(data[i].split('\t')[3]))
        a5.append(np.float(data[i].split('\t')[4]))
        a6.append(np.float(data[i].split('\t')[5]))
        az0.append(np.float(data[i].split('\t')[6]))
        f0.append(np.float(data[i].split('\t')[7]))
        f1.append(np.float(data[i].split('\t')[8]))
        f2.append(np.float(data[i].split('\t')[9]))
        f3.append(np.float(data[i].split('\t')[10]))
        f4.append(np.float(data[i].split('\t')[11]))
        f5.append(np.float(data[i].split('\t')[12]))

    mjd = np.array(mjd)
    a0 = np.array(a0)
    a1 = np.array(a1)
    a4 = np.array(a4)
    a5 = np.array(a5)
    a6 = np.array(a6)
    az0 = np.array(az0)
    f0 = np.array(f0)
    f1 = np.array(f1)
    f2 = np.array(f2)
    f3 = np.array(f3)
    f4 = np.array(f4)
    f5 = np.array(f5)

    #Make sure everything is sorted.
    date_length = len(mjd)
    sorted_indices = sorted(range(date_length), key=lambda k: mjd[k])
    mjd = mjd[sorted_indices]
    a0 = a0[sorted_indices]
    a1 = a1[sorted_indices]
    a4 = a4[sorted_indices]
    a5 = a5[sorted_indices]
    a6 = a6[sorted_indices]
    az0 = az0[sorted_indices]
    f0 = f0[sorted_indices]
    f1 = f1[sorted_indices]
    f2 = f2[sorted_indices]
    f3 = f3[sorted_indices]
    f4 = f4[sorted_indices]
    f5 = f5[sorted_indices]

    return mjd,a0,a1,a4,a5,a6,az0,f0,f1,f2,f3,f4,f5

def write_tilt_config(tilt_fit_file, config_file=None, plotting=False, skip_last=0,
                      writename_prefix='sptpol_offline_tilts'):
    """
    Read in a list of az tilt fits, get rid of bad fits, and write the
    remaining az tilt fits to a config file.

    Arguments
    ---------
    tilt_fit_file : string
        The absolute path to a configuration file containing azimuth
        tilt fits.
    config_file : string, optional
        The absolute path to a (different?) configuration file containing
        (different?) azimuth tilt fits.
    plotting : boolean, optional
        If `True` the function will generate summary plots of the fits.
    skip_last : int, optional
        The number of dates to be skipped and disregarded in the fitting
        calculation.
    writename_prefix : string, optional
        A descriptive screen that may be prepended to the output
        configuration file name.

    Returns
    -------
    writename : string
        The name of the saved configuration file.

    TODO
    ----
    * Clarify what the config_file argument is supposed to be
    """
    
    mjd, a2, a3, tilt_ang = read_tilt_config(config_file=tilt_fit_file, file_type='txt', skip_last=skip_last)

    #Look for spurious points (negative fits, nans, angles very large or very small).
    mjd_cleaned = []
    a2_cleaned = []
    a3_cleaned = []
    tilt_ang_cleaned = []
    bad_counter = 0
    for i in range(len(mjd)):
        #if a2[i] < 0. or a3[i] < 0. or a2[i]!=a2[i] or \
        #a3[i]!=a3[i] or np.abs(tilt_ang[i]-60) > 30. or \
        #tilt_ang[i] < 0.:  # This was fine for normal operations, but not after we re-level the telescope!
        if not np.isfinite(a2[i]) or not np.isfinite(a3[i]):  # equiv to a2!=a2 type checks for nan above
            print('Skipping sample', i, 'at MJD', mjd[i])
            bad_counter += 1
        else:
            mjd_cleaned.append(mjd[i])
            a2_cleaned.append(a2[i])
            a3_cleaned.append(a3[i])
            tilt_ang_cleaned.append(tilt_ang[i])

    mjd_cleaned = np.array(mjd_cleaned)
    a2_cleaned = np.array(a2_cleaned)
    a3_cleaned = np.array(a3_cleaned)
    tilt_ang_cleaned = np.array(tilt_ang_cleaned)

    print('Number of cleaned tilt fits: ', len(mjd_cleaned))
    print('Number of bad tilt fits: ', bad_counter)

    #Find the earliest and latest dates
    date1 = core.G3Time(min(mjd_cleaned))
    date2 = core.G3Time(max(mjd_cleaned))

    print('Earliest date: ', date1)
    print('Latest date: ', date2)

    #Write a config file that contains parameters from the new fit.
    if config_file!=None:
        old_mjd,old_a2,old_a3,old_tilt_ang = read_tilt_config(config_file=config_file, \
                                                              file_type=config_file.split('.')[-1])

        startdate = core.G3Time(np.min([np.min(old_mjd),core.G3Time(date1).mjd]))
    else: startdate = date1

    writename = writename_prefix+'_'+str(lt().tm_year)+ \
                str(lt().tm_mon).zfill(2)+str(lt().tm_mday).zfill(2)+'_000000.config'
    f = open(writename,'w')
    f.write('VALIDITYBEG '+str(startdate)[:-7]+'\n')
    f.write('VALIDITYEND NEXT\n')
    f.write('\n')
    f.write('#Tilt parameters for offline pointing model.\n')
    f.write('#Last valid observation fit: '+str(date2)[:-7]+'\n')
    f.write('\n')
    f.write('#mjd	         a2	 a3     tilt_ang\n')

    if config_file!=None:
        for i in range(len(old_mjd)):
            f.write(str(old_mjd[i])+'\t'+str(old_a2[i])+'\t'+str(old_a3[i])+'\t'+\
                    str(old_tilt_ang[i])+'\n')
    
    for i in range(len(mjd_cleaned)):
        if config_file != None:
            if mjd_cleaned[i] <= old_mjd[-1]: continue
        f.write(str(mjd_cleaned[i])+'\t'+str(a2_cleaned[i])+'\t'+str(a3_cleaned[i])+'\t'+\
                str(tilt_ang_cleaned[i])+'\n')
    f.close()

    if plotting==True:
        import pylab as py
        print('Plotting summaries...')
        py.clf()
        py.figure(1)
        py.plot(mjd_cleaned, a2_cleaned*60., 'b+', label='Parameter a2')
        py.plot(mjd_cleaned, a3_cleaned*60., 'r+', label='Parameter a3')
        py.title('Tilt Parameters:  ' + str(date1)[0:11] + ' - ' + str(date2)[0:11])
        py.xlabel('Date [MJD]')
        py.ylabel('Tilt Parameter Fits [arcmin]')
        py.legend(loc='upper left')
        py.savefig('az_tilt_summary'+str(date1)[:-7]+'_'+str(date2)[:-7]+'.png')

        py.figure(2)
        py.clf()
        py.plot(mjd_cleaned, tilt_ang_cleaned, 'b+')
        py.title('Azimuth Bearing Tilt:  ' + str(date1)[0:11] + ' - ' + str(date2)[0:11])
        py.xlabel('Date [MJD]')
        py.ylabel('Az Bearing Tilt [deg]')
        py.savefig('az_tilt_angle_summary'+str(date1)[:-7]+'_'+str(date2)[:-7]+'.png')
        py.clf()

    print('Config file ' + writename + ' written.')

    return writename


def write_hii_config(fit_file, config_file=None, focus=None):
    '''
    Read in HII param fits and write a config file for only those
    matching the given f0 bench offset position.

    Arguments
    ---------
    fit_file : string
        The absolute path to a file containing HII parameter fits
    config_file : string, optional
        The absolute path to an existing configuration file containing
        HII region fits.
    focus : type?, optional
        Description

    Returns
    -------
    writename : string
        The name of the saved HII fits configuration file.

    TODO
    ----
    * Track down argument types
    * Write argument descriptions
    '''

    mjd,a0,a1,a4,a5,a6,az0,f0,f1,f2,f3,f4,f5 = get_pointing_fit_params(fit_file)

    if focus != None:
        good_focus = np.where(f0 == focus)[0]
    else:
        good_focus = np.arange(len(mjd))

    mjd = mjd[good_focus]
    a0 = a0[good_focus]
    a1 = a1[good_focus]
    a4 = a4[good_focus]
    a5 = a5[good_focus]
    a6 = a6[good_focus]
    az0 = az0[good_focus]

    startdate = core.G3Time(mjd[0])
    print('Start date for new observations: ', startdate)

    #Write a config file that contains parameters from the new fit.
    if config_file!=None:
        old_mjd,old_a0,old_a1,old_a4,old_a5,old_a6,old_az0 = read_hii_config(config_files=config_file)

        startdate = core.G3Time(np.min([np.min(old_mjd),core.G3Time(startdate).mjd]))
        print('Startdate considering old observations in config file:', startdate)

    writename = 'sptpol_offline_hii_params_'+str(lt().tm_year)+ \
                str(lt().tm_mon).zfill(2)+str(lt().tm_mday).zfill(2)+'_000000.config'

    with open(writename, 'w') as out_file:
        out_file.write('VALIDITYBEG '+str(startdate).split('.')[0]+'\n')
        out_file.write('VALIDITYEND NEXT\n')
        out_file.write('\n')
        out_file.write('##mjd	a0	a1     a4     a5     a6     az0\n')

        if config_file!=None:
            for i in range(len(old_mjd)):
                out_file.write(str(old_mjd[i])+'\t'+str(old_a0[i])+'\t'+str(old_a1[i])+'\t'+\
                    str(old_a4[i])+'\t'+str(old_a5[i])+'\t'+str(old_a6[i])+'\t'+str(old_az0[i])+'\n')

        for i in range(len(mjd)):
            if config_file != None:
                if mjd[i] <= old_mjd[-1]: continue
                if a5[i] < -0.3: continue
                
            out_file.write('%.4f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n' %
                           (mjd[i], a0[i], a1[i], a4[i], a5[i], a6[i], az0[i]))

    print('Config file ' + writename + ' written.')

    return writename


def read_tilt_config(config_file, file_type='config', skip_last=0):
    """
    Read in a list of az tilt fits, make sure they increase in MJD, and
    return parameters.

    Arguments
    ---------
    config_file : string
        The absolute path to a configuration file containing azimuth
        tilt fits.
    file_type : string, optional
        Must be either 'config' or 'text'. Defines how the function
        should parse the file located at the input path.
    skip_last : int, optional
        Description

    Returns
    -------
    several_things? : tuple, 4 elements
        Description

    TODO
    ----
    * Write argument descriptions
    * Write returned value description
    """

    data = open(config_file, 'r').read().split('\n')[:-1]

    mjd = []
    a2 = []
    a3 = []
    tilt_ang = []

    if file_type == 'config':
        for i in range(7,len(data)):
            mjd.append(np.float(data[i].split('\t')[0]))
            a2.append(np.float(data[i].split('\t')[1]))
            a3.append(np.float(data[i].split('\t')[2]))
            tilt_ang.append(np.float(data[i].split('\t')[3]))

    if file_type == 'txt':
        for i in range(len(data)):
            mjd.append(np.float(data[i].split('\t')[0]))
            a2.append(np.float(data[i].split('\t')[1]))
            a3.append(np.float(data[i].split('\t')[2]))
            tilt_ang.append(np.float(data[i].split('\t')[3]))

    mjd = np.array(mjd)
    a2 = np.array(a2)
    a3 = np.array(a3)
    tilt_ang = np.array(tilt_ang)

    #Sort the data into a list of increasing dates so we can easily append to earlier config files.
    date_length = len(mjd)
    sorted_indices = sorted(range(date_length), key=lambda k: mjd[k])
    mjd = mjd[sorted_indices]
    a2 = a2[sorted_indices]
    a3 = a3[sorted_indices]
    tilt_ang = tilt_ang[sorted_indices]

    #Remove the last few dates, if requested.
    if skip_last != 0:
        mjd = mjd[0:-skip_last]
        a2 = a2[0:-skip_last]
        a3 = a3[0:-skip_last]
    return mjd, a2, a3, tilt_ang

def read_hii_config(config_files):
    """
    Read in a list of old hii pointing parameters from a .config file.

    Arguments
    ---------
    config_files : type?
        Description

    Returns
    -------
    several_things? : tuple, 7 elements
        Description

    TODO
    ----
    * Track down argument type
    * Write argument description
    * Write returned value description
    """
    if type(config_files) == str:
        config_files = [config_files]

    mjd = []
    a0 = []
    a1 = []
    a4 = []
    a5 = []
    a6 = []
    az0 = []

    for j in range(len(config_files)):
        data = open(config_files[j], 'r').read().split('\n')[:-1]
        for i in range(4,len(data)):
            mjd.append(np.float(data[i].split('\t')[0]))
            a0.append(np.float(data[i].split('\t')[1]))
            a1.append(np.float(data[i].split('\t')[2]))
            a4.append(np.float(data[i].split('\t')[3]))
            a5.append(np.float(data[i].split('\t')[4]))
            a6.append(np.float(data[i].split('\t')[5]))
            az0.append(np.float(data[i].split('\t')[6]))

    mjd = np.array(mjd)
    a0 = np.array(a0)
    a1 = np.array(a1)
    a4 = np.array(a4)
    a5 = np.array(a5)
    a6 = np.array(a6)
    az0 = np.array(az0)

    #Sort the data into a list of increasing dates so we can easily append to earlier config files.
    date_length = len(mjd)
    sorted_indices = sorted(range(date_length), key=lambda k: mjd[k])
    mjd = mjd[sorted_indices]
    a0 = a0[sorted_indices]
    a1 = a1[sorted_indices]
    a4 = a4[sorted_indices]
    a5 = a5[sorted_indices]
    a6 = a6[sorted_indices]
    az0 = az0[sorted_indices]

    return mjd, a0, a1, a4, a5, a6, az0

def get_online_tilt_correction(obs_id, tilts_dir='/poleanalysis/sptdaq/azTilts'):
    """
    Convert the saved a2 and a3 az tilt pointing parameters to the form
    needed to correct the online pointing for az tilt. The line in GCP
    that corrects the az tilt reads:
    
    tilts tilt_ha*0.88, -tilt_lat*0.88, 0.00000

    where tilt_ha is saved as a3 and tilt_lat is saved as a2 in the 
    sptpol_offline_tilt config files.


    Arguments
    ---------
    obs_id : string or int
        A GCP observation ID.
    tilts_dir : string, optional
        The absolute path to a directory containing azimuth tilt
        observation IDs.

    Returns
    -------
    gcp_tilt_correction : string
        A string containing the new tilt correction for online pointing.

    TODO
    ----
    * Make sure the following original docstring block has been fully
        incorporated where appropriate:
        
        config_file: current sptpol_offline_tilts config file, containing
        mjd, a2, a3, and tilt_ang inputs.

    OUTPUTS:
        gcp_tilt_correction: A string containing the new tilt correction
                             for online pointing.
    """

    d = [f for f in core.G3File(os.path.join(tilts_dir,str(int(obs_id)),'0001.g3'))][0]
         
    a2 = d['tiltLat']
    a3 = d['tiltHA']

  
    
    new_a2 = 0.88*a3/core.G3Units.degrees
    new_a3 = -0.88*a2/core.G3Units.degrees
   
#    gcp_tilt_correction = 'tilts ' + str(new_a2) + ', ' + str(new_a3) + ', -0.01358062  #From get_online_tilt_correction from tilt '+str(int(obs_id))

# we are now fitting for el tilt, so get a4 from current file
    gcp_dir = os.getenv('GCP_DIR')
    initfile = os.path.join(gcp_dir,'config','init','pointing.init')
    f = open(initfile,'r')
    for line in f:
        if 'tilts' in line:
            if line[0] != '#':
                tiltline = line
    f.close
    el_tilt_str = tiltline.split(',')[-1].split('#')[0].strip()
    el_tilt = np.float(el_tilt_str)
    gcp_tilt_correction = ('tilts %.14f, %.14f, %.14f # From get_online_tilt_correction from tilt %d (%s)'
                           % (new_a2, new_a3, el_tilt, int(obs_id), obsid_to_g3time(int(obs_id))))
    # gcp_tilt_correction = 'tilts ' + str(new_a2) + ', ' + str(new_a3) + ', ' + el_tilt_str + '  #From get_online_tilt_correction from tilt '+str(int(obs_id))

    print('Append the following to gcp/config/init/pointing.init')
    print('and comment out the previous command:\n')
    print(gcp_tilt_correction)

    return gcp_tilt_correction


def plot_tilt_config(config_file, filename=None):
    """
    Plot the contents of a tilt configuration file.

    Arguments
    ---------
    config_file : string
        The absolute path to a configuration file containing fit
        parameters.
    filename : string, optional
        If `filename` is given, the plots generated by `plot_tilt_config()`
        are save to `filename`.
    """
    mjd, a2, a3, tilt_ang = read_tilt_config(config_file, file_type='config')

    mjd = mjd[1:]
    a2 = a2[1:]
    a3 = a3[1:]
    tilt_ang = tilt_ang[1:]

    import pylab as py
    py.figure(figsize=(8,9))
    
    py.subplot(3,1,1)
    # on the plot, they are in this order:
    py.plot(mjd, np.sqrt(a2**2+a3**2), label='mag')
    py.plot(mjd, a3, label='a3')
    py.plot(mjd, a2, label='a2')
    py.xlabel('mjd')
    py.ylabel('tilt component [deg]')
    py.legend(loc='center')
    
    py.subplot(3,1,2)
    py.plot(mjd, np.rad2deg(np.arctan2(a3,a2)))
    py.xlabel('mjd')
    py.ylabel('tilt arctan [deg]')
    
    py.subplot(3,1,3)
    py.scatter(a2, a3, c=mjd)
    py.colorbar()
    py.xlabel('a2')
    py.ylabel('a3')

    py.suptitle(config_file)

    if filename:
        print('Writing to:', filename)
        py.savefig(filename)



def get_point_source_offsets(source_file, snr_cut=6., offset_cut=3.):
    '''
    Read point source lists generated by Lindsey's IDL procedure and spit
    out delta_ra, delta_dec. This will be read in and fitted to the
    offline pointing model. Note list of offsets should be generated from
    maps made only with the 'tilts', 'linsen', and 'just_refraction' model
    flags.

    Returns delta_ra, delta_dec in units of degrees.

    Arguments
    ---------
    source_file : string
        The absolute path to a file containing (point source lists?)
    snr_cut : float, optional
        Description
    offset_cut : float, optional
        Description

    Returns
    -------
    several_things? : tuple, 6 elements
        A six element tuple where each element is a list

    TODO
    ----
    * Review function description
    * Write argument descriptions
    * Write better description of the returned values
    '''

    #Open the file and skip the first line, which is a header.
    d = open(source_file).read().split('\n')[1:-1]

    ra = []
    dec = []
    ara = []
    adec = []
    snr = []
    offset = []

    for i in range(len(d)):
        this_offset = float(filter(None, d[i].split(' '))[-1])
        this_snr = float(filter(None, d[i].split(' '))[3])
        if (this_offset <= offset_cut) and (this_snr >= snr_cut):
            ra.append(float(str(np.round(float(filter(None, d[i].split(' '))[1]),4))))
            dec.append(float(str(np.round(float(filter(None, d[i].split(' '))[2]),4))))
            ara.append(filter(None, d[i].split(' '))[4])
            adec.append(filter(None, d[i].split(' '))[5])
            snr.append(float(filter(None, d[i].split(' '))[3]))
            offset.append(float(filter(None, d[i].split(' '))[-1]))

    #Let's wrap the RAs to be +/- 180
    ra = [d-360. if d > 180 else d for d in ra]
    ara1 = [float(d)-360. if float(d) > 180 else float(d) for d in ara]
    adec1 = [float(d) for d in adec]

    ara = [str(np.round(d,5)) for d in ara1]

    ra = np.array(ra)
    dec = np.array(dec)
    ara1 = np.array(ara1)
    adec = np.array(adec)
    snr = np.array(snr)

    
    delta_ra = ra - ara1
    delta_dec = dec - adec1

    #Scale error by S/N.
    err_ra = np.array(offset)/snr
    err_dec = np.array(offset)/snr

    return list(ara), list(adec), list(delta_ra), list(delta_dec), list(err_ra), list(err_dec)
    

def plot_bench_positions(source_scan_file='/data/jhenning/pointing/source_scan_structure_20161209_201600.pkl'):
    '''
    Plot the f0 bench position for RCW38 observations to determine
    periods of constant focus.

    Arguments
    ---------
    source_file_scan : string, optional
        The absolute path to a pickle file containing (something?)

    Returns
    -------
    Nothing?

    TODO
    ----
    * Improve argument description
    '''

    d = pk.load(open(source_scan_file,'r'))

    mjd = []
    f = []

    for i in range(len(d)):
        if d[i]['source'] == 'rcw38':
            f.append(d[i]['focus_position'][0])
            mjd.append(d[i]['starttime_mjd'])

    f = np.array(f)
    mjd = np.array(mjd)

    import pylab as py
    py.plot(mjd, f, 'k+')

