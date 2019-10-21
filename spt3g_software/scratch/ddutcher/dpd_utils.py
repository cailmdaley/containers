import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle as pkl
import datetime
import pandas as pd
import scipy
from glob import glob
from copy import copy
import astropy.stats as astrostats
from spt3g import core, dfmux, calibration, coordinateutils, gcp, std_processing, mapmaker, mapspectra
from spt3g.util.extractdata import extract_keys, MultiAccumulator
from spt3g.mapspectra import map_analysis
from spt3g.util import fitting, framecombiner, stats

def check_completed(obs_id, source='ra0hdec-57.5'):
    '''
    Given a field observation id, will find the corresponding
    GCP log, and check for 'abort_schedule'.
    Will also check the minimum and maximum ra, dec in observation,
    and the number of bolometers with cal s/n above some threshold.
    '''
    # As an example, use obs_id 20495272 from Aug. 26th 2017
    # on scott+amundsen
    obs_id=str(obs_id)
    data_dir = '/spt/data/bolodata/downsampled/' + source
    gcplog_dir = '/spt/data/gcplogs/' # on scott+amundsen

    check_numbolos(obs_id,source=source)

    obs_datafiles = sorted(glob(os.path.join(data_dir, obs_id, '0*.g3')))
    #if len(obs_datafiles)==1:
    #    print('Only one .g3 file in directory. Cannot be full observation')
    #    return
    first = core.G3File(obs_datafiles[0])
    for fr in first:
        if 'WiringMap' in fr:
            wiring_map = fr['WiringMap']
        if 'ObservationStart' in fr:
            start_time = fr['ObservationStart']
            print('Obs started at '+str(start_time))
        if 'RawBoresightEl' in fr:
            el_i = fr['RawBoresightEl'][0]/core.G3Units.degrees
            hkmap = fr['DfMuxHousekeeping']
            tuned = 0
            for bolo in fr['RawTimestreams_I'].keys():
                status = dfmux.HousekeepingForBolo(hkmap, wiring_map, bolo)
                if status.state=='tuned' and status.carrier_amplitude!=0:
                    tuned += 1
            break
    print('Number of bolos tuned at start of obs: '+str(tuned))
    print('Initial Elevation = '+str(el_i))
    last = core.G3File(obs_datafiles[-1])
    for fr in last:
        if 'RawBoresightEl' in fr:
            el_f = fr['RawBoresightEl'][-1]/core.G3Units.degrees
            stop_time = fr['TrackerStatus'].time[-1]
            date_str = stop_time.GetFileFormatString()
    print('Final Elevation = '+str(el_f))
    print('Obs ended at '+str(stop_time))

    # Get list of all GCP logs, find where our obs fits in
    # Since I'm doing this for how the observation ended, check
    # for GCP log containing the END of the observation (not the start).
    gcp_logs = sorted(glob(gcplog_dir+'*.log')+
                      [gcplog_dir+date_str+'.log'])
    log_ind = gcp_logs.index(gcplog_dir+date_str+'.log') - 1
    print('Found GCP log '+gcp_logs[log_ind])
    logfile = open(gcp_logs[log_ind])
    tmp = 0
    aborted = False
    for line in logfile:
        if 'abort_schedule' in line:
            t = line.split(' ')[1]
            if abs(int(t.replace(':','')) -
                   int(stop_time.GetFileFormatString().split('_')[-1])) < 2:
                print('Schedule was aborted: ')
                aborted = True
                tmp += 1
                print(line)
        if tmp != 0:
            if 'Exiting schedule' in line:
                print(line)
                tmp -= 1
    if not aborted:
        print('Schedule completed without error.')
    logfile.close()

def check_numbolos(obs_id, source='ra0hdec-57.5', sn_cutoff = 20) :
    cal_dir = '/spt/user/production/calibration/calframe/' + source
    offlinecal_file = glob(os.path.join(cal_dir, str(obs_id)+'.g3'))
    print(str(obs_id)+':')
    if len(offlinecal_file)==0:
        print('No offline calibration file for this observation')
        return
    cal  = list(core.G3File(offlinecal_file))[0]
    
    if not 'CalibratorResponseSN' in cal:
        print('No calibrator response data in calframe')
    else:
        num_bolos = 0
        for bolo in cal['CalibratorResponseSN'].keys():
            if cal['CalibratorResponseSN'][bolo] > 20:
                num_bolos += 1
        print('Number of bolos with CalSN above '+str(sn_cutoff)+': '
              + str(num_bolos))
    if not 'ElnodSNSlopes' in cal:
        print('No Elnod response data in calframe')
    else:
        num_bolos = 0
        for bolo in cal['ElnodSNSlopes'].keys():
            if cal['ElnodSNSlopes'][bolo] > 20:
                num_bolos += 1
        print('Number of bolos with ElnodSN above '+str(sn_cutoff)+': '
              + str(num_bolos))

class JumpToScan(object):
    ''' 
    Does no Scan frame processing until specified Scan is reached.
    Proceeds to process total of num_scans Scans.
    '''
    def __init__(self, start_scan, num_scans=1):
        self.desired_scan = start_scan
        self.num_scans = num_scans
        self.current_scan = -1
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            if self.current_scan < self.desired_scan:
                self.current_scan += 1
                return False
            elif self.current_scan == self.desired_scan + self.num_scans:
                core.G3Pipeline.halt_processing()
                return []
            else:
                self.current_scan += 1
                return
        else:
            return
    
def GrabFridgeLog(start_time = None, stop_time = None,
                  obsid = None, source='ra0hdec-57.5'):
    """
    Given a start and stop time, or an observation id,
    will find the relevant fridge log and return it as a dictionary.
    """
    assert any((obsid != None , (start_time != None and
                                 stop_time != None)))

    if obsid != None:
        # Grab start and end times of obs
        if not isinstance(obsid, str):
            obsid = str(obsid)
        data_dir = '/spt/data/bolodata/downsampled/' + source
        obs_datafiles = sorted(glob(os.path.join(data_dir, obsid, '0*.g3')))
        first = core.G3File(obs_datafiles[0])
        for fr in first:
            if 'ObservationStart' in fr:
                start_time = fr['ObservationStart']
                stop_time = fr['ObservationStop']
                break

    # Set up the He10 fridge cryoboard register indices
    registers = {0:'UCHEAD', 1:'ICHEAD', 2:'HE4HEAD', 3:'HE4FB',
                4:'HE4PUMP', 5:'ICPUMP', 6:'UCPUMP', 7:'HE4SW',
                8:'ICSW', 9:'UCSW', 10:'UC_STAGE', 11:'LC_TOWER',
                12:'IC_STAGE', 13:'4K_HEAD', 14:'4K_SQUID_STRAP',
                15:'50K_HEAD'}
    keys = {}
    callback = {}
    # Use `callback` to convert returned G3Times to mjd,
    # don't do anything to thermometer values.
    for reg, therm in registers.items():
        keys[therm] = ['array','cryo', 'temperature', 0, reg]
        callback[therm] = None
        keys['Time'] = ['array','cryo','utc']
        callback['Time'] = lambda t: g3time_to_datetime(t)

    # Initialize the MultiAccumulator object
    fridge_temps = MultiAccumulator(keys = keys)#, callback = callback)

    # Start the pipeline
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader,
             start_time = start_time, stop_time = stop_time,
             basedir = '/spt/data/arc')
    pipe.Add(fridge_temps)
    pipe.Run()

    # Return dictionary
    return fridge_temps.extract_values()

def GrabWeatherLog(start_time = None, stop_time = None,
                  obsid = None, source='ra0hdec-57.5'):
    """
    Given a start and stop time, or an observation id,
    will find the relevant weather log and return it as a dictionary.
    """
    assert any((obsid != None , (start_time != None and
                                 stop_time != None)))

    if obsid != None:
        # Grab start and end times of obs
        if not isinstance(obsid, str):
            obsid = str(obsid)
        data_dir = '/spt/data/bolodata/downsampled/' + source
        obs_datafiles = sorted(glob(os.path.join(data_dir, obsid, '0*.g3')))
        first = core.G3File(obs_datafiles[0])
        for fr in first:
            if 'ObservationStart' in fr:
                start_time = fr['ObservationStart']
                stop_time = fr['ObservationStop']
                break

    # Set up weather registers of interest
    registers = ['airTemperature','internalTemperature','pressure',
                'relativeHumidity','windDirection','windSpeed']
    keys = {}
    callback = {}
    # Use `callback` to convert returned G3Times to mjd,
    # don't do anything to thermometer values.
    for reg in registers:
        keys[reg] = ['array','weather', reg]
        callback[reg] = None
        keys['Time'] = ['array','weather','utc']
        callback['Time'] = lambda t: t.mjd

    # Initialize the MultiAccumulator object
    weather = MultiAccumulator(keys = keys, callback = callback)

    # Start the pipeline
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader,
             start_time = start_time, stop_time = stop_time,
             basedir = '/spt/data/arc')
    pipe.Add(weather)
    pipe.Run()

    # Return dictionary
    return weather.extract_values()


def unit_relations(quant, unit=None, source=None, dec=None, 
                   az_speed=None, sky_speed=None, width=None):
    '''
    Perform a variety of quantity relations useful for maps.
    Available units:
    ----------------
    ell
    freq, Hz
    poly
    deg
    
    sky_speed = az_speed*cos(dec)
    For 1500d pt 1, sky_speed is 0.73 deg/s for all subfields.
    For 1500d pt 2, az_speed is 1.0 deg/s for all subfields.
    For 2018 500d, sky_speed is 0.54 deg/s
    '''
    units = ['ell','freq','poly','deg']
    quant=float(quant)
    if source==None and any((dec==None, width == None)):
        source = 'ra0hdec-57.5'
        print("You didn't full specify a sky patch (source or dec and width)")
        print("Using default ra0hdec-57.5")
    # Specify declination, field width, and sky speed
    source_dict={'ra0hdec-57.5':{'dec':57.5,'width':60,'sky_speed':0.54},
                'ra0hdec-44.75':{'dec':44.25,'width':100,'sky_speed':0.73},
                'ra0hdec-52.25':{'dec':52.25,'width':100,'sky_speed':0.73},
                'ra0hdec-59.75':{'dec':59.75,'width':100,'sky_speed':0.73},
                'ra0hdec-67.25':{'dec':67.25,'width':100,'sky_speed':0.73}
               }
    
    if source in source_dict.keys():
        dec = source_dict[source]['dec']
        width = source_dict[source]['width']
        sky_speed = source_dict[source]['sky_speed']
    if sky_speed==None and az_speed !=None:
        sky_speed = az_speed*np.cos(dec*np.pi/180.)
    # Should I create two variables? Name differently?
    width = width*np.cos(dec*np.pi/180.)
    if unit==None:
        print('Specify the type of unit you entered. Choices are')
        print(units)
    if unit=='ell':
        ell = quant
        deg = 180./ell
        freq = ell*sky_speed/360.
        poly = width/deg
    if unit =='freq' or unit =='Hz':
        freq = quant
        ell = 360.*freq/sky_speed
        deg = 180./ell
        poly = width/deg
    if unit == 'deg':
        deg = quant
        ell = 180./deg
        freq = ell*sky_speed/360.
        poly = width/deg
    if unit == 'poly':
        poly = quant
        deg = width/poly
        ell = 180./deg
        freq = ell*sky_speed/360.
        
    print('Degree-scale: %.3f\nEll: %.1f\nTimestream Freq: %.2f Hz\nPoly-order: %d'%(
        deg,ell,freq,poly))
    
def GrabGCPLog(start_time = None, stop_time = None,
                  obsid = None, source='ra0hdec-57.5'):
    """
    Given a start and stop time, or an observation id,
    will find the relevant gcp log and return a buncha values as a dictionary.
    """
    assert any((obsid != None , (start_time != None and
                                 stop_time != None)))

    if obsid != None:
        # Grab start and end times of obs
        if not isinstance(obsid, str):
            obsid = str(obsid)
        data_dir = '/spt/data/bolodata/downsampled/' + source
        obs_datafiles = sorted(glob(os.path.join(data_dir, obsid, '0*.g3')))
        first = core.G3File(obs_datafiles[0])
        for fr in first:
            if 'ObservationStart' in fr:
                start_time = fr['ObservationStart']
                stop_time = fr['ObservationStop']
                break

    # Set up weather registers of interest
    weather_reg = ['airTemperature','pressure',
                'relativeHumidity','windDirection','windSpeed']
    keys = {}
    callback = {}
    # Use `callback` to convert returned G3Times to mjd,
    # don't do anything to other values.
    for reg in weather_reg:
        keys[reg] = ['array','weather', reg]
        callback[reg] = None
        keys['Time'] = ['array','weather','utc']
        callback['Time'] = lambda t: t.mjd
        
    # Initialize the weather MultiAccumulator object
    weather = MultiAccumulator(keys = keys, callback = callback)
    
    # Set up the He10 fridge cryoboard register indices
    fridge_reg = {0:'UCHEAD', 1:'ICHEAD', 2:'HE4HEAD', 3:'HE4FB',
                4:'HE4PUMP', 5:'ICPUMP', 6:'UCPUMP', 7:'HE4SW',
                8:'ICSW', 9:'UCSW', 10:'UC_STAGE', 11:'LC_TOWER',
                12:'IC_STAGE', 13:'4K_HEAD', 14:'4K_SQUID_STRAP',
                15:'50K_HEAD'}
    keys = {}
    callback = {}
    # Use `callback` to convert returned G3Times to mjd,
    # don't do anything to thermometer values.
    for reg, therm in fridge_reg.items():
        keys[therm] = ['array','cryo', 'temperature', 0, reg]
        callback[therm] = None
        keys['Time'] = ['array','cryo','utc']
        callback['Time'] = lambda t: t.mjd

    # Initialize the fridge MultiAccumulator object
    fridge = MultiAccumulator(keys = keys, callback = callback)
    
    # Set up acu registers of interest. (-_rate registers are garbage)
    acu_reg = ['az_err','az_pos', 'el_err','el_pos']
    keys = {}
    callback = {}
    # Use `callback` to convert returned G3Times to mjd,
    # don't do anything to other values.
    for reg in acu_reg:
        keys[reg] = ['antenna0','acu', reg]
        callback[reg] = None
        
    # Initialize the acu MultiAccumulator object
    acu = MultiAccumulator(keys = keys, callback = callback)

    # Start the pipeline
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader,
             start_time = start_time, stop_time = stop_time,
             basedir = '/spt/data/arc')
    pipe.Add(weather)
    pipe.Add(fridge)
    pipe.Add(acu)
    pipe.Run()

    return {'fridge': fridge.extract_values(),
            'weather': weather.extract_values(),
            'acu': acu.extract_values(),}

def g3time_to_datetime(dt):
    dt = str(dt)
    dt_start,dt_end=dt.split('.')
    dt_end = dt_end[0:6]#need to truncate seconds decimal places
    this_dt = datetime.datetime.strptime(dt_start,"%d-%b-%Y:%H:%M:%S")
    this_dt_out = this_dt.replace(microsecond = int(dt_end))
    return mpl.dates.date2num(this_dt_out)

def read_camb_cls(filename):
    ''' Read in camb cl file'''
    with open(filename) as f:
        header = f.readline().strip('\n').strip('#').split(' ')
    header = list(filter(('').__ne__,header))
    df = pd.read_csv(filename, names=header, sep='\s+',comment='#')
    return df.to_dict(orient='list')

def knox_formula(beam_fwhm = 1.2*core.G3Units.arcmin, fsky=1500./41253, 
                 noise_level=21.*core.G3Units.uK*core.G3Units.arcmin,
                cl_file = None):
    if cl_file is None:
        spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
        cl_file = os.path.join(spt3g_software,
                               'simulations/camb/base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat')
    cls = read_camb_cls(cl_file)
    dl = np.array(cls['EE'])
    ell = np.array(cls['L'])
    cl = 2*np.pi*dl/(ell**2+ell)

    noise_level /= (core.G3Units.rad*core.G3Units.uK)
    w_inv = noise_level**2

    beam_sig = beam_fwhm/(2*np.sqrt(2*np.log(2)))
    beam = np.exp(((beam_sig/core.G3Units.rad)**2)*ell**2)

    sensitivity = np.sqrt(2/(fsky*(2*ell+1)))*(1 + beam*(w_inv*ell**2)/(cl*ell**2))

    return cls['L'], sensitivity

def ps_plot(coadds, band='150GHz', waf = '',
            plot_maps = False, plot_depth = False, plot_ps = False,
            plot_weight = False, return_data = False,
            ylim = (10,1000), xlim = (300, 5000)):
    '''
    Plot autospectra and depth for T,Q,U in uK-arcmin
    '''
    obs = core.G3File(coadds)
    right = core.G3Frame(core.G3FrameType.Map)
    left = core.G3Frame(core.G3FrameType.Map)
    for fr in obs:
        if 'Id' not in fr:
            continue
        if 'Right'+band+waf in fr['Id'].replace('-',''):
            for k in fr.keys():
                if k not in right:
                    right[k] = fr[k]
                else:
                    val = right[k]
                    del right[k]
                    right[k]=val+fr[k]
        if 'Left'+band+waf in fr['Id'].replace('-',''):
            for k in fr.keys():
                if k not in left:
                    left[k] = fr[k]
                else:
                    val = left[k]
                    del left[k]
                    left[k]=val+fr[k]

    combined = core.G3Frame(core.G3FrameType.Map)
    combined['T'] = left['T']+right['T']
    combined['Q'] = left['Q']+right['Q']
    combined['U'] = left['U']+right['U']
    combined['Wpol'] = left['Wpol']+right['Wpol']

    t,q,u = mapmaker.mapmakerutils.remove_weight(combined['T'], combined['Q'], combined['U'], combined['Wpol'])
    
    lt,lq,lu = mapmaker.mapmakerutils.remove_weight(left['T'], left['Q'], left['U'], left['Wpol'])
    rt,rq,ru = mapmaker.mapmakerutils.remove_weight(right['T'], right['Q'], right['U'], right['Wpol'])

    difference = core.G3Frame(core.G3FrameType.Map)
    difference['T'] = (lt - rt)/2.
    difference['Q'] = (lq - rq)/2.
    difference['U'] = (lu - ru)/2.
    difference['Wpol'] = left['Wpol']+right['Wpol']

    apod = mapspectra.apodmask.make_border_apodization(
            combined['Wpol'],apod_type='cos', radius_arcmin=30.)
#             ,zero_border_arcmin=10, smooth_weights_arcmin=5)

    ptsrc = mapspectra.apodmask.make_apodized_ptsrc_mask(
        t,'/home/ddutcher/scratch/1500d_3band_10sigma_ptsrc.txt')

    apod *= ptsrc
    
    if plot_weight:
        mid = int(np.shape(combined['Wpol'].TT)[1]/2)
        tot_w = np.asarray(combined['Wpol'].TT)
        plt.figure()
        plt.plot(tot_w[:,mid]/np.max(tot_w[:,mid]))
        plt.figure()
        plt.imshow(tot_w/tot_w.max(), origin='lower', cmap = plt.cm.gray,
                   vmin=0, vmax=1)
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        
    if plot_maps:
        mask = copy(apod)
        mask[np.where(mask==0.0)] = np.nan
        
        plt.figure(figsize=(8,6))
        plt.imshow(mask*t/core.G3Units.uK,cmap = plt.cm.gray, vmin = -300, vmax = 300, origin='lower')
        plt.title('T')
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        
        plt.figure(figsize=(8,6))
        plt.imshow(mask*scipy.ndimage.gaussian_filter(q,sigma=2)/core.G3Units.uK,
                   cmap = plt.cm.gray, vmin = -20, vmax = 20, origin='lower')
        plt.title('Q')
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()

        plt.figure(figsize=(8,6))
        plt.imshow(mask*scipy.ndimage.gaussian_filter(u,sigma=2)/core.G3Units.uK,
                   cmap = plt.cm.gray, vmin = -20, vmax = 20, origin='lower')
        plt.title('U')
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
    
    if plot_depth:
        combined_cls = map_analysis.calculateCls(combined, apod_mask = apod,
                                                 ell_min=50,ell_max=8000,delta_ell=50,
                                                 qu = True, flatten_pol=True)
        difference_cls = map_analysis.calculateCls(difference, apod_mask = apod,
                                                   ell_min=50,ell_max=8000,delta_ell=50,
                                                   qu = True, flatten_pol=True)

        plt.figure()
        ls=['-','--']
        for i, cls in enumerate([combined_cls, difference_cls]):
            plt.plot(cls['ell'], (cls['TT']**0.5 / (core.G3Units.arcmin * core.G3Units.uK)),
                     color = 'b',linestyle=ls[i], label='TT')
            plt.plot(cls['ell'], (cls['QQ']**0.5 / (core.G3Units.arcmin * core.G3Units.uK)),
                     color = 'g',linestyle=ls[i], label='QQ')
            plt.plot(cls['ell'], (cls['UU']**0.5 / (core.G3Units.arcmin * core.G3Units.uK)),
                     color = 'r', linestyle=ls[i], label='UU')
            if i == 1:
                tidx = np.where((cls['ell']>3000) & (cls['ell']<4000))[0]
                pidx = np.where((cls['ell']>3000) & (cls['ell']<4000))[0]
                qnoise = np.sqrt(np.mean(cls['QQ'][pidx])) / (core.G3Units.arcmin * core.G3Units.uK)
                unoise = np.sqrt(np.mean(cls['UU'][pidx])) / (core.G3Units.arcmin * core.G3Units.uK)
                print('3000 < l < 4000:')
                print('T ',np.sqrt(np.mean(cls['TT'][tidx])) / (core.G3Units.arcmin * core.G3Units.uK))
                print('3000 < l < 4000:')
                print('Pol', (qnoise+unoise)/2.)

        plt.loglog()
        plt.legend(fontsize='small')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('Ell')
        plt.grid(which='both')
        plt.ylabel('Depth (uK-arcmin)')
        
    if plot_ps:
        eb_cls = map_analysis.calculateCls(left, cross_map = right, apod_mask = apod,
                               ell_min=50,ell_max=5000,delta_ell=50,
                               qu = False, flatten_pol=True)
        fig, ax = plt.subplots()
        ax.loglog(eb_cls['ell'], (1/core.G3Units.uK**2)*((eb_cls['TT']))*(eb_cls['ell']*(eb_cls['ell']+1)/(2*np.pi)),
         color = 'b', label='TT')
        ax.loglog(eb_cls['ell'],  (1/core.G3Units.uK**2)*((eb_cls['EE']))*(eb_cls['ell']*(eb_cls['ell']+1)/(2*np.pi)),
                 color = 'g', label='EE')
        ax.loglog(eb_cls['ell'],  (1/core.G3Units.uK**2)*((eb_cls['BB']))*(eb_cls['ell']*(eb_cls['ell']+1)/(2*np.pi)),
                 color = 'r', label='BB')

        ax.legend(fontsize='small')
        ax.set_xlim(100,3000)
        ax.set_ylim(.01,10000)
        ax.set_xlabel('Ell')
        ax.grid(which='both',linestyle='--')
        ax.set_ylabel('$D_\ell$ $[\mu K ^2]$',fontsize=12)
        ax.set_title("E, B")
        plt.tight_layout()
        
    if return_data:
        return left, right
    
def plot_weight_hist(weight_data, title='', return_data = False):

    if isinstance(weight_data, str):
        map_stats = sorted(glob(weight_data))
        weight_data = dict()
        for pth in map_stats:
            field = pth.split('/')[-1].split('_')[0]
            obsid = pth.split('/')[-1].split('_')[1]

            f = open(pth, 'rb')
            data = pkl.load(f)
            f.close()

            if field not in weight_data.keys():
                weight_data[field] = dict()
            weight_data[field][obsid] = data['med_w']

    elif isinstance(weight_data, dict):
        pass
    else:
        raise ValueError('input not recognized')

    new_colors = ['C0', 'C1', 'C2', 'C3']
    fig, ax = plt.subplots(nrows = 4, figsize=(9,6),sharex=True)
    for ind,field in enumerate(sorted(weight_data.keys())):
        to_plot = []
        for obsid, med_w in weight_data[field].items():
            if not np.isnan(med_w):
                to_plot.append(med_w)
            else:
                print(obsid+' has nan in weights')
        std = astrostats.mad_std(to_plot)
        h = ax[ind].hist(to_plot,label=field, bins=20,
                        color = new_colors[ind]),#range=(0,0.3))
        ax[ind].axvline(np.median(to_plot),color='k',linestyle='--')
        ax[ind].axvline(np.median(to_plot)-cutsig*std,
                        color='r',linestyle='-')
        ax[ind].axvline(np.median(to_plot)+cutsig*std,
                        color='r',linestyle='-')
        ax[ind].legend(loc=1)
        ax[ind].set_xlim(0,5000)
        num_cut = len(np.where(abs(np.array(to_plot) - np.median(to_plot))
                               > cutsig*std)[0])
        print(field+': %s of %s cut'%(num_cut,len(to_plot)))
    plt.suptitle(title)

    if return_data:
        return weight_data


def plot_spt3g_map(input_map, map_id = None, pixel_mask = None,
                   combine_map_frames = False, fig = None,
                   title = '', suptitle = ''):
    if isinstance(input_map, str):
        g3file = core.G3File(input_map)
        if combine_map_frames == True:
            fc = framecombiner.MapFrameCombiner(fr_id = map_id)
            for frame in g3file:
                # This is to stop from combining bsmaps or ptsrc maps
                if 'Wunpol' in frame or 'Wpol' in frame:
                    fc(frame)
            input_map = fc.output_frame
        elif map_id is not None:
            found = False
            for frame in g3file:
                if 'Id' in frame and frame['Id']==map_id:
                    input_map = frame
                    found = True
                    break
            if not found:
                print("No Map frame with Id matching %s found."%map_id)
                return
        else:
            print("Using first valid Map frame found in .g3 file")
            for frame in g3file:
                if frame.type == core.G3FrameType.Map:
                    if 'T' in frame and ('Wpol' in frame or 'Wunpol' in frame):
                        input_map = frame
                        break

    if isinstance(input_map, core.G3Frame):
        wgt = None
        if 'Wunpol' in input_map:
            wgt = input_map['Wunpol']
        elif 'Wpol' in input_map:
            wgt = input_map['Wpol']
        if 'T' not in input_map or wgt is None:
            raise KeyError("'T' or weight keys not in map frame")
        if input_map['T'].is_weighted:
            tmap = mapmaker.mapmakerutils.remove_weight_t(
                input_map['T'], wgt)
        else:
            tmap = input_map['T']

        if pixel_mask is None:
            pixel_mask = np.asarray(wgt.TT) != 0

    if isinstance(input_map, coordinateutils.FlatSkyMap):
        if input_map.is_weighted:
            print('Input map is marked as weighted, '
                  'but plotting will proceed')
        tmap = input_map
    
        if pixel_mask is None:
            pixel_mask = np.isfinite(np.asarray(tmap))
    
    #plot in uK
    tmap /= core.G3Units.uK
    
    # Figure out an appropriate colorscale
    sdtemp = np.nanstd(np.array(tmap)[np.where(pixel_mask)])
    if sdtemp == 0 or not np.isfinite(sdtemp):
        raise ValueError(
            "Map has zero or non-finite RMS.")
    atemp = fitting.gaussfit_hist(
        np.array(tmap)[np.where(pixel_mask)], 1000,
        -5.*sdtemp, 5.*sdtemp, do_plot = False)
    
    # Now that we're done using it as an indexing bool,
    # use pixel_mask as plotting mask
    pixel_mask = np.array(pixel_mask, dtype = float)
    pixel_mask[np.where(pixel_mask == 0.)] = np.nan

    # Make the plot
    plt.figure(fig)
    ax = plt.gca()
    im = ax.imshow(tmap*pixel_mask,
                   origin='lower', cmap=plt.cm.gray,
                   interpolation = None,
                   vmin = -5.*atemp['params'][2],
                   vmax = 5.*atemp['params'][2])
    plt.xticks([])
    plt.yticks([])
    plt.title(title)
    plt.suptitle(suptitle)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0)

    cbar = plt.colorbar(im, cax = cax)
    cbar.set_label("$\mu K_{CMB}$", labelpad = -10)
    cbar.ax.tick_params(axis='y',width=0,length=0)
    plt.tight_layout()
