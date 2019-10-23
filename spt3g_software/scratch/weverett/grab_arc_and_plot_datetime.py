from spt3g import core, gcp
import numpy as np
import pylab as pl
import glob as glob
import cPickle as pkl
from scipy.ndimage import gaussian_filter1d as gf1d
import matplotlib as mpl
import datetime as datetime

def grab_and_append_arc_frame_data(arc_file_list_txt='',arc_file_time_range=[],dat_dir='/spt_data/arc/',year='y2'):

    #usage, for example:
    #data_out = grab_and_append_arc_frame_data(arf_file_time_range=['20180306_000000','20180307_000000'],data_dir = 'path/to/arc/files')

    arc_file_list = []
    if arc_file_list_txt:
        print('getting arc files from txt file')
        arc_files = open(arc_file_list_txt,'r')
        arcs = arc_files.readlines()
        for aa in arcs:
            arc_file_list.append(aa.split('\n')[0])

    if arc_file_time_range:
        print('getting arc files from time range')
        if year == 'y1':
            all_tstamps = glob.glob(dat_dir+'2017*.dat')
        if year == 'y2':
            all_tstamps = glob.glob(dat_dir+'2018*.dat')
        all_tstamps_just_time = []
        for date in all_tstamps:
            all_tstamps_just_time.append(np.float(str(date.split('/')[-1].split('_')[0])+'.'+str(date.split('/')[-1].split('_')[1].split('.')[0])))
        all_tstamps_just_time = np.array(all_tstamps_just_time)
        this_time_range_start = np.float(str(arc_file_time_range[0].split('_')[0])+'.'+str(arc_file_time_range[0].split('_')[1]))
        this_time_range_stop = np.float(str(arc_file_time_range[1].split('_')[0])+'.'+str(arc_file_time_range[1].split('_')[1]))

        choose_tstamps_within_start_stop = all_tstamps_just_time[(all_tstamps_just_time >= this_time_range_start) & (all_tstamps_just_time <= this_time_range_stop)]
        choose_tstamps_within_start_stop = np.sort(choose_tstamps_within_start_stop)
        print(choose_tstamps_within_start_stop)

        for tstamp in choose_tstamps_within_start_stop:
            arc_file_list.append(dat_dir+str("%.6f" % tstamp).split('.')[0]+'_'+str("%.6f" % tstamp).split('.')[1]+'.dat')

    temperaturesHe10 = {}
    temperaturesHe10['0'] = []
    temperaturesHe10['1'] = []
    temperaturesHe10['2'] = []
    temperaturesHe10['3'] = []
    temperaturesHe10['4'] = []
    temperaturesHe10['5'] = []
    temperaturesHe10['6'] = []
    temperaturesHe10['7'] = []
    temperaturesHe10['8'] = []
    temperaturesHe10['9'] = []
    temperaturesHe10['10'] = []
    temperaturesHe10['11'] = []
    temperaturesHe10['12'] = []
    temperaturesHe10['13'] = []
    temperaturesHe10['14'] = []
    temperaturesHe10['15'] = []
    temperaturesRx = {}
    temperaturesRx['0'] = []
    temperaturesRx['1'] = []
    temperaturesRx['2'] = []
    temperaturesRx['3'] = []
    temperaturesRx['4'] = []
    temperaturesRx['6'] = []
    temperaturesRx['7'] = []
    temperaturesRx['8'] = []
    temperaturesRx['9'] = []
    temperaturesRx['10'] = []
    temperaturesRx['11'] = []
    temperaturesRx['12'] = []
    temperaturesRx['13'] = []
    temperaturesRx['14'] = []
    temperaturesOp = {}
    temperaturesOp['0'] = []
    temperaturesOp['1'] = []
    temperaturesOp['2'] = []
    temperaturesOp['3'] = []
    temperaturesOp['4'] = []
    temperaturesOp['5'] = []
    temperaturesOp['6'] = []
    temperaturesOp['7'] = []
    temperaturesOp['8'] = []
    temperaturesOp['9'] = []
    temperaturesOp['10'] = []
    temperaturesOp['11'] = []
    temperaturesOp['12'] = []
    temperaturesOp['13'] = []
    temperaturesOp['14'] = []
    temperaturesOp['15'] = []
    utc_time_cryo = []
    utc_time_tracker = []
    utc_time_cryo_date_time = []
    utc_time_tracker_date_time = []
    az = {}
    az['actual'] = []
    az['expected'] = []
    el = {}
    el['actual'] = []
    el['expected'] = []

    motor_current = {}
    motor_current['0'] = []
    motor_current['1'] = []
    motor_current['2'] = []
    motor_current['3'] = []
    bench = {}
    bench['actual'] = {}
    bench['actual']['0'] = []
    bench['actual']['1'] = []
    bench['actual']['2'] = []
    bench['actual']['3'] = []
    bench['actual']['4'] = []
    bench['actual']['5'] = []
    bench['expected'] = {}
    bench['expected']['0'] = []
    bench['expected']['1'] = []
    bench['expected']['2'] = []
    bench['expected']['3'] = []
    bench['expected']['4'] = []
    bench['expected']['5'] = []
    utc_time_bench = []
    dc_currents = {}
    dc_currents['0'] = []
    dc_currents['1'] = []
    dc_currents['2'] = []
    dc_currents['3'] = []
    dc_currents['4'] = []
    dc_currents['5'] = []
    dc_currents['6'] = []
    dc_currents['7'] = []
    tracker_time_status = []
    acu_status_byte_int = []
    acu_status_byte = {}
    acu_status_byte['0'] = []
    acu_status_byte['1'] = []
    acu_status_byte['2'] = []
    acu_status_byte['3'] = []
    acu_status_byte['4'] = []
    acu_status_byte['5'] = []
    acu_status_byte['6'] = []
    acu_status_byte['7'] = []
    cabin_temperature = [] #antenna,scu,temp[20]

    weather_station = {}
    weather_station['utc'] = []
    weather_station['utc_date_time'] = []
    weather_station['windSpeed'] = []
    weather_station['windDirection'] = []
    weather_station['airTemperature'] = []
    
    for aa in arc_file_list:
        print('loading ' + aa)
        af = gcp.ARCFileReader(aa)

        count_frames = 0
        while count_frames < 1000.:
            try:
                fr = af(None)[0]
                for key in temperaturesHe10:
                    temperaturesHe10[key].append(fr['array']['cryo']['temperature'][0][int(key)])
                for key in temperaturesOp:
                    temperaturesOp[key].append(fr['array']['cryo']['temperature'][1][int(key)])
                for key in temperaturesRx:
                    temperaturesRx[key].append(fr['array']['cryo']['temperature'][2][int(key)])

                bits = np.zeros(8)
                this_acu_status_byte = fr['antenna0']['acu']['acu_status'].value
                acu_status_byte_int.append(this_acu_status_byte)
                for i in range(8):
                    bits[i] = this_acu_status_byte & 1
                    this_acu_status_byte = this_acu_status_byte >> 1
                for bb in range(len(bits)):
                    acu_status_byte[str(bb)].append(bits[bb])

                cabin_temperature.append(fr['antenna0']['scu']['temp'][20])
                this_utc_time_cryo = fr['array']['cryo']['utc']
                utc_time_cryo.append(this_utc_time_cryo.mjd)
                utc_time_cryo_date_time.append(str(this_utc_time_cryo))
                
                weather_station['utc'].append(fr['array']['weather']['utc'].mjd)
                weather_station['utc_date_time'].append(str(fr['array']['weather']['utc']))
                weather_station['windSpeed'].append(fr['array']['weather']['windSpeed'].value)
                weather_station['windDirection'].append(fr['array']['weather']['windDirection'].value)
                weather_station['airTemperature'].append(fr['array']['weather']['airTemperature'].value)

                this_utc_time_tracker = fr['antenna0']['tracker']['utc'][0]
                this_az_actual = fr['antenna0']['tracker']['actual'][0]
                this_az_expected = fr['antenna0']['tracker']['expected'][0]
                this_el_actual = fr['antenna0']['tracker']['actual'][1]
                this_el_expected = fr['antenna0']['tracker']['expected'][1]
                this_tracker_time_status = fr['antenna0']['tracker']['time_status'][0]
                this_utc_time_bench = fr['antenna0']['scu']['benchSampleTime'][0]
                for key in dc_currents:
                    for ii in np.linspace(0,99,100):
                        dc_currents[key].append(fr['array']['dc']['currents'][int(key)][int(ii)])
            
                for ii in np.linspace(0,99,100):
                    utc_time_tracker.append(this_utc_time_tracker[int(ii)].mjd)
                    utc_time_tracker_date_time.append(str(this_utc_time_tracker[int(ii)]))
                    az['actual'].append(this_az_actual[int(ii)])
                    az['expected'].append(this_az_expected[int(ii)])
                    el['actual'].append(this_el_actual[int(ii)])
                    el['expected'].append(this_el_expected[int(ii)])
                    tracker_time_status.append(this_tracker_time_status[int(ii)])
                    for key in motor_current:
                        motor_current[key].append(fr['antenna0']['tracker']['motor_current'][int(key)][int(ii)])
                    utc_time_bench.append(this_utc_time_bench[int(ii)].mjd)
                    for key in bench['actual']:
                        bench['actual'][key].append(fr['antenna0']['scu']['benchActual'][int(key)][int(ii)])
                        bench['expected'][key].append(fr['antenna0']['scu']['benchExpected'][int(key)][int(ii)])
                count_frames += 1.
            except IndexError:
                #print 'this file has fewer than 1000 frames'
                count_frames += 1.

    utc_time_cryo = np.array(utc_time_cryo)*24.*3600.
    utc_time_tracker = np.array(utc_time_tracker)*24.*3600.
    utc_time_weather = np.array(weather_station['utc'])*24.*3600.
    #utc_time_bench = np.array(utc_time_bench)*24*3600.
    utc_time_cryo_shift = utc_time_cryo - utc_time_tracker[0]
    utc_time_tracker_shift = utc_time_tracker - utc_time_tracker[0]
    utc_time_weather_shift = utc_time_weather - utc_time_tracker[0]

    #make plot times for plot_date:
    utc_tracker_date_time = []
    for dt in utc_time_tracker_date_time:
        dt_start,dt_end=dt.split('.')
        dt_end = dt_end[0:6]#need to truncate seconds decimal places
        this_dt = datetime.datetime.strptime(dt_start,"%d-%b-%Y:%H:%M:%S")
        this_dt_out = this_dt.replace(microsecond = int(dt_end))
        utc_tracker_date_time.append(this_dt_out)
    tracker_dates = mpl.dates.date2num(utc_tracker_date_time)

    utc_cryo_date_time = []
    for dt in utc_time_cryo_date_time:
        dt_start,dt_end=dt.split('.')
        dt_end = dt_end[0:6]#need to truncate seconds decimal places
        this_dt = datetime.datetime.strptime(dt_start,"%d-%b-%Y:%H:%M:%S")
        this_dt_out = this_dt.replace(microsecond = int(dt_end))
        utc_cryo_date_time.append(this_dt_out)
    cryo_dates = mpl.dates.date2num(utc_cryo_date_time)

    el['actual'] = np.array(el['actual'])*180./np.pi
    el['expected'] = np.array(el['expected'])*180./np.pi
    az['actual'] = np.array(az['actual'])*180./np.pi
    az['expected'] = np.array(az['expected'])*180./np.pi

    az_vel_actual = np.diff(az['actual'])/np.diff(utc_time_tracker)
    el_vel_actual = np.diff(el['actual'])/np.diff(utc_time_tracker)

    az_accel_actual = np.diff(az_vel_actual)/np.diff(utc_time_tracker[0:-1])
    el_accel_actual = np.diff(el_vel_actual)/np.diff(utc_time_tracker[0:-1])

    az_accel_actual_sm = gf1d(az_accel_actual,5)
    el_accel_actual_sm = gf1d(el_accel_actual,5)

    data_out = {}
    data_out['time_cryo'] = utc_time_cryo_shift
    data_out['utc_cryo_date_time'] = cryo_dates
    data_out['time_tracker'] = utc_time_tracker_shift
    data_out['utc_tracker_date_time'] = tracker_dates
    data_out['tracker_time_status'] = tracker_time_status
    data_out['az'] = az
    data_out['el'] = el
    data_out['bench'] = bench
    data_out['bench_time'] = utc_time_bench
    data_out['dc_currents'] = dc_currents
    data_out['cabin_temperature'] = cabin_temperature
    data_out['temperaturesHe10'] = temperaturesHe10
    data_out['temperaturesRx'] = temperaturesRx
    data_out['temperaturesOp'] = temperaturesOp
    data_out['time_weather'] = utc_time_weather_shift
    data_out['weather_station'] = weather_station
    data_out['acu_status_byte'] = acu_status_byte
    data_out['acu_status_byte_int'] = acu_status_byte_int
    data_out['velocity'] = {}
    data_out['velocity']['az'] = az_vel_actual
    data_out['velocity']['el'] = el_vel_actual
    data_out['accel'] = {}
    data_out['accel']['az'] = az_accel_actual
    data_out['accel']['az_smooth'] = az_accel_actual_sm
    data_out['accel']['el'] = el_accel_actual
    data_out['accel']['el_smooth'] = el_accel_actual_sm

    return data_out

def plot_temps_az_el(data_out,plot_label='',time_range=[],az_range=[],el_range=[],temp_range=[],azeldiff_range=[],vel_range=[],accel_range=[],smooth_dc_currents=False,save_dir=''):
    
    # Generates temperatures and telescope movement png from data loaded into data_out.  
    # Optional:
    # time_range = [] UTC time range over which to plot 
    #      (format: time_range = ['02/06/2018 01:00:00','02/07/2018 02:00:00'])
    # temp_range = [] temperature range (K) for the UC/IC temperatures (subplot 1)
    # other *_range = [] plot range for az, el, vel, etc.

    utc_time_cryo = data_out['time_cryo']
    utc_time_tracker = data_out['time_tracker']
    utc_time_weather = data_out['time_weather']
    utc_tracker_date_time = data_out['utc_tracker_date_time']
    utc_cryo_date_time = data_out['utc_cryo_date_time']
    tracker_time_status = data_out['tracker_time_status']
    temperaturesHe10 = data_out['temperaturesHe10']
    temperaturesRx = data_out['temperaturesRx']
    az = data_out['az']
    el = data_out['el']
    bench = data_out['bench']
    dc_currents = data_out['dc_currents']
    cabin_temperature = data_out['cabin_temperature']
    weather_station = data_out['weather_station']
    acu_status_byte = data_out['acu_status_byte']
    az_vel_actual = data_out['velocity']['az']
    el_vel_actual = data_out['velocity']['el']
    az_accel_actual = data_out['accel']['az']
    el_accel_actual = data_out['accel']['el']
    az_accel_actual_sm = data_out['accel']['az_smooth']
    el_accel_actual_sm = data_out['accel']['el_smooth']

    if time_range:
        start_time = datetime.datetime.strptime(time_range[0],'%m/%d/%Y %H:%M:%S')
        start_time = mpl.dates.date2num(start_time)
        stop_time = datetime.datetime.strptime(time_range[1],'%m/%d/%Y %H:%M:%S')
        stop_time = mpl.dates.date2num(stop_time)
    
    if el_range == []:
        el_range = [0,90.]
    if az_range == []:
        az_range = [0,360.]

    pl.figure('Plot ALL the things',figsize=(18,25))
    pl.subplot(12,1,1)
    pl.title('Compare fridge temps with az/el for '+plot_label)
    # UC/IC temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['11'],color='turquoise',label = 'LC tower',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['10'],color='indigo',label = 'UC stage',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['0'],color='indigo',ls='--',label = 'UC head',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['12'],color='tomato',label = 'IC stage',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['1'],color='tomato',ls='--',label = 'IC head',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    if temp_range:
        pl.ylim(temp_range)
    else:
        pl.ylim([0.25,0.6])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,2)
    # selected He4 and 4K temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['13'],color='chartreuse',label = '4K head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['14'],color='navy',label = '4K SQ strap',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['4'],color='darkorange',label='4K pl near PTC',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['2'],color='turquoise',label = 'He4 head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['3'],color='violet',label = 'He4 FB',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,6.])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,3)
    # Cabin temperature
    pl.plot_date(utc_cryo_date_time,cabin_temperature,color='purple',label = 'cabin',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [C]',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])
     
    pl.subplot(12,1,4)
    #az position
    pl.plot_date(utc_tracker_date_time,az['expected'],color='k',ls='--',label='expected',marker='')
    pl.plot_date(utc_tracker_date_time,az['actual'],color='tomato',zorder=5,label='actual',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('az [deg]',fontsize=10)
    if az_range:
        pl.ylim([az_range[0],az_range[1]])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,5)
    # el position
    pl.plot_date(utc_tracker_date_time,el['expected'],color='k',ls='--',label = 'expected',marker='')
    pl.plot_date(utc_tracker_date_time,el['actual'],color='cyan',zorder=5,label = 'actual',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('el [deg]',fontsize=10)
    if el_range:
        pl.ylim([el_range[0],el_range[1]])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,6)
    # az/el actual - expected
    pl.plot_date(utc_tracker_date_time,az['actual']-az['expected'],color='tomato',label='az',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,el['actual']-el['expected'],color='cyan',label='el',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('act-exp [deg]',fontsize=10)
    if azeldiff_range:
        pl.ylim(azeldiff_range)
    else:
        pl.ylim([-5,5])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,7)
    # az/el velocity
    pl.plot_date(utc_tracker_date_time[0:-1],az_vel_actual,color='darkorange',label='az velocity',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[0:-1],el_vel_actual,color='blue',label='el velocity',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    if vel_range:
        pl.ylim(vel_range)
    else:
        pl.ylim([-3,3])
    pl.ylabel('[deg/s]',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,8)
    # az/el acceleration
    pl.plot_date(utc_tracker_date_time[0:-2],az_accel_actual,color='gold',label='az accel',alpha=0.3,ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[0:-2],el_accel_actual,color='cyan',label='el accel',alpha=0.3,ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[0:-2],az_accel_actual_sm,color='red',label='az accel smooth',zorder=3,ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[0:-2],el_accel_actual_sm,color='navy',label='el accel smooth',zorder=3,ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('[deg/s$^2$]',fontsize=10)
    if accel_range:
        pl.ylim(accel_range)
    else:
        pl.ylim([-1.5,1.5])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,9)
    # motor currents
    if smooth_dc_currents:
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['0'])-1000.)*(1./100.),5),color='tomato',label = 'motor voltage El1 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['1'])-1000.)*(1./100.),5),color='darkorange',label = 'motor voltage El2 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['2'])-1000.)*(1./100.),5),color='gold',label = 'motor voltage El3 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['3'])-1000.)*(1./100.),5),color='lime',label = 'motor voltage El4 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['4'])-1000.)*(1./100.),5),color='turquoise',label = 'motor voltage Az1 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['5'])-1000.)*(1./100.),5),color='teal',label = 'motor voltage Az2 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['6'])-1000.)*(1./100.),5),color='violet',label = 'motor voltage Az3 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,gf1d((np.array(dc_currents['7'])-1000.)*(1./100.),5),color='indigo',label = 'motor voltage Az4 smooth',alpha=0.6,ls='-',marker='')
        pl.plot_date([utc_tracker_date_time[0],utc_tracker_date_time[-1]],[0,0],ls=':',color='grey',alpha=0.6,marker='')
    else:
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['0'])-1000.)*(1./100.),color='tomato',label = 'motor voltage El1',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['1'])-1000.)*(1./100.),color='darkorange',label = 'motor voltage El2',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['2'])-1000.)*(1./100.),color='gold',label = 'motor voltage El3',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['3'])-1000.)*(1./100.),color='lime',label = 'motor voltage El4',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['4'])-1000.)*(1./100.),color='turquoise',label = 'motor voltage Az1',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['5'])-1000.)*(1./100.),color='teal',label = 'motor voltage Az2',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['6'])-1000.)*(1./100.),color='violet',label = 'motor voltage Az3',alpha=0.6,ls='-',marker='')
        pl.plot_date(utc_tracker_date_time,(np.array(dc_currents['7'])-1000.)*(1./100.),color='indigo',label = 'motor voltage Az4',alpha=0.6,ls='-',marker='')
        pl.plot_date([utc_tracker_date_time[0],utc_tracker_date_time[-1]],[0,0],ls=':',color='grey',alpha=0.6,marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('[V]',fontsize=10)
    #pl.ylim([800,1300])
    pl.ylim([-2.5,3.5])
    #pl.ylim([-4,5])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,10)
    # bench location (actual)
    pl.plot_date(utc_tracker_date_time,bench['actual']['0'],color='tomato',label = 'bench 0 act',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,bench['actual']['1'],color='gold',label = 'bench 1 act',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,bench['actual']['2'],color='lime',label = 'bench 2 act',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,bench['actual']['3'],color='turquoise',label = 'bench 3 act',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,bench['actual']['4'],color='violet',label = 'bench 4 act',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,bench['actual']['5'],color='indigo',label = 'bench 5 act',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.ylim([30000,50000])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('[um]',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,11)
    # faults: tracker time diff, tracker status, ACU status bits
    pl.plot_date(utc_tracker_date_time[0:-1],np.diff(utc_time_tracker),label='tracker time diff',color='indigo',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time,tracker_time_status,label='tracker time status',color='lime',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['0'])+10,label = 'El motor fault+10',color='tomato',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['1'])+12,label = 'AZ motor fault+12',color='darkorange',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['2'])+14,label = 'TB inactive+14',color='gold',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['3'])+16,label = 'TB par out of range+16',color='green',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['4'])+18,label = 'TB Xmit fault+18',color='turquoise',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['5'])+20,label = 'Tel on target+20',color='teal',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['6'])+22,label = 'El enabled+22',color='navy',ls='-',marker='')
    pl.plot_date(utc_tracker_date_time[::100],np.array(acu_status_byte['7'])+24,label = 'AZ enabled+24',color='violet',ls='-',marker='')
    pl.gca().axes.get_xaxis().set_ticklabels([])
    pl.grid()
    pl.ylim([-2,26])
    pl.legend(loc=2,prop={'size':7})
    pl.ylabel('faults',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(12,1,12)
    # weather data
    #this is slightly wrong (by <= 1sec) -- need to interpolate weather time and data onto cryo time
    pl.plot_date(utc_cryo_date_time,weather_station['windSpeed'],label='wind speed [m/s]',color='darkorange',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,weather_station['windDirection'],label='wind direction [deg]',color='magenta',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,weather_station['airTemperature'],label='air temp [C]',color='turquoise',ls='-',marker='')
    pl.gca().xaxis.set_major_formatter(mpl.dates.DateFormatter('%m/%d/%Y %H:%M:%S'))
    pl.xticks(rotation=90)
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('weather',fontsize=10)
    pl.ylim([min(weather_station['airTemperature'])*1.2,360.])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplots_adjust(hspace=0.05)
    pl.savefig(save_dir+'Plot_temps_az_el_'+plot_label+'.png',bbox_inches='tight')

    pl.close('all')
    return
            
def analyze_dc_current(data_out,time_range=[],motor_number=''):
    # not finished...
    time_tracker = data_out['time_tracker']
    this_dc_current = np.array(data_out['dc_currents'][motor_number])
    this_dc_current_timerange = this_dc_current[(time_tracker >= time_range[0]) & (time_tracker <= time_range[1])]
    zero_pad_dc_current_init = np.insert(this_dc_current,0,np.zeros(len(this_dc_current)/2.))
    zero_pad_dc_current = np.append(zero_pad_dc_current_init,np.zeros(len(this_dc_current)/2.))
    dc_current_fft = np.fft.fft(zero_pad_dc_current)
    time_step = 0.01
                                 
    freqs=np.fft.fftfreq(len(zero_pad_dc_current),d=time_step)    
    return dc_current_fft,freqs


def plot_temps_only(data_out,plot_label='',time_range=[],temp_range=[],save_dir=''):
    
    # Generates temperatures-only png from data loaded into data_out.  
    # Optional:
    # time_range = [] UTC time range over which to plot 
    #      (format: time_range = ['02/06/2018 01:00:00','02/07/2018 02:00:00'])
    # temp_range = [] temperature range (K) for the UC/IC temperatures (subplot 1)

    utc_cryo_date_time = data_out['utc_cryo_date_time']
    temperaturesHe10 = data_out['temperaturesHe10']
    cabin_temperature = data_out['cabin_temperature']
    temperaturesRx = data_out['temperaturesRx']
    temperaturesOp = data_out['temperaturesOp']

    if time_range:
        start_time = datetime.datetime.strptime(time_range[0],'%m/%d/%Y %H:%M:%S')
        start_time = mpl.dates.date2num(start_time)
        stop_time = datetime.datetime.strptime(time_range[1],'%m/%d/%Y %H:%M:%S')
        stop_time = mpl.dates.date2num(stop_time)

    pl.figure('Plot ALL the temps',figsize=(18,25))

    pl.subplot(11,1,1)
    pl.title('Plot all temperatures for '+plot_label)
    # UC/IC temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['11'],color='turquoise',label = 'LC tower',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['10'],color='indigo',label = 'UC stage',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['0'],color='indigo',ls='--',label = 'UC head',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['12'],color='tomato',label = 'IC stage',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['1'],color='tomato',ls='--',label = 'IC head',marker='')
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.grid()
    if temp_range:
        pl.ylim(temp_range)
    else:
        pl.ylim([0,5])
    if time_range:
        pl.xlim([start_time,stop_time])
        
    pl.subplot(11,1,2)
    # He4 temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['2'],color='turquoise',label = 'He4 head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['3'],color='violet',label = 'He4 FB',ls='-',marker='')
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.grid()
    pl.ylim([0,6.])
    #pl.ylim([0,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,3)
    # He10 fridge pump temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['6'],color='indigo',label = 'UC pump',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['5'],color='tomato',label = 'IC pump',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['4'],color='gold',label = 'He4 pump',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,60])
    #pl.ylim([0,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,4)
    #He10 fridge switch temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['9'],color='indigo',label = 'UC switch',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['8'],color='tomato',label = 'IC switch',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['7'],color='gold',label = 'He4 switch',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,30.])
    #pl.ylim([0,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,5)
    # He10 and Rx 4K temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['13'],color='chartreuse',label = '4K head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['14'],color='navy',label = '4K SQ strap',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['0'],color='purple',label='4K plate far',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['1'],color='teal',label='4K strap opt side',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['3'],color='gold',label='4K plate top',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['4'],color='darkorange',label='4K plate near PTC',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([2.5,6.])
    #pl.ylim([2,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,6)
    # SQ: WH temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesRx['2'],color='lime',label='SQ: WH1 slot 2',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['8'],color='turquoise',label='SQ: WH2 Slot 7',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['9'],color='magenta',label='SQ: WH5 Slot 2',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['10'],color='indigo',label='SQ: WH3 Slot 5',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([2.,18.])
    #pl.ylim([2,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,7)
    # He10 and Rx 50K temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesHe10['15'],color='turquoise',label = '50K head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['6'],color='violet',label = '50K harness3 middle',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['7'],color='navy',label='50K strap',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([30,60])
    #pl.ylim([30,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,8)
    # Optics cryostat 4K temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesOp['8'],color='turquoise',label = 'G1 4K head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['9'],color='violet',label = 'G2 4K strap',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['10'],color='navy',label='G3 4K lens tab',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['11'],color='chartreuse',label='G4 4K lens tab far',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['14'],color='gold',label='R3 4K Lyot flange',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['12'],color='tomato',label='R1 4K top top PTC',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['15'],color='darkorange',label='R4 4K lyot',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,15])
    #pl.ylim([0,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,9)
    # Optics cryostat 50K temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesOp['0'],color='gold',label = 'B1 50K WBP near',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['1'],color='darkorange',label = 'B2 50K WBP far',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['2'],color='tomato',label='B3 50K diving board',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['3'],color='purple',label='B4 50K top bot PTC',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['4'],color='navy',label='Y1 50K head',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['5'],color='turquoise',label='Y2 50K window strap near',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['6'],color='teal',label='Y3 50K tube strap near',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['7'],color='chartreuse',label='Y4 50K tube',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesOp['13'],color='green',label='R2 50K mid-op bot PTC',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,100])
    #pl.ylim([40,150])
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,10)
    # Calibrator temperatures
    pl.plot_date(utc_cryo_date_time,temperaturesRx['11'],color='gold',label = 'Cal fil TC',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['12'],color='turquoise',label = 'Cal ambient 1',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['13'],color='chartreuse',label = 'Cal ambient 2',ls='-',marker='')
    pl.plot_date(utc_cryo_date_time,temperaturesRx['14'],color='navy',label = 'Cal ambient 3',ls='-',marker='')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [K]',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplot(11,1,11)
    # Cabin temperature
    pl.plot_date(utc_cryo_date_time,cabin_temperature,color='purple',label = 'cabin',ls='-',marker='')
    pl.gca().xaxis.set_major_formatter(mpl.dates.DateFormatter('%m/%d/%Y %H:%M:%S'))
    pl.xticks(rotation=90)
    pl.grid()
    pl.legend(loc=2,prop={'size':9})
    pl.ylabel('Temps [C]',fontsize=10)
    if time_range:
        pl.xlim([start_time,stop_time])

    pl.subplots_adjust(hspace=0.05)

    pl.savefig(save_dir+'Plot_temps_only_'+plot_label+'.png',bbox_inches='tight')

    pl.close('all')
    return

def save_data_out(data_out,pkl_label='',save_dir=''):
    pkl.dump(data_out,open(save_dir+pkl_label+'_arc_data_out.pkl','wb'))
    return
