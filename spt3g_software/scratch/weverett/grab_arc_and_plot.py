from spt3g import core, gcp
import numpy as np
import pylab as pl
import glob as glob
import cPickle as pkl
from scipy.ndimage import gaussian_filter1d as gf1d

def grab_and_append_arc_frame_data_kludge_way(arc_file_list_txt='',arc_file_time_range=[],dat_dir='/spt_data/arc/',year='y2'):
    #usage, for example:
    #data_out = gaap.grab_and_append_arc_frame_data_kludge_way(arf_file_time_range=['20180306_000000','20180307_000000'],data_dir = 'path/to/arc/files')

    arc_file_list = []
    if arc_file_list_txt:
        arc_files = open(arc_file_list_txt,'r')
        arcs = arc_files.readlines()
        for aa in arcs:
            arc_file_list.append(aa.split('\n')[0])
    if arc_file_time_range:
        print 'getting arc files from time range'
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
        print choose_tstamps_within_start_stop

        for tstamp in choose_tstamps_within_start_stop:
            arc_file_list.append(dat_dir+str("%.6f" % tstamp).split('.')[0]+'_'+str("%.6f" % tstamp).split('.')[1]+'.dat')

    temperatures = {}
    temperatures['0'] = []
    temperatures['1'] = []
    temperatures['2'] = []
    temperatures['3'] = []
    temperatures['4'] = []
    temperatures['5'] = []
    temperatures['6'] = []
    temperatures['7'] = []
    temperatures['8'] = []
    temperatures['9'] = []
    temperatures['10'] = []
    temperatures['11'] = []
    temperatures['12'] = []
    temperatures['13'] = []
    temperatures['14'] = []
    temperatures['15'] = []
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
    weather_station['windSpeed'] = []
    weather_station['windDirection'] = []
    weather_station['airTemperature'] = []
    
    for aa in arc_file_list:
        print 'loading ', aa
        af = gcp.ARCFileReader(aa)

        count_frames = 0
        while count_frames < 1000.:
            try:
                fr = af(None)[0]
                for key in temperatures:
                    temperatures[key].append(fr['array']['cryo']['temperature'][0][int(key)])
                for key in temperaturesOp:
                    temperaturesOp[key].append(fr['array']['cryo']['temperature'][1][int(key)])
                for key in temperaturesRx:
                    temperaturesRx[key].append(fr['array']['cryo']['temperature'][2][int(key)])
                bits = np.zeros(8)
                this_acu_status_byte = fr['antenna0']['acu']['acu_status'].value
                for i in range(8):
                    bits[i] = this_acu_status_byte & 1
                    this_acu_status_byte = this_acu_status_byte >> 1
                for bb in range(len(bits)):
                    acu_status_byte[str(bb)].append(bits[bb])

                cabin_temperature.append(fr['antenna0']['scu']['temp'][20])
                this_utc_time_cryo = fr['array']['cryo']['utc']
                utc_time_cryo.append(this_utc_time_cryo.mjd)
                
                weather_station['utc'].append(fr['array']['weather']['utc'].mjd)
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
                print 'this file has fewer than 1000 frames'
                count_frames += 1.

    utc_time_cryo = np.array(utc_time_cryo)*24.*3600.
    utc_time_tracker = np.array(utc_time_tracker)*24.*3600.
    utc_time_weather = np.array(weather_station['utc'])*24.*3600.
    #utc_time_bench = np.array(utc_time_bench)*24*3600.
    utc_time_cryo_shift = utc_time_cryo - utc_time_tracker[0]
    utc_time_tracker_shift = utc_time_tracker - utc_time_tracker[0]
    utc_time_weather_shift = utc_time_weather - utc_time_tracker[0]

    el['actual'] = np.array(el['actual'])*180./np.pi
    el['expected'] = np.array(el['expected'])*180./np.pi
    az['actual'] = np.array(az['actual'])*180./np.pi
    az['expected'] = np.array(az['expected'])*180./np.pi

    data_out = {}
    data_out['temperatures'] = temperatures
    data_out['time_cryo'] = utc_time_cryo_shift
    data_out['time_tracker'] = utc_time_tracker_shift
    data_out['tracker_time_status'] = tracker_time_status
    data_out['az'] = az
    data_out['el'] = el
    data_out['bench'] = bench
    data_out['bench_time'] = utc_time_bench
    data_out['dc_currents'] = dc_currents
    data_out['cabin_temperature'] = cabin_temperature
    data_out['temperaturesRx'] = temperaturesRx
    data_out['temperaturesOp'] = temperaturesOp
    data_out['time_weather'] = utc_time_weather_shift
    data_out['weather_station'] = weather_station
    data_out['acu_status_byte'] = acu_status_byte

    #return temperatures,utc_time_cryo_shift,utc_time_tracker_shift,az,el,motor_current,bench,dc_currents,tracker_time_status
    return data_out

def plot_kludge(data_out,plot_label='',time_range=[],az_range=[],el_range=[],temp_range=[],azeldiff_range=[],vel_range=[],accel_range=[],smooth_dc_currents=False,save_pkl=False):
#def plot_kludge(temperatures,utc_time_cryo,utc_time_tracker,az,el,motor_current,bench,dc_currents,tracker_time_status,plot_label='',time_range=[],az_range=[],el_range=[],temp_range=[],azeldiff_range=[],vel_range=[],accel_range=[],smooth_dc_currents=False):

    utc_time_cryo = data_out['time_cryo']
    utc_time_tracker = data_out['time_tracker']
    utc_time_weather = data_out['time_weather']
    tracker_time_status = data_out['tracker_time_status']
    temperatures = data_out['temperatures']
    az = data_out['az']
    el = data_out['el']
    bench = data_out['bench']
    dc_currents = data_out['dc_currents']
    cabin_temperature = data_out['cabin_temperature']
    temperaturesRx = data_out['temperaturesRx']
    weather_station = data_out['weather_station']
    acu_status_byte = data_out['acu_status_byte']

    if time_range == []:
        time_range = [utc_time_tracker[0],utc_time_tracker[-1]]
    if el_range == []:
        el_range = [0,90]
    az_vel_actual = np.diff(az['actual'])/np.diff(utc_time_tracker)
    el_vel_actual = np.diff(el['actual'])/np.diff(utc_time_tracker)

    az_accel_actual = np.diff(az_vel_actual)/np.diff(utc_time_tracker[0:-1])
    el_accel_actual = np.diff(el_vel_actual)/np.diff(utc_time_tracker[0:-1])

    az_accel_actual_sm = gf1d(az_accel_actual,5)
    el_accel_actual_sm = gf1d(el_accel_actual,5)

    data_out['velocity'] = {}
    data_out['velocity']['az'] = az_vel_actual
    data_out['velocity']['el'] = el_vel_actual
    data_out['accel'] = {}
    data_out['accel']['az'] = az_accel_actual
    data_out['accel']['az_smooth'] = az_accel_actual_sm
    data_out['accel']['el'] = el_accel_actual
    data_out['accel']['el_smooth'] = el_accel_actual_sm

    pl.figure('Plot ALL the things',figsize=(12,18))
    pl.subplot(13,1,1)
    pl.title('Compare fridge temps with az/el for '+plot_label)
    pl.plot(utc_time_cryo,temperatures['11'],color='teal',label = 'LC tower')
    pl.plot(utc_time_cryo,temperatures['10'],color='indigo',label = 'UC stage')
    pl.plot(utc_time_cryo,temperatures['0'],color='indigo',linestyle='--',label = 'UC head')
    pl.plot(utc_time_cryo,temperatures['12'],color='tomato',label = 'IC stage')
    pl.plot(utc_time_cryo,temperatures['1'],color='tomato',linestyle='--',label = 'IC head')
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.grid()
    if temp_range:
        pl.ylim(temp_range)
    else:
        pl.ylim([0.25,0.6])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,2)
    pl.plot(utc_time_cryo,temperatures['13'],color='gold',label = '4K head')
    pl.plot(utc_time_cryo,temperatures['14'],color='navy',label = '4K SQ strap')
    pl.plot(utc_time_cryo,temperatures['6'],color='lime',label = 'UC pump')
    pl.plot(utc_time_cryo,temperatures['5'],color='olive',label = 'IC pump')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([2.5,6.])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,3)
    pl.plot(utc_time_cryo,cabin_temperature,color='coral',label = 'Cabin temp')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [C]',fontsize=10)
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,4)
    pl.plot(utc_time_cryo,temperatures['2'],color='turquoise',label = 'He4Head')
    pl.plot(utc_time_cryo,temperatures['3'],color='violet',label = 'He4FB')
    pl.plot(utc_time_cryo,temperaturesRx['1'],color='navy',label='4K strap Optics side')
    pl.plot(utc_time_cryo,temperaturesRx['4'],color='darkorange',label='4K plate near PTC')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.ylim([0,5])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])
     
    pl.subplot(13,1,5)
    pl.plot(utc_time_tracker,az['expected'],color='k',linestyle='--',label='expected')
    pl.plot(utc_time_tracker,az['actual'],color='tomato',zorder=5,label='actual')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('az [deg]')
    if az_range:
        pl.ylim([az_range[0],az_range[1]])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,6)
    pl.plot(utc_time_tracker,el['expected'],color='k',linestyle='--',label = 'expected')
    pl.plot(utc_time_tracker,el['actual'],color='cyan',zorder=5,label = 'actual')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('el [deg]')
    if el_range:
        pl.ylim([el_range[0],el_range[1]])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,7)
    pl.plot(utc_time_tracker,az['actual']-az['expected'],color='tomato',label='az')
    pl.plot(utc_time_tracker,el['actual']-el['expected'],color='cyan',label='el')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('act-exp [deg]')
    if azeldiff_range:
        pl.ylim(azeldiff_range)
    else:
        pl.ylim([-5,5])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,8)
    pl.plot(utc_time_tracker[0:-1],az_vel_actual,color='orange',label='az velocity')
    pl.plot(utc_time_tracker[0:-1],el_vel_actual,color='blue',label='el velocity')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    if vel_range:
        pl.ylim(vel_range)
    else:
        pl.ylim([-3,3])
    pl.ylabel('[deg/sec]')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,9)
    pl.plot(utc_time_tracker[0:-2],az_accel_actual,color='gold',label='az accel',alpha=0.3)
    pl.plot(utc_time_tracker[0:-2],el_accel_actual,color='cyan',label='el accel',alpha=0.3)
    pl.plot(utc_time_tracker[0:-2],az_accel_actual_sm,color='red',label='az accel smooth',zorder=3)
    pl.plot(utc_time_tracker[0:-2],el_accel_actual_sm,color='indigo',label='el accel smooth',zorder=3)
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('[deg/sec^2]')
    if accel_range:
        pl.ylim(accel_range)
    else:
        pl.ylim([-1.5,1.5])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,10)
    if smooth_dc_currents:
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['0'])-1000.)*(1./100.),5),color='tomato',label = 'motor voltage El1 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['1'])-1000.)*(1./100.),5),color='darkorange',label = 'motor voltage El2 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['2'])-1000.)*(1./100.),5),color='gold',label = 'motor voltage El3 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['3'])-1000.)*(1./100.),5),color='lime',label = 'motor voltage El4 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['4'])-1000.)*(1./100.),5),color='turquoise',label = 'motor voltage Az1 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['5'])-1000.)*(1./100.),5),color='teal',label = 'motor voltage Az2 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['6'])-1000.)*(1./100.),5),color='violet',label = 'motor voltage Az3 smooth',alpha=0.6)
        pl.plot(utc_time_tracker,gf1d((np.array(dc_currents['7'])-1000.)*(1./100.),5),color='indigo',label = 'motor voltage Az4 smooth',alpha=0.6)
        pl.plot([time_range[0],time_range[1]],[0,0],linestyle=':',color='grey',alpha=0.6)
    else:
        pl.plot(utc_time_tracker,(np.array(dc_currents['0'])-1000.)*(1./100.),color='tomato',label = 'motor voltage El1',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['1'])-1000.)*(1./100.),color='darkorange',label = 'motor voltage El2',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['2'])-1000.)*(1./100.),color='gold',label = 'motor voltage El3',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['3'])-1000.)*(1./100.),color='lime',label = 'motor voltage El4',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['4'])-1000.)*(1./100.),color='turquoise',label = 'motor voltage Az1',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['5'])-1000.)*(1./100.),color='teal',label = 'motor voltage Az2',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['6'])-1000.)*(1./100.),color='violet',label = 'motor voltage Az3',alpha=0.6)
        pl.plot(utc_time_tracker,(np.array(dc_currents['7'])-1000.)*(1./100.),color='indigo',label = 'motor voltage Az4',alpha=0.6)
        pl.plot([time_range[0],time_range[1]],[0,0],linestyle=':',color='grey',alpha=0.6)
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':6})
    pl.ylabel('Motor voltage[V]')
    #pl.ylim([800,1300])
    pl.ylim([-2.5,3.5])
    #pl.ylim([-4,5])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,11)
    #pl.plot(utc_time_bench,bench['actual']['0'],color='tomato',label = 'bench 0 actual')
    pl.plot(utc_time_tracker,bench['actual']['0'],color='tomato',label = 'bench 0 actual')
    #pl.plot(utc_time,bench['expected']['0'],color='tomato',linestyle='--')
    #pl.plot(utc_time_bench,bench['actual']['1'],color='gold',label = 'bench 1 actual')
    pl.plot(utc_time_tracker,bench['actual']['1'],color='gold',label = 'bench 1 actual')
    #pl.plot(utc_time_bench,bench['expected']['1'],color='gold',linestyle='--')
    #pl.plot(utc_time_bench,bench['actual']['2'],color='green',label = 'bench 2 actual')
    pl.plot(utc_time_tracker,bench['actual']['2'],color='green',label = 'bench 2 actual')
    #pl.plot(utc_time_bench,bench['expected']['2'],color='green',linestyle='--')
    #pl.plot(utc_time_bench,bench['actual']['3'],color='turquoise',label = 'bench 3 actual')
    pl.plot(utc_time_tracker,bench['actual']['3'],color='turquoise',label = 'bench 3 actual')
    #pl.plot(utc_time_bench,bench['expected']['3'],color='turquoise',linestyle='--')
    #pl.plot(utc_time_bench,bench['actual']['4'],color='violet',label = 'bench 4 actual')
    pl.plot(utc_time_tracker,bench['actual']['4'],color='violet',label = 'bench 4 actual')
    #pl.plot(utc_time_bench,bench['expected']['4'],color='violet',linestyle='--')
    #pl.plot(utc_time_bench,bench['actual']['5'],color='indigo',label = 'bench 5 actual')
    pl.plot(utc_time_tracker,bench['actual']['5'],color='indigo',label = 'bench 5 actual')
    #pl.plot(utc_time_bench,bench['expected']['5'],color='indigo',linestyle='--')
    pl.ylabel('Bench pos [um]')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.ylim([30000,50000])
    pl.legend(loc=2,prop={'size':8})
    #pl.xlabel('UTC time [sec]')
    #pl.xlabel('Time [sec]')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,12)
    pl.plot(utc_time_tracker[0:-1],np.diff(utc_time_tracker),label='tracker time diff',color='indigo')
    pl.plot(utc_time_tracker,tracker_time_status,label='tracker time status',color='lime')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['0'])+10,label = 'El motor fault',color='tomato')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['1'])+12,label = 'AZ motor fault',color='darkorange')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['2'])+14,label = 'TB inactive',color='gold')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['3'])+16,label = 'TB par out of range',color='green')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['4'])+18,label = 'TB Xmit fault',color='turquoise')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['5'])+20,label = 'Tel on target',color='teal')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['6'])+22,label = 'El enabled',color='navy')
    pl.plot(utc_time_tracker[::100],np.array(acu_status_byte['7'])+24,label = 'AZ enabled',color='violet')
    #pl.xlabel('Time [sec]')
    pl.ylim([-2,26])
    pl.legend(loc=2,prop={'size':5})
    pl.grid()
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(13,1,13)
    pl.plot(utc_time_weather,weather_station['windSpeed'],label='wind speed [m/s]',color='darkorange')
    pl.plot(utc_time_weather,weather_station['windDirection'],label='wind direction [deg]',color='magenta')
    pl.plot(utc_time_weather,weather_station['airTemperature'],label='air temp [C]',color='turquoise')
    pl.xlabel('Time [sec]')
    #pl.ylim([-5,10])
    pl.legend(loc=2,prop={'size':8})
    pl.grid()
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplots_adjust(hspace=0.05)
    pl.savefig('Plot_temps_and_az_el_'+plot_label+'.png',bbox_inches='tight')

    #pl.show()

    pl.close('all')

    if save_pkl == True:
        pkl.dump(data_out,str(time_range[0])+'_'+str(time_range[1])+'_arc_data_out.pkl')

    return
            
def analyze_dc_current(data_out,time_range=[],motor_number=''):
    time_tracker = data_out['time_tracker']
    this_dc_current = np.array(data_out['dc_currents'][motor_number])
    this_dc_current_timerange = this_dc_current[(time_tracker >= time_range[0]) & (time_tracker <= time_range[1])]
    zero_pad_dc_current_init = np.insert(this_dc_current,0,np.zeros(len(this_dc_current)/2.))
    zero_pad_dc_current = np.append(zero_pad_dc_current_init,np.zeros(len(this_dc_current)/2.))
    dc_current_fft = np.fft.fft(zero_pad_dc_current)
    time_step = 0.01
                                 
    freqs=np.fft.fftfreq(len(zero_pad_dc_current),d=time_step)


    
    return dc_current_fft,freqs

        
def plot_temps_only(data_out,plot_label='',time_range=[],temp_range=[],save_pkl=False):
    
    utc_time_cryo = data_out['time_cryo']
    utc_time_tracker = data_out['time_tracker']
    tracker_time_status = data_out['tracker_time_status']
    temperatures = data_out['temperatures']
    az = data_out['az']
    el = data_out['el']
    bench = data_out['bench']
    dc_currents = data_out['dc_currents']
    cabin_temperature = data_out['cabin_temperature']
    temperaturesRx = data_out['temperaturesRx']
    temperaturesOp = data_out['temperaturesOp']

    if time_range == []:
        time_range = [utc_time_tracker[0],utc_time_tracker[-1]]

    pl.figure('Plot ALL the temps',figsize=(12,18))
    pl.subplot(10,1,1)
    pl.title('Plot all temps for '+plot_label)
    pl.plot(utc_time_cryo,temperatures['11'],color='teal',label = 'LC tower')
    pl.plot(utc_time_cryo,temperatures['10'],color='indigo',label = 'UC stage')
    pl.plot(utc_time_cryo,temperatures['0'],color='indigo',linestyle='--',label = 'UC head')
    pl.plot(utc_time_cryo,temperatures['12'],color='tomato',label = 'IC stage')
    pl.plot(utc_time_cryo,temperatures['1'],color='tomato',linestyle='--',label = 'IC head')
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.grid()
    if temp_range:
        pl.ylim(temp_range)
    else:
        pl.ylim([0,4.5])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,2)
    pl.plot(utc_time_cryo,temperatures['2'],color='gold',label = 'He4 head')
    pl.plot(utc_time_cryo,temperatures['3'],color='chartreuse',label = 'He4FB')
    this_axes=pl.gca()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    pl.grid()
    #pl.ylim([0,6.])
    pl.ylim([0,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,3)
    pl.plot(utc_time_cryo,temperatures['6'],color='indigo',label = 'UC pump')
    pl.plot(utc_time_cryo,temperatures['5'],color='tomato',label = 'IC pump')
    pl.plot(utc_time_cryo,temperatures['4'],color='gold',label = 'He4 pump')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([0,60])
    pl.ylim([0,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,4)
    pl.plot(utc_time_cryo,temperatures['9'],color='indigo',label = 'UC switch')
    pl.plot(utc_time_cryo,temperatures['8'],color='tomato',label = 'IC switch')
    pl.plot(utc_time_cryo,temperatures['7'],color='gold',label = 'He4 switch')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([0,30.])
    pl.ylim([0,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,5)
    pl.plot(utc_time_cryo,temperatures['13'],color='chartreuse',label = '4K head')
    pl.plot(utc_time_cryo,temperatures['14'],color='navy',label = '4K SQ strap')
    pl.plot(utc_time_cryo,temperaturesRx['0'],color='purple',label='4K plate Far')
    pl.plot(utc_time_cryo,temperaturesRx['1'],color='violet',label='4K strap Optics side')
    pl.plot(utc_time_cryo,temperaturesRx['3'],color='orange',label='4K plate Top')
    pl.plot(utc_time_cryo,temperaturesRx['4'],color='darkorange',label='4K plate near PTC')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([2.5,7.])
    pl.ylim([0,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,6)
    pl.plot(utc_time_cryo,temperaturesRx['2'],color='olive',label='SQ:WH1 slot 2')
    pl.plot(utc_time_cryo,temperaturesRx['8'],color='turquoise',label='SQ:WH2 Slot 7')
    pl.plot(utc_time_cryo,temperaturesRx['9'],color='blue',label='SQ:WH5 Slot 2')
    pl.plot(utc_time_cryo,temperaturesRx['10'],color='teal',label='SQ:WH3 Slot 5')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([2.,18.])
    pl.ylim([0,150])
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,7)
    pl.plot(utc_time_cryo,temperatures['15'],color='turquoise',label = '50Khead')
    pl.plot(utc_time_cryo,temperaturesRx['6'],color='violet',label = '50K Harness3 Middle')
    pl.plot(utc_time_cryo,temperaturesRx['7'],color='navy',label='50K strap')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([30,60])
    pl.ylim([30,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,8)
    pl.plot(utc_time_cryo,temperaturesOp['8'],color='turquoise',label = 'G1 4K head')
    pl.plot(utc_time_cryo,temperaturesOp['9'],color='violet',label = 'G2 4K strap')
    pl.plot(utc_time_cryo,temperaturesOp['10'],color='navy',label='G3 4K lens tab')
    pl.plot(utc_time_cryo,temperaturesOp['11'],color='chartreuse',label='G4 4K lens tab far')
    pl.plot(utc_time_cryo,temperaturesOp['14'],color='gold',label='R3 4K Lyot flange')
    pl.plot(utc_time_cryo,temperaturesOp['12'],color='tomato',label='R1 4K Top Top PTC')
    pl.plot(utc_time_cryo,temperaturesOp['15'],color='darkorange',label='R4 4K Lyot')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [K]',fontsize=10)
    #if temp_range:
    #    pl.ylim(temp_range)
    #else:
    #pl.ylim([3,6])
    pl.ylim([0,150])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,9)
    pl.plot(utc_time_cryo,temperaturesOp['0'],color='gold',label = 'B1 50K WBP near')
    pl.plot(utc_time_cryo,temperaturesOp['1'],color='darkorange',label = 'B2 50K WBP far')
    pl.plot(utc_time_cryo,temperaturesOp['2'],color='red',label='B3 50K Diving board')
    pl.plot(utc_time_cryo,temperaturesOp['3'],color='purple',label='B4 50K Top Bot PTC')
    pl.plot(utc_time_cryo,temperaturesOp['4'],color='navy',label='Y1 50K Head')
    pl.plot(utc_time_cryo,temperaturesOp['5'],color='turquoise',label='Y2 50K Window strap near')
    pl.plot(utc_time_cryo,temperaturesOp['6'],color='teal',label='Y3 50K Tube strap near')
    pl.plot(utc_time_cryo,temperaturesOp['7'],color='chartreuse',label='Y4 50K Tube')
    pl.plot(utc_time_cryo,temperaturesOp['13'],color='green',label='R2 50K Mid-op Bot PTC')
    this_axes=pl.gca()
    pl.grid()
    this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':7})
    pl.ylabel('Temps [K]',fontsize=10)
    #pl.ylim([35,65])
    pl.ylim([35,200])
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])

    pl.subplot(10,1,10)
    pl.plot(utc_time_cryo,cabin_temperature,color='coral',label = 'Cabin temp')
    this_axes=pl.gca()
    pl.grid()
    #this_axes.axes.get_xaxis().set_ticklabels([])
    pl.legend(loc=2,prop={'size':8})
    pl.ylabel('Temps [C]',fontsize=10)
    #pl.yscale('log')
    if time_range:
        pl.xlim([time_range[0],time_range[1]])
    pl.xlabel('time [sec] after start of first arc file loaded')
    pl.subplots_adjust(hspace=0.05)
    pl.savefig('Plot_all_temps_'+plot_label+'.png',bbox_inches='tight')

    #pl.show()

    pl.close('all')

    if save_pkl == True:
        pkl.dump(data_out,str(time_range[0])+'_'+str(time_range[1])+'_arc_data_out_temps_only.pkl')

    return
