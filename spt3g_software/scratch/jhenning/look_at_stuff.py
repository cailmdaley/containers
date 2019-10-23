from spt3g import core, gcp, std_processing
import numpy as np
import pylab as py
from scipy import ndimage

def grabScanRate(filename, sigma=20):
    d = core.G3File(filename)

    az_actual = []
    el_actual = []
    temperatures = []
    currents = []
    utc = []


    for f in d:
        utc.append(f['antenna0']['tracker']['utc'][0])
        az_actual.append(f['antenna0']['tracker']['actual'][0])
        el_actual.append(f['antenna0']['tracker']['actual'][1])
        temperatures.append(f['array']['cryo']['temperature'][0])#[0])
        currents.append(np.array(f['array']['dc']['currents']))

    utc = np.array(utc).flatten()
    utc = np.array([t.mjd for t in utc])*24*3600
    az_actual = np.array(az_actual).flatten()
    el_actual = np.array(el_actual).flatten()
    temperature = np.array(temperatures).flatten()

    smooth_az_actual = ndimage.gaussian_filter(az_actual, sigma)*180/np.pi
    smooth_el_actual = ndimage.gaussian_filter(el_actual, sigma)*180/np.pi

    az_rate = (smooth_az_actual[1:] - smooth_az_actual[:-1])* 100
    el_rate = (smooth_el_actual[1:] - smooth_el_actual[:-1])* 100
    az_accel = (az_rate[1:] - az_rate[:-1])*100
    el_accel = (el_rate[1:] - el_rate[:-1])*100

    return smooth_az_actual, smooth_el_actual, az_rate, el_rate, az_accel, el_accel, currents, np.array(temperatures), utc


def parseTemperatureLogs(tempfile='log_2017-02-01T14h34.dat'):
    d = open(tempfile,'r').read().split('\n')[:-1]

    time = []
    uc_head = []

    for i in range(2,len(d)):
        date = filter(None, d[i].split(' '))[0].split('/')
        time.append(core.G3Time(date[0]+'-'+date[1]+'-'+date[2]+':'+filter(None, d[i].split(' '))[1]).mjd*3600*24)
        uc_head.append(float(filter(None,d[i].split(' '))[2]))
        
    time = np.array(time)
    #time -= time[0]

    uc_head = np.array(uc_head)

    return time, uc_head
    
