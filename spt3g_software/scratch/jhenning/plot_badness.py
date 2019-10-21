import numpy as np
import pylab as py
from spt3g import core, gcp
import look_at_stuff as ut


#test_files = ['saturn_fast_point.g3',
#              'saturn_vfast_point.g3']

test_files = ['rcw38_fast_point_20170202_173400.g3',
              'saturn_fast_point_20170202_140630.g3',
              'saturn_focus_20170202_154911.g3']

              #'saturn_focus_20170203.g3']#,
              #'saturn_fast_point_20170203.g3']

#temp_time, uc_head = ut.parseTemperatureLogs(tempfile='log_2017-02-01T14h25.dat')

sigma = 5

#Plot acceleration versus stage temperature
for i in range(len(test_files)):
    az_actual, el_actual, az_rate, el_rate, az_accel, el_accel, motor, temperatures, utc = ut.grabScanRate(filename=test_files[i], sigma=sigma)

    fig = py.figure(i)
    fig.clf()
    ax1 = fig.add_axes((0.1,0.1,0.8,0.8))
    ax1.plot(utc[1:-1]-utc[0], az_accel, 'k-', label='Acceleration')
    #ax1.plot(utc[1:]-utc[0], az_rate, 'b--', label='Velocity')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Acceleration [Deg/s^2] / Velocity [Deg/s]')
    ax1.set_xlim((0,len(utc)/100))
    ax1.set_ylim((-2.5,2.5))
    py.legend(loc='lower right')
    ax2 = ax1.twinx()
    ax2.plot(utc[::100]-utc[0], temperatures[:,0])#, label='UC Head')
    ax2.plot(utc[::100]-utc[0], temperatures[:,10])#, label='UC Stage')
    ax2.plot(utc[::100]-utc[0], temperatures[:,1])#, label='IC Head')
    ax2.plot(utc[::100]-utc[0], temperatures[:,12])#, label='IC Stage')
    ax2.set_ylabel('Temperature [K]')
    ax2.set_xlim((0,len(utc)/100))
    ax2.set_ylim((0.280,0.345))
    py.legend(loc='upper right')
    py.title(test_files[i].split('.')[0]+' - '+str(sigma)+' Sigma Smoothing')
    py.savefig('az_'+test_files[i].split('.')[0]+'_'+str(sigma)+'sigma.png')

    fig.clf()
    ax1 = fig.add_axes((0.1,0.1,0.8,0.8))
    ax1.plot(utc[1:-1]-utc[0], el_accel, 'k-', label='Acceleration')
    #ax1.plot(utc[1:]-utc[0], el_rate, 'b--', label='Velocity')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Acceleration [Deg/s^2] / Velocity [Deg/s]')
    ax1.set_xlim((0,len(utc)/100))
    ax1.set_ylim((-2.5,2.5))     
    py.legend(loc='lower right')
    ax2 = ax1.twinx()
    #ax2.plot(utc[::100]-utc[0], temperatures, 'r-', label='Temperature')
    ax2.plot(utc[::100]-utc[0], temperatures[:,0])#, label='UC Head')
    ax2.plot(utc[::100]-utc[0], temperatures[:,10])#, label='UC Stage')
    ax2.plot(utc[::100]-utc[0], temperatures[:,1])#, label='IC Head')
    ax2.plot(utc[::100]-utc[0], temperatures[:,12])#, label='IC Stage')
    #ax2.plot(temp_time-13*3600 - utc[0], uc_head, 'r-', label='Temperature')
    ax2.set_ylabel('Temperature [K]')
    ax2.set_xlim((0,len(utc)/100))
    ax2.set_ylim((0.280,0.345))
    py.legend(loc='upper right')
    py.title(test_files[i].split('.')[0]+' - '+str(sigma)+' Sigma Smoothing')
    py.savefig('el_'+test_files[i].split('.')[0]+'_'+str(sigma)+'sigma.png')



    #fig = py.figure(i)
    #fig.clf()
    #ax1 = fig.add_axes((0.1,0.1,0.8,0.8))
    #ax1.plot(utc -utc[0], az_actual, 'k-', label='Az Actual')
    #ax1.set_xlabel('Time [s]')
    #ax1.set_ylabel('Azimuth [Deg]')
    #ax1.set_xlim((0,200))
    #ax1.set_ylim((152,157))
    #py.legend(loc='lower right')
    #ax2 = ax1.twinx()
    #ax2.plot(utc[::100]-utc[0], temperatures, 'r-', label='Temperature')
    #ax2.set_ylabel('Temperature [K]')
    #ax2.set_xlim((0,200))
    #ax2.set_ylim((0.270,0.325))
    #py.legend(loc='upper right')
    #py.title(test_files[i].split('.')[0]+' - '+str(sigma)+' Sigma Smoothing')
    #py.savefig('azActual_'+test_files[i].split('.')[0]+'_'+str(sigma)+'sigma.png')

    #fig = py.figure(i)
    #fig.clf()
    #ax1 = fig.add_axes((0.1,0.1,0.8,0.8))
    #ax1.plot(utc-utc[0], el_actual, 'k-', label='El Actual')
    #ax1.set_xlabel('Time [s]')
    #ax1.set_ylabel('Elevation [Deg]')
    #ax1.set_xlim((0,200))
    #py.legend(loc='lower right')
    #ax2 = ax1.twinx()
    #ax2.plot(utc[::100]-utc[0], temperatures, 'r-', label='Temperature')
    #ax2.set_ylabel('Temperature [K]')
    #ax2.set_xlim((0,200))
    #ax2.set_ylim((0.270,0.325))
    #py.legend(loc='upper right')
    #py.title(test_files[i].split('.')[0]+' - '+str(sigma)+' Sigma Smoothing')
    #py.savefig('elActual_'+test_files[i].split('.')[0]+'_'+str(sigma)+'sigma.png')
