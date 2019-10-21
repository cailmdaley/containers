from spt3g import core, dfmux, std_processing, util, gcp
import numpy as np
import glob, pickle
import get_temperatures
from spt3g.util import tctools
from matplotlib.pyplot import *

#code mostly shamelessly stolen from TC

''' 
USAGE:

step 1: edit cycle_start and cycle_stop at the top of the script to reflect the start and stop times you want.  Edit arcdir if necessry.  

step 2: in a terminal type:
python /path/to/code/pole_fridge_plots.py

plots will be saved in the directory from which you ran the script
'''

def grab_reg(f, reg='temperature', data1=[]):
    key = 'array'
    key2 = 'cryo'
    try:
        data1.append(f[key][key2][reg])
    except:
        pass

arcdir = '/spt_data/arc/'

#-54
#cycle_start = core.G3Time('190130 06:07:58')
#cycle_end = core.G3Time('190130 10:39:03')

#-55*
#cycle_start = core.G3Time('190130 10:42:21')
#cycle_end = core.G3Time('190130 15:06:03')

#-56
#cycle_start = core.G3Time('190131 10:44:26')
#cycle_end = core.G3Time('190131 15:28:07')

#-58
#cycle_start = core.G3Time('190131 15:47:21')
#cycle_end = core.G3Time('190131 20:02:27')


#-59*
#cycle_start = core.G3Time('190131 20:08:04')
#cycle_end = core.G3Time('190201 00:35:40')

#-61*
#cycle_start = core.G3Time('190202 08:53:34')
#cycle_end = core.G3Time('190202 13:18:27')

#-62
#cycle_start = core.G3Time('190203 07:01:12')
#cycle_end = core.G3Time('190203 11:34:41')

#-63
#cycle_start = core.G3Time('190203 11:49:54')
#cycle_end = core.G3Time('190203 16:12:34')

#-64
#cycle_start = core.G3Time('190204 12:18:54')
#cycle_end = core.G3Time('190204 17:12:32')

#-65
#cycle_start = core.G3Time('190204 23:31:58')
#cycle_end = core.G3Time('190205 04:43:53')

#-66
#cycle_start = core.G3Time('190205 05:26:33')
#cycle_end = core.G3Time('190205 09:48:10')

#-67
#cycle_start = core.G3Time('190205 09:56:05')
#cycle_end = core.G3Time('190205 14:53:24')

#-69*
#cycle_start = core.G3Time('190206 11:30:05')
#cycle_end = core.G3Time('190206 16:09:26') #fridge blows at 190207 09:11:08


#-73*
#cycle_start = core.G3Time('190207 22:21:35')
#cycle_end = core.G3Time('190208 02:42:06')

#-74*
#cycle_start = core.G3Time('190208 21:52:49')
#cycle_end = core.G3Time('190209 02:25:57')

#-76*
#cycle_start = core.G3Time('190209 19:15:48')
#cycle_end = core.G3Time('190209 23:59:07')


#77*
#cycle_start = core.G3Time('190211 19:10:23')
#cycle_end = core.G3Time('190211 23:41:35') #fridge blows at ~ 190212 15:30:00

#79*
cycle_start = core.G3Time('190212 15:53:24')
cycle_end = core.G3Time('190212 20:17:25')



data1 = []

pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time=cycle_start, stop_time=cycle_end+1.*core.G3Units.hour, basedir=arcdir)
pipe1.Add(gcp.ARCExtract)
pipe1.Add(grab_reg, data1=data1)
pipe1.Run()

ucstage = []
icstage = []
uchead = []
ichead = []
he4head = []
ucpump = []
icpump = []
he4pump = []
head_4k = []
head_50k = []
lctower = []

for data in data1:
    ucstage.append(data[0][10])
    icstage.append(data[0][12])
    uchead.append(data[0][0])
    ichead.append(data[0][1])
    he4head.append(data[0][2])
    ucpump.append(data[0][6])
    icpump.append(data[0][5])
    he4pump.append(data[0][4])
    head_4k.append(data[0][13])
    head_50k.append(data[0][15])
    lctower.append(data[0][11])
    
ucstage = np.asarray(ucstage)
icstage = np.asarray(icstage)
uchead = np.asarray(uchead)
ichead = np.asarray(ichead)
he4head = np.asarray(he4head)
ucpump = np.asarray(ucpump)
icpump = np.asarray(icpump)
he4pump = np.asarray(he4pump)
head_4k = np.asarray(head_4k)
head_50k = np.asarray(head_50k)
lctower = np.asarray(lctower)
mjds = np.arange(len(ucstage))/86400. + cycle_start.mjd

figname = 'Cycle_'+cycle_start.Summary()[0:20]+'.png'
figure(2)
clf()
hrs_after_cycle = (mjds - cycle_end.mjd)*24.
plot(hrs_after_cycle,uchead,label='UC head',color='g')
plot(hrs_after_cycle,ichead,label='IC head',color='r')
plot(hrs_after_cycle,he4head,label='He4 head',color='k')
plot(hrs_after_cycle,ucstage,label='UC stage',color='b')
plot(hrs_after_cycle,lctower,label='LC tower',color='y')
plot(hrs_after_cycle,ucpump,'--',label='UC pump',color='g')
plot(hrs_after_cycle,icpump,'--',label='IC pump',color='r')
plot(hrs_after_cycle,he4pump,'--',label='He4 pump',color='k')
xlim(-6,1)
yscale('log')
ylim(0.25,100)
ylabel('T [K]')
xlabel('hours after drop_bolos finishes')
grid(which='both') 
legend(loc=1)
savefig(figname)

figname2 = 'Cycle_'+cycle_start.Summary()[0:20]+'_lowT.png'
figure(3)
clf()
hrs_after_cycle = (mjds - cycle_end.mjd)*24.
plot(hrs_after_cycle,uchead,label='UC head',color='g')
plot(hrs_after_cycle,ichead,label='IC head',color='r')
plot(hrs_after_cycle,ucstage,label='UC stage',color='b')
plot(hrs_after_cycle,lctower,label='LC tower',color='y')
xlim(-6,1)
yscale('linear')
ylim(0.25,0.4)
ylabel('T [K]')
xlabel('hours after drop_bolos finishes')
grid(which='both') 
legend(loc=1)
savefig(figname2)
