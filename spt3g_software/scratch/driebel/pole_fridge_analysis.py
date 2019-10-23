from spt3g import core, std_processing
import numpy as np
from spt3g.util import extractdata
from matplotlib.pyplot import *
import datetime as dt
import time as Time
import os
import matplotlib.dates as mdates
from read_sched_table import *
import argparse
from dsr_tools import *
from last_eight_cycles import *

program_start_time = Time.time()

#  Winter 2019 DSR

''' 
USAGE:
This file has been heavily modified since the summer.

Run this file from the command line as such:
python pole_fridge_analysis_working.py <optional start time>

If no start time is passed, it will analyze the last completed observation run.
The program pulls data from 30min after the start of the next cycle_tune.  Thus, you must wait until
that data has been written to an archive file, 30-45 minutes after the start of the next cycle_tune.
Give it a good hour just to make sure.

To analyze another cycle of your choice, pass the start time of the cycle_tune.sch exactly as
it appears in the schedule log table, eg:

python pole_fridge_analysis_working.py '190514 02:19:27'

Note the time stamp is a single string.  The program takes a minute or two to run.

All plots will be placed in the directory specified in plot_dir, below
'''

arcdir = '/spt_data/arc/'
spt_fmt = '%y%m%d %H:%M:%S'
g3_fmt =  '%d-%b-%Y:%H:%M:%S' # The date format output by the Summary() method of a G3Time() object
plot_dir = '/home/driebel/code/spt3g_software/scratch/driebel/cycle_plots/'

if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

parser = argparse.ArgumentParser()
parser.add_argument('target', nargs = "?", default = 'latest_complete_cycle')
args = parser.parse_args()

cycle = read_sched_table(args.target)

print('Cycle {}\\\\'.format(cycle.cycle_num))
print('[[Image(cycle{}.png, 400px)]]'.format(cycle.cycle_num))
print('[[Image(cycle{}_lowT.png​, 400px)]]'.format(cycle.cycle_num))
print('[[Image(cycle{}_hold_time_temps.png, 400px)]]'.format(cycle.cycle_num))
print('[[Image(cycle{}_UC_pump_power.png, 400px)]]\\\\'.format(cycle.cycle_num))
print()

#Define the single large accum_dict.  A bit messy in the definition, b/c it's copy/pasted, but it works.
accum_dict = {}

keys = ['ucstage','icstage','uchead','ichead','he4head','ucpump','icpump',
	'he4pump','4k_head','50k_head','lctower', '4k_strap']
base = ['array', 'cryo', 'temperature', 0]
registers = [10, 12, 0, 1, 2, 6, 5, 4, 13, 15, 11, 14]

for i, val in enumerate(keys):
    accum_dict[val] = base + [registers[i]]

accum_dict['ucswitch'] = ['array','cryo','heater_dac',0,2] #UC switch power
accum_dict['icswitch'] = ['array','cryo','heater_dac',0,1] #IC switch power
accum_dict['4k_lyot']  = ['array','cryo','temperature', 1, 15]
accum_dict['50k_wbp'] = ['array','cryo','temperature',1,0]
accum_dict['uc_pump_pwr'] = ['array', 'cryo', 'heater_dac', 0, 5]
accum_dict['spt_time'] = ['array','frame','utc']

archive_start = Time.time()

mult_accumulator = extractdata.MultiAccumulator(accum_dict)
pipe1 = core.G3Pipeline()
pipe1.Add(std_processing.ARCTimerangeReader, start_time = cycle.spt_cycle_start,
          stop_time = cycle.spt_schedule_end + 30.*core.G3Units.minute, basedir=arcdir)
pipe1.Add(mult_accumulator)
pipe1.Run()
data = mult_accumulator.extract_values()
data['time'] = np.asarray([spt_to_dt(i) for i in data['spt_time']])

archive_time = Time.time() - archive_start


# Data now contains all register data @1Hz for entire schedule until 30 minutes after schedule ends.  
# Timestamps must be found via numpy here.

# find when stages blow.  Search from 1 hour before schedule end to 30 min after (end of array)

search_window_start = cycle.dt_schedule_end - dt.timedelta(hours=1)
search_window_end = cycle.dt_schedule_end + dt.timedelta(minutes=30)

start_index = np.where(data['time'] == search_window_start)[0][0]
stop_index = np.where(data['time'] == search_window_end)[0][0]

uc_blown_index = np.where(movingAvg(data['ucstage'][start_index:stop_index],5) >= 0.317)[0][0]
ic_blown_index = np.where(movingAvg(data['icstage'][start_index:stop_index],5) >= 0.367)[0][0]

cycle.uc_blown = data['time'][start_index:stop_index][uc_blown_index].strftime(spt_fmt)
cycle.ic_blown = data['time'][start_index:stop_index][ic_blown_index].strftime(spt_fmt)

# stage blown times done

print('* Cycle {0} - {1} dither {2}'.format(cycle.cycle_num, cycle.cadence, cycle.dither))
print(' * Cycle Time: {}'.format(str(cycle.cycle_time)))
print(' * Hold Time: {}'.format(str(cycle.hold_time)))
print(' * Total Time: {}'.format(str(cycle.total_time)))
print(' * Efficiency: {0:.1f}%'.format(cycle.efficiency*100))

if ((cycle.dt_schedule_end < cycle.dt_uc_blown) and
    (cycle.dt_schedule_end < cycle.dt_ic_blown)):
    print(" * Neither Stage blew before the schedule ended")

if ((cycle.dt_schedule_end > cycle.dt_uc_blown) and
    (cycle.dt_schedule_end < cycle.dt_ic_blown)):
    print(" * The UC Stage blew before the schedule ended")

if ((cycle.dt_schedule_end < cycle.dt_uc_blown) and
    (cycle.dt_schedule_end > cycle.dt_ic_blown)):
    print(" * The IC Stage blew before the schedule ended")
    ic_duration = str(cycle.dt_ic_blown - cycle.dt_cycle_end)
    print(" * IC hold time: {}".format(ic_duration))

if ((cycle.dt_schedule_end > cycle.dt_uc_blown) and
    (cycle.dt_schedule_end > cycle.dt_ic_blown)):
    print(" * Both Stages blew before the schedule ended")
    ic_duration = str(cycle.dt_ic_blown - cycle.dt_cycle_end)
    print(" * IC hold time: {}".format(ic_duration))

print(" * 4K Head Temp at start of cycle: {0:.2f} K".format(np.average(data['4k_head'][0:5])))
print(" * 4K Strap temp at start of cycle: {0:.2f} K".format(np.average(data['4k_strap'][0:5])))

print()


# The switches activate beteen -3 and -2 hours after the end of the cycle.
# Focus on this "window" of time, and find when the switch turns on (power > 0).
# pull out -3 to -2 hour window
window_start = cycle.dt_cycle_end - dt.timedelta(hours = 3)
window_end = cycle.dt_cycle_end - dt.timedelta(hours = 1.5)
window = [np.asarray(data['time'] > window_start) & np.asarray(data['time'] <= window_end)][0]


right_before_ucswitch = np.where(data['ucswitch'][window] > 0)[0][0] - 1
print('Right before the UC switch is activated:')
print('UC Head temp = {0:.2f} K'.format(data['uchead'][window][right_before_ucswitch]))
print('UC Stage temp = {0:.2f} K'.format(data['ucstage'][window][right_before_ucswitch]))
print('LC Tower temp = {0:.2f} K'.format(np.average(data['lctower'][window][right_before_ucswitch-2:right_before_ucswitch+3])))
print()
print('5 minutes before that:')
UC5 = right_before_ucswitch - (5 * 60)
print('UC Head temp = {0:.2f} K'.format(data['uchead'][window][UC5]))
print('UC Stage temp = {0:.2f} K'.format(data['ucstage'][window][UC5]))
average_start = right_before_ucswitch - (5 * 60) - 2
average_end = right_before_ucswitch - (5 * 60) + 3
print('LC Tower temp = {0:.2f} K'.format(np.average(data['lctower'][window][average_start:average_end])))
print()
right_before_icswitch = np.where(data['icswitch'][window] > 0)[0][0] - 1
print('Right before the IC switch is activated:')
print('IC Head temp = {0:.2f} K'.format(data['ichead'][window][right_before_icswitch]))
print('IC Stage temp = {0:.2f} K'.format(data['icstage'][window][right_before_icswitch]))
print()
print('5 minutes before that:')
IC5 = right_before_icswitch - (5 * 60)
print('IC Head temp = {0:.2f} K'.format(data['ichead'][window][IC5]))
print('IC Stage temp = {0:.2f} K'.format(data['icstage'][window][IC5]))

mjds = np.arange(len(data['ucstage']))/86400. + cycle.spt_cycle_start.mjd
hrs_after_cycle = (mjds - cycle.spt_cycle_end.mjd)*24.

#Define first plot window: from cycle_start to 1 hour after cycle_end
window_start = cycle.dt_cycle_start
window_end = cycle.dt_cycle_end + dt.timedelta(hours = 1)
plot_window = [np.asarray(data['time'] >= window_start) & np.asarray(data['time'] <= window_end)][0]



figname1 = plot_dir+'cycle{}.png'.format(cycle.cycle_num)
figure(1)
clf()
plot(hrs_after_cycle[plot_window],data['uchead'][plot_window],label='UC head',color='g')
plot(hrs_after_cycle[plot_window],data['ichead'][plot_window],label='IC head',color='r')
plot(hrs_after_cycle[plot_window],data['he4head'][plot_window],label='He4 head',color='k')
plot(hrs_after_cycle[plot_window],data['ucstage'][plot_window],label='UC stage',color='b')
plot(hrs_after_cycle[plot_window],data['lctower'][plot_window],label='LC tower',color='y')
plot(hrs_after_cycle[plot_window],data['ucpump'][plot_window],'--',label='UC pump',color='g')
plot(hrs_after_cycle[plot_window],data['icpump'][plot_window],'--',label='IC pump',color='r')
plot(hrs_after_cycle[plot_window],data['he4pump'][plot_window],'--',label='He4 pump',color='k')
xlim(-6,1)
yscale('log')
ylim(0.25,100)
ylabel('T [K]')
xlabel('hours after drop_bolos finishes')
grid(which='both') 
legend(loc=1)
savefig(figname1)
close(1)


figname2 = plot_dir+'cycle{}_lowT.png'.format(cycle.cycle_num)
figure(2)
clf()
plot(hrs_after_cycle[plot_window],data['uchead'][plot_window],label='UC head',color='g')
plot(hrs_after_cycle[plot_window],data['ichead'][plot_window],label='IC head',color='r')
plot(hrs_after_cycle[plot_window],data['ucstage'][plot_window],label='UC stage',color='b')
plot(hrs_after_cycle[plot_window],data['lctower'][plot_window],label='LC tower',color='y')
xlim(-6,1)
yscale('linear')
ylim(0.25,0.4)
ylabel('T [K]')
xlabel('hours after drop_bolos finishes')
grid(which='both') 
legend(loc=1)
savefig(figname2)
close(2)


#Define second plot window, from cycle_end to 6 hours after cycle_end
window_start = cycle.dt_cycle_end
window_end = cycle.dt_cycle_end + dt.timedelta(hours = 6)
plot_window = [np.asarray(data['time'] >= window_start) & np.asarray(data['time'] <= window_end)][0]

timestamp = np.arange(len(data['50k_head'][plot_window]))/60 # minutes


figname3 = plot_dir + 'Cycle{}_50K.png'.format(cycle.cycle_num)
figure(3)
plot(timestamp[::60], data['50k_head'][plot_window][::60], 'b.', markersize=2, 
     label='PT 50K Head')
plot(timestamp[::60], data['50k_wbp'][plot_window][::60], 'g.', markersize=2, 
     label='B1 50K WBP Near')
ylabel('Temperature [K]')
xlabel('Time Since Cycle End [min]')
title('Cycle {}: 50K'.format(cycle.cycle_num))
legend(loc=0)
savefig(figname3)
close(3)

figname4 = plot_dir + 'Cycle{}_4K.png'.format(cycle.cycle_num)
figure(4)
plot(timestamp[::60], data['4k_strap'][plot_window][::60], 'k.', markersize=2, 
     label='PT 4K Strap')
plot(timestamp[::60], data['4k_lyot'][plot_window][::60], 'r.', markersize=2, 
     label='R4 4K Lyot')
ylabel('Temperature [K]')
xlabel('Time Since Cycle End [min]')
title('Cycle {}: 4K'.format(cycle.cycle_num))
legend(loc=0)
savefig(figname4)
close(4)


# Define the thrid plot window, from 5 minutes before the cycle end to 5 minutes after UC blows.
window_start = cycle.dt_cycle_end - dt.timedelta(minutes = 5)
window_end = cycle.dt_uc_blown + dt.timedelta(minutes = 5)
plot_window = [np.asarray(data['time'] >= window_start) & np.asarray(data['time'] <= window_end)][0]


#xtick_labels=[spt_to_dt(i) for i in data['time'][plot_window]]
xtick_labels = data['time'][plot_window]

figname5 = plot_dir + 'cycle{}_hold_time_temps.png'.format(cycle.cycle_num)
figure(5)
fig, ax = subplots(1)
ax.plot(xtick_labels[::60], data['uchead'][plot_window][::60],'g',label='UC Head')
ax.plot(xtick_labels[::60], data['ucstage'][plot_window][::60],'b',label='UC Stage')
ax.plot(xtick_labels[::60], data['ichead'][plot_window][::60],'r',label='IC Head')
ax.plot(xtick_labels[::60], data['lctower'][plot_window][::60],'y',label='LC Tower')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
yscale('linear')
ylim(0.28,0.35)
ylabel('T [K]')
xlabel('Time')
#fig.autofmt_xdate() 
#xticks(rotation=30)
grid(which='both') 
legend(loc=0)
savefig(figname5)
close(5)

figname6 = plot_dir + 'cycle{}_UC_pump_power.png'.format(cycle.cycle_num)
figure(6)
fig2, ax2 = subplots(1)
ax2.plot(xtick_labels[::60], data['uc_pump_pwr'][plot_window][::60],'k',label='UC Pump Power')
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ylabel('Power [mW]')
xlabel('Time')
grid(which='both') 
ylim(0,60)
legend(loc=0)
savefig(figname6)
close('all')

print()
print('The last eight completed cycles were: {}'.format(last_eight_cycles()))
print()
program_length = Time.time() - program_start_time
print('Analysis took: {0:.2f} seconds'.format(program_length))
print('Accessing the archive files took: {0:.2f} sec'.format(archive_time))
print('Archive access is {:.2f}% of runtime'.format(archive_time/program_length * 100))
