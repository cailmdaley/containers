import numpy as np
from spt3g import core

dir1 = '/home/sptdaq/spt3g_temp_logs/pole_cooldown2/'
file1 = dir1 + 'log_2017-01-28T02h31.dat'
file2 = dir1 + 'log_2017-01-29T20h10.dat'

f1 = open(file1,'r')
lines1 = f1.readlines()
f1.close()
fridgenames = filter(None,lines1[1].split(' '))
fridgenames.pop(0)
fridgenames.pop(len(fridgenames)-1)
sfdates = []
sftemps = []
for line in lines1:
    if line[0] == '#':
        continue
    line = line.replace('\n','')
    fstrings = filter(None,line.split(' '))
    sfdates.append(fstrings[0:2])
    sftemps.append(fstrings[2:])
fridge = {}
fridge['dates'] = []
for sfdate in sfdates:
    sdate = filter(None,sfdate[0].split('/'))
    fridge['dates'].append(core.G3Time(sdate[0]+'-'+sdate[1]+'-'+sdate[2]+':'+sfdate[1]))
fridge['mjds'] = np.asarray([date.mjd for date in fridge['dates']])
for name in fridgenames:
    fridge[name] = []
for stemp in sftemps:
    for name, temp in zip(fridgenames,stemp):
        fridge[name].append(np.float(temp))
for name in fridgenames:
    fridge[name] = np.asarray(fridge[name])

mjd0 = core.G3Time('01-Feb-2017:00:00:00').mjd




