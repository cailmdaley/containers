import numpy as np
#from spt3g.examples import convert-sptpol-data as csd
from subprocess import check_call
import sptpol_software.util.time as time


obs = open('scans_ra0hdec-57.5.txt','r').read().split('\n')[5:-1]

start_times = []
stop_times = []

nobs = 20

#Extract start and stop times
for i in range(len(obs)):
    this_obs = filter(None,obs[i].split(' '))

    year = this_obs[1].split('-')[-1].split(':')[0]

    if year == '2015':
        start_times.append(filter(None,obs[i].split(' '))[1])
        stop_times.append(filter(None,obs[i].split(' '))[2])


#Truncate the list for quicklook stuff.
skip = int(len(start_times)/nobs)
starts = start_times[::skip]
stops = stop_times[::skip]

#Now save a G3 file for this observations.
for i in range(len(starts)):
    start = time.SptDatetime(starts[i])
    stop = time.SptDatetime(stops[i])
    start_tag = str(start.year)+str(start.month).zfill(2)+str(start.day).zfill(2)+'_'+str(start.hour).zfill(2)+str(start.minute).zfill(2)+str(start.second).zfill(2)
    stop_tag = str(stop.year)+str(stop.month).zfill(2)+str(stop.day).zfill(2)+'_'+str(stop.hour).zfill(2)+str(stop.minute).zfill(2)+str(stop.second).zfill(2)

    check_call(["python", '/home/jhenning/sptpol_code/spt3g_software/examples/convert-sptpol-data.py', str(starts[i]), str(stops[i]), '/data57/jhenning/spt3g_pointing/data/500d_'+start_tag+'-'+stop_tag+'.g3'])
