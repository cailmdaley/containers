import numpy as np
import os
from spt3g import core, gcp, std_processing

def grab_antenna0(f, boardlist = []):
    key = 'antenna0'
    try:
        boardlist.append(f[key])
    except:
        pass

arcdir='/buffer/arc/'
date = '23-Jan-2017'
times = ['12:16:06',
         '12:50:13',
         '13:18:16',
         '13:48:34',
         '14:21:23',
         '14:49:58',
         '15:18:21',
         '15:46:35',
         '16:45:17',
         '17:19:58',
         '17:51:38']

gtimes = [core.G3Time(date + ':' + time) for time in times]
benchpos = []

for gtime in gtimes:
    data1 = []
#    start_time = gtime
#    stop_time = gtime.__add__(1.*core.G3Units.s)
    start_time = gtime.__add__(120.*core.G3Units.s)
    stop_time = gtime.__add__(120.*core.G3Units.s)
    pipe1 = core.G3Pipeline()
    pipe1.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe1.Add(gcp.ARCExtract)
    pipe1.Add(grab_antenna0, boardlist=data1)
    pipe1.Run()
    thisbp = []
    for i in range(6):
        thisbp.append(data1[0]['scu']['benchActual'][i][0])
    benchpos.append(thisbp)

benchfoc = [3.,-3.,5.,7.,9.,11.,13.,15.,17.,19.]
nbenchpos = np.zeros([len(times),6])
for i in range(len(benchpos)):
    nbenchpos[i,:] = np.asarray(benchpos[i])
